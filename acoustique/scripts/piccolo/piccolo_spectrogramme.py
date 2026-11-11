"""Spectrogramme du Piccolo II: acquisition en direct et visualisation de .wav.

Deux modes dans la même fenêtre:

- Direct: « Démarrer » lance le flux audio du Piccolo II et affiche un
  spectrogramme qui défile (chaque colonne est la FFT d'une seule trame avec
  fenêtre de Hann, sans moyenne). L'audio est écrit sur le disque au fur et à
  mesure dans un .wav temporaire, donc la durée n'est limitée que par l'espace
  disque. « Arrêter » affiche tout l'enregistrement, et « Enregistrer .wav »
  le déplace à l'endroit voulu.
- Fichier: « Ouvrir » affiche le spectrogramme complet d'un .wav. Le fichier
  est lu directement sur le disque, seule la zone visible est calculée. La
  molette zoome, le glisser déplace la vue. Le curseur bleu (à glisser) donne
  le spectre instantané à cet instant dans le panneau du bas.

Exemples:
    python piccolo_spectrogramme.py
    python piccolo_spectrogramme.py ../../data/mesure.wav
    python piccolo_spectrogramme.py --nfft 8192 --fmax 5000 --fenetre 60

Fermer Piccolo Utility (PIC2_*.exe) avant une acquisition en direct.

En direct, avec une fenêtre longue, le spectrogramme affiche au plus
COLS_DIRECT colonnes: une trame sur plusieurs est montrée, mais toutes les
données sont enregistrées.

Niveaux en dBFS (0 dB = sinus pleine échelle) plus --cal. Pour avoir des
dB SPL, ajuster --cal contre la lecture LZeq de l'appareil (dépend de la
plage Low/High).
"""

import argparse
import gc
import math
import os
import queue
import shutil
import struct
import sys
import tempfile
import threading
import wave
from datetime import datetime
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from numpy.lib.stride_tricks import sliding_window_view
from pyqtgraph.Qt import QtCore, QtWidgets
from scipy import fft as sfft
from scipy.io import wavfile

from piccolo_usb import FS, PiccoloII, PiccoloStream, piccolo_present

MSG_ABSENT = ("Le Piccolo II n'est pas détecté.\n\n"
              "Vérifier qu'il est allumé et branché au PC par USB, et que "
              "Piccolo Utility (PIC2_*.exe) est fermé, puis cliquer de nouveau sur Démarrer.")

DATA_DIR = Path(__file__).resolve().parents[2] / 'data'
NFFT_CHOIX = [512, 1024, 2048, 4096, 8192, 16384]
FENETRES = [10, 30, 60, 120, 300, 600]   # durée visible en direct [s]
COLS_DIRECT = 2500     # colonnes max du spectrogramme en direct
MAX_COLS = 1500        # colonnes calculées pour la zone visible d'un fichier
F_MIN_PIC = 20.0       # le pic affiché ignore le bruit sous cette fréquence [Hz]
ENV_BLOC = 1024        # échantillons par point de l'enveloppe de la forme d'onde
MAX_PTS_BRUTS = 200_000


# --------------------------------------------------------------------------
# Calcul des spectres (signal en float, pleine échelle = 1)
# --------------------------------------------------------------------------
def spectre_db(trames, fenetre, cal):
    # la moyenne de chaque trame est retirée: le Piccolo a un petit décalage DC
    # qui, sinon, déborde dans les premiers canaux de fréquence
    trames = np.asarray(trames, np.float32)
    trames = trames - trames.mean(axis=-1, keepdims=True)
    spec = sfft.rfft(trames * fenetre.astype(np.float32), axis=-1, workers=-1)
    mag = np.abs(spec) * (2.0 / fenetre.sum())
    return 20 * np.log10(mag + 1e-12) + cal


class StftContinue:
    """STFT d'un flux: chaque appel à push() retourne les nouvelles trames."""

    def __init__(self, nfft, hop, cal):
        self.nfft, self.hop, self.cal = nfft, hop, cal
        self.fen = np.hanning(nfft)
        self.buf = np.zeros(0, np.float32)

    def push(self, x):
        self.buf = np.concatenate([self.buf, x])
        n = (len(self.buf) - self.nfft) // self.hop + 1
        if n <= 0:
            return None
        idx = np.arange(self.nfft)[None, :] + self.hop * np.arange(n)[:, None]
        out = spectre_db(self.buf[idx], self.fen, self.cal)
        self.buf = self.buf[n * self.hop:]
        return out


def niveau_db(x, cal):
    rms = np.sqrt(np.mean((x - x.mean()) ** 2))
    return 20 * np.log10(rms + 1e-12) + cal + 3.01   # sinus pleine échelle = 0 dB


def pic(freqs, db):
    k = np.searchsorted(freqs, F_MIN_PIC)
    i = k + int(np.argmax(db[k:]))
    return freqs[i], db[i]


def fmt_duree(s):
    if s < 60:
        return f"{s:.2f} s"
    h, reste = divmod(s, 3600)
    m, s = divmod(reste, 60)
    return (f"{int(h)} h " if h else "") + f"{int(m)} min {s:04.1f} s"


# --------------------------------------------------------------------------
# Signal lu sur le disque
# --------------------------------------------------------------------------
class Signal:
    """Signal mono. raw peut être un memmap int16 plus grand que la mémoire."""

    def __init__(self, raw, fs, scale=1.0):
        self.raw, self.fs, self.scale = raw, fs, scale
        self.n = len(raw)
        self._env = None

    @property
    def duree(self):
        return self.n / self.fs

    def trames(self, centres, nfft):
        """Trames de nfft échantillons centrées sur `centres` (zéros hors du signal)."""
        idx = np.asarray(centres)[:, None] + (np.arange(nfft) - nfft // 2)
        ok = (idx >= 0) & (idx < self.n)
        v = np.asarray(self.raw[np.clip(idx, 0, self.n - 1)], np.float32) * self.scale
        v[~ok] = 0
        return v

    def seg(self, a, b):
        """Échantillons a..b en float32, complétés par des zéros hors du signal."""
        out = np.zeros(b - a, np.float32)
        lo, hi = max(a, 0), min(b, self.n)
        if hi > lo:
            out[lo - a:hi - a] = np.asarray(self.raw[lo:hi], np.float32) * self.scale
        return out

    def stft_zone(self, nfft, i0, i1, cal):
        """Au plus MAX_COLS trames centrées entre i0 et i1. Retourne (dB, t0, dt)."""
        # pas minimal nfft/8: des trames plus serrées n'ajoutent pas de détail
        hop = max(1, nfft // 8, math.ceil((i1 - i0) / MAX_COLS))
        centres = np.arange(i0, i1, hop)
        fen = np.hanning(nfft)
        out = np.empty((len(centres), nfft // 2 + 1), np.float32)
        if i1 - i0 + nfft <= 8_000_000:
            # zone courte: un seul segment contigu lu sur le disque
            seg = self.seg(i0 - nfft // 2, i1 + nfft - nfft // 2)
            vues = sliding_window_view(seg, nfft)
            idx = centres - i0
            for k in range(0, len(idx), 256):
                out[k:k + 256] = spectre_db(vues[idx[k:k + 256]], fen, cal)
        else:
            # zone longue: seulement les trames nécessaires
            for k in range(0, len(centres), 256):
                out[k:k + 256] = spectre_db(self.trames(centres[k:k + 256], nfft), fen, cal)
        return out, centres[0] / self.fs, hop / self.fs

    def enveloppe(self):
        """Min et max par bloc de ENV_BLOC échantillons, calculés une seule fois."""
        if self._env is None:
            nb = math.ceil(self.n / ENV_BLOC)
            mn, mx = np.empty(nb, np.float32), np.empty(nb, np.float32)
            pas = ENV_BLOC * 1024
            for a in range(0, self.n, pas):
                blk = np.asarray(self.raw[a:a + pas], np.float32) * self.scale
                m = math.ceil(len(blk) / ENV_BLOC)
                blk = np.pad(blk, (0, m * ENV_BLOC - len(blk)), mode='edge').reshape(m, -1)
                k = a // ENV_BLOC
                mn[k:k + m], mx[k:k + m] = blk.min(1), blk.max(1)
            self._env = (mn, mx)
        return self._env

    def forme_onde(self, i0, i1, max_pts=4000):
        """(t, y) à tracer pour la zone i0..i1: brut si court, sinon enveloppe."""
        if i1 - i0 <= MAX_PTS_BRUTS:
            y = np.asarray(self.raw[i0:i1], np.float32) * self.scale
            return np.arange(i0, i1) / self.fs, y
        mn, mx = self.enveloppe()
        k0, k1 = i0 // ENV_BLOC, min(len(mn), math.ceil(i1 / ENV_BLOC))
        g = max(1, math.ceil((k1 - k0) / max_pts))
        a, b = mn[k0:k1], mx[k0:k1]
        pad = (-len(a)) % g
        a = np.pad(a, (0, pad), mode='edge').reshape(-1, g).min(1)
        b = np.pad(b, (0, pad), mode='edge').reshape(-1, g).max(1)
        t = (k0 + (np.arange(len(a)) + 0.5) * g) * ENV_BLOC / self.fs
        y = np.empty(2 * len(a), np.float32)
        y[0::2], y[1::2] = a, b
        return np.repeat(t, 2), y


def info_wav(chemin):
    """Parcourt les blocs RIFF. Retourne (fmt, offset des données, nb d'octets)."""
    with open(chemin, 'rb') as f:
        riff, _, wav = struct.unpack('<4sI4s', f.read(12))
        if riff != b'RIFF' or wav != b'WAVE':
            raise ValueError("ce n'est pas un fichier WAV")
        fmt = None
        while True:
            h = f.read(8)
            if len(h) < 8:
                raise ValueError("bloc de données introuvable")
            cid, taille = struct.unpack('<4sI', h)
            if cid == b'fmt ':
                fmt = struct.unpack_from('<HHIIHH', f.read(taille))
                f.seek(taille & 1, 1)
            elif cid == b'data':
                off = f.tell()
                reste = os.path.getsize(chemin) - off
                if taille == 0 or taille > reste:   # fichier interrompu
                    taille = reste
                return fmt, off, taille
            else:
                f.seek(taille + (taille & 1), 1)


def lire_wav_complet(chemin):
    """Formats autres que PCM 16 bits: lecture complète en mémoire."""
    try:
        with wave.open(str(chemin)) as w:
            fs, nc, sw = w.getframerate(), w.getnchannels(), w.getsampwidth()
            raw = w.readframes(w.getnframes())
        if sw == 1:
            x = (np.frombuffer(raw, np.uint8).astype(np.float32) - 128) / 128
        elif sw == 3:
            b = np.frombuffer(raw, np.uint8).reshape(-1, 3)
            v = (b[:, 0].astype(np.int32) | (b[:, 1].astype(np.int32) << 8)
                 | (b[:, 2].astype(np.int32) << 16))
            x = np.where(v >= 1 << 23, v - (1 << 24), v) / 2.0 ** 23
        elif sw == 4:
            x = np.frombuffer(raw, '<i4') / 2.0 ** 31
        else:
            x = np.frombuffer(raw, '<i2') / 32768.0
        return fs, np.asarray(x.reshape(-1, nc)[:, 0], np.float32)
    except wave.Error:   # virgule flottante
        fs, x = wavfile.read(chemin)
        if x.ndim > 1:
            x = x[:, 0]
        if x.dtype == np.int32:
            x = x / 2.0 ** 31
        return fs, np.asarray(x, np.float32)


def ouvrir_signal(chemin):
    """PCM 16 bits: lu sur le disque sans être chargé (memmap). Sinon en mémoire.

    Le module wave n'est pas utilisé pour l'en-tête: il refuse celui des
    fichiers de l'Audio Recorder du Piccolo (débit d'octets incohérent).
    """
    fmt, off, taille = info_wav(chemin)
    if fmt and fmt[0] in (1, 0xFFFE) and fmt[5] == 16:
        nc, fs = fmt[1], fmt[2]
        n = taille // (2 * nc)
        if n == 0:
            raise ValueError("le fichier ne contient aucun échantillon")
        raw = np.memmap(chemin, '<i2', 'r', off, shape=(n * nc,))
        if nc > 1:
            raw = raw.reshape(n, nc)[:, 0]
        return Signal(raw, fs, 1 / 32768)
    fs, x = lire_wav_complet(chemin)
    return Signal(x, fs, 1.0)


# --------------------------------------------------------------------------
# Acquisition en arrière-plan: écrit le .wav et alimente l'affichage
# --------------------------------------------------------------------------
class PiccoloSource(threading.Thread):
    def __init__(self, q, chemin):
        super().__init__(daemon=True)
        self.q, self.chemin = q, chemin
        self.stop_evt = threading.Event()
        self.error = None
        self.version = ''
        self.debordements = 0
        self.n = 0

    def run(self):
        pic_ = w = None
        try:
            pic_ = PiccoloII()
            self.version = pic_.version
            # wave remet l'en-tête à jour à chaque écriture: le fichier reste
            # lisible même si le programme est interrompu
            w = wave.open(str(self.chemin), 'wb')
            w.setnchannels(1)
            w.setsampwidth(2)
            w.setframerate(FS)
            with PiccoloStream(pic_) as s:
                while not self.stop_evt.is_set():
                    x = s.read_chunk()
                    w.writeframes(x.tobytes())
                    self.n += len(x)
                    self.debordements = s.overflows
                    self.q.put(x)
        except Exception as e:
            self.error = e
        finally:
            if w:
                w.close()
            if pic_:
                pic_.close()

    def stop(self):
        self.stop_evt.set()
        self.join(timeout=3)


# --------------------------------------------------------------------------
# Interface
# --------------------------------------------------------------------------
class Fenetre(QtWidgets.QMainWindow):
    def __init__(self, args):
        super().__init__()
        self.args = args
        self.cal = args.cal
        self.nfft = args.nfft
        self.fenetre = args.fenetre
        self.fs = FS
        self.mode = 'vide'          # 'vide', 'direct' ou 'fichier'
        self.src = None
        self.q = queue.Queue()
        self.sig = None             # Signal affiché en mode fichier
        self.fichier_temp = None    # .wav de la dernière acquisition
        self.non_sauve = False
        self.nom = ''

        pg.setConfigOptions(imageAxisOrder='row-major')
        self.setWindowTitle("Piccolo II")
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        vl = QtWidgets.QVBoxLayout(central)

        # Barre de commandes
        barre = QtWidgets.QHBoxLayout()
        vl.addLayout(barre)
        self.b_acq = QtWidgets.QPushButton("Démarrer")
        self.b_acq.setMinimumWidth(110)
        self.b_acq.clicked.connect(self.basculer_acq)
        self.b_ouvrir = QtWidgets.QPushButton("Ouvrir .wav...")
        self.b_ouvrir.clicked.connect(self.ouvrir)
        self.b_sauver = QtWidgets.QPushButton("Enregistrer .wav...")
        self.b_sauver.clicked.connect(self.sauver)
        self.b_tout = QtWidgets.QPushButton("Tout afficher")
        self.b_tout.clicked.connect(self.tout_afficher)
        self.c_nfft = QtWidgets.QComboBox()
        self.c_nfft.addItems([str(n) for n in NFFT_CHOIX])
        self.c_nfft.setCurrentText(str(self.nfft))
        self.c_nfft.currentTextChanged.connect(self.changer_nfft)
        self.c_fen = QtWidgets.QComboBox()
        self.c_fen.addItems([f"{f} s" for f in FENETRES])
        self.c_fen.setCurrentText(f"{self.fenetre} s")
        self.c_fen.currentTextChanged.connect(self.changer_fenetre)
        for w in (self.b_acq, self.b_ouvrir, self.b_sauver, self.b_tout):
            barre.addWidget(w)
        barre.addSpacing(20)
        barre.addWidget(QtWidgets.QLabel("nfft"))
        barre.addWidget(self.c_nfft)
        self.l_res = QtWidgets.QLabel()
        barre.addWidget(self.l_res)
        barre.addSpacing(20)
        barre.addWidget(QtWidgets.QLabel("Fenêtre en direct"))
        barre.addWidget(self.c_fen)
        barre.addStretch()

        gw = pg.GraphicsLayoutWidget()
        vl.addWidget(gw)
        self.l_info = QtWidgets.QLabel()
        self.l_info.setTextFormat(QtCore.Qt.TextFormat.RichText)
        vl.addWidget(self.l_info)

        # Forme d'onde
        self.p_wave = gw.addPlot(row=0, col=0)
        self.p_wave.setLabel('left', 'Amplitude', units='FS')
        self.p_wave.setLabel('bottom', 'Temps', units='s')
        self.p_wave.setMaximumHeight(170)
        self.c_wave = self.p_wave.plot(pen=pg.mkPen('#4fc3f7', width=1))
        self.curseur_wave = pg.InfiniteLine(angle=90, pen=pg.mkPen('#ffb74d'))
        self.p_wave.addItem(self.curseur_wave)
        self.curseur_wave.hide()

        # Spectrogramme
        self.p_spec = gw.addPlot(row=1, col=0)
        self.p_spec.setLabel('left', 'Fréquence', units='Hz')
        self.p_spec.setLabel('bottom', 'Temps', units='s')
        self.img = pg.ImageItem()
        self.p_spec.addItem(self.img)
        cmap = pg.colormap.get('inferno')
        self.img.setColorMap(cmap)
        self.img.setImage(np.full((2, 2), args.min_db, np.float32), autoLevels=False)
        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img)
        self.hist.gradient.setColorMap(cmap)
        self.hist.setLevels(args.min_db, args.max_db)
        self.hist.setHistogramRange(args.min_db - 20, args.max_db + 20)
        gw.addItem(self.hist, row=1, col=1)
        self.curseur = pg.InfiniteLine(angle=90, movable=True,
                                       pen=pg.mkPen('#4fc3f7', width=2),
                                       hoverPen=pg.mkPen('#ffffff', width=3))
        self.curseur.setZValue(10)
        self.p_spec.addItem(self.curseur)
        self.curseur.hide()
        self.curseur.sigPositionChanged.connect(self.maj_curseur)

        # Spectre instantané
        self.p_inst = gw.addPlot(row=2, col=0)
        self.p_inst.setLabel('left', 'Niveau', units='dB')
        self.p_inst.setLabel('bottom', 'Fréquence', units='Hz')
        self.p_inst.setYRange(args.min_db, args.max_db + 10)
        self.p_inst.showGrid(x=True, y=True, alpha=0.3)
        self.p_inst.setMaximumHeight(230)
        self.c_inst = self.p_inst.plot(pen=pg.mkPen('#ffb74d', width=1))

        # Recalcul de la zone visible (mode fichier), avec un léger délai
        self.vue_timer = QtCore.QTimer(singleShot=True)
        self.vue_timer.timeout.connect(self.rendre_vue)
        self.p_spec.sigXRangeChanged.connect(lambda *_: self.vue_timer.start(120))

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.maj_direct)

        self.set_boutons(False)
        self.maj_resolution()
        self.info("Cliquer sur <b>Démarrer</b> pour une acquisition, "
                  "ou sur <b>Ouvrir</b> pour voir un fichier.")

    # -- utilitaires -------------------------------------------------------
    @property
    def fmax(self):
        return min(self.args.fmax or self.fs / 2, self.fs / 2)

    @property
    def freqs(self):
        return np.fft.rfftfreq(self.nfft, 1 / self.fs)

    def info(self, txt):
        self.l_info.setText(txt.replace('  ', '&nbsp;&nbsp;'))

    def erreur(self, txt):
        self.info(f"<span style='color:#c62828'>{txt}</span>")

    def maj_resolution(self):
        self.l_res.setText(f"  df = {self.fs / self.nfft:.1f} Hz, "
                           f"trame = {1000 * self.nfft / self.fs:.0f} ms")

    def confirmer_perte(self):
        if not self.non_sauve:
            return True
        r = QtWidgets.QMessageBox.question(
            self, "Acquisition non enregistrée",
            "La dernière acquisition n'est pas enregistrée. Continuer quand même?")
        return r == QtWidgets.QMessageBox.StandardButton.Yes

    def set_boutons(self, direct):
        self.b_acq.setText("Arrêter" if direct else "Démarrer")
        self.b_ouvrir.setEnabled(not direct)
        self.c_nfft.setEnabled(not direct)
        self.b_tout.setEnabled(not direct and self.mode == 'fichier')
        self.b_sauver.setEnabled(not direct and self.non_sauve)

    def liberer_temp(self):
        """Supprime le .wav temporaire de la dernière acquisition."""
        if self.fichier_temp is None:
            return
        self.sig = None      # ferme le memmap avant la suppression
        gc.collect()
        try:
            os.remove(self.fichier_temp)
        except OSError:
            pass
        self.fichier_temp = None
        self.non_sauve = False

    # -- mode direct -------------------------------------------------------
    def basculer_acq(self):
        if self.mode == 'direct':
            self.arreter()
        else:
            self.demarrer()

    def alerte(self, titre, txt):
        self.erreur(txt.replace('\n', ' '))
        QtWidgets.QMessageBox.warning(self, titre, txt)

    def demarrer(self):
        if not piccolo_present():
            self.alerte("Piccolo II introuvable", MSG_ABSENT)
            return
        if not self.confirmer_perte():
            return
        self.liberer_temp()
        self.sig = None
        self.fichier_temp = Path(tempfile.gettempdir()) / \
            f"piccolo_{datetime.now():%Y%m%d_%H%M%S}.wav"
        self.mode = 'direct'
        self.fs = FS
        self.q = queue.Queue()
        self.curseur.hide()
        self.curseur_wave.hide()
        self.p_wave.setXLink(None)
        self.wave = np.zeros(int(0.1 * self.fs), np.float32)
        self.t_wave = np.arange(len(self.wave)) / self.fs
        self.c_wave.setData(self.t_wave, self.wave)
        self.p_wave.setXRange(0, 0.1, padding=0)
        self.p_wave.enableAutoRange(axis='y')
        self.p_inst.setXRange(0, self.fmax)
        self.maj_resolution()
        self.configurer_direct()
        self.src = PiccoloSource(self.q, self.fichier_temp)
        self.src.start()
        self.set_boutons(True)
        self.info("Connexion au Piccolo II...")
        self.timer.start(30)

    def configurer_direct(self):
        """(Re)crée l'image défilante pour la fenêtre choisie."""
        fen = self.fenetre
        self.hop = max(self.nfft // 8, math.ceil(fen * self.fs / COLS_DIRECT))
        self.stft = StftContinue(self.nfft, self.hop, self.cal)
        self.ncols = int(fen * self.fs / self.hop)
        nb = self.nfft // 2 + 1
        df = self.fs / self.nfft
        self.img_data = np.full((nb, self.ncols), self.args.min_db, np.float32)
        self.img.setImage(self.img_data, autoLevels=False)
        self.img.setRect(QtCore.QRectF(-fen, -df / 2, fen, nb * df))
        self.p_spec.setLimits(xMin=-fen, xMax=0, yMin=0, yMax=self.fs / 2, minXRange=None)
        self.p_spec.setXRange(-fen, 0, padding=0)
        self.p_spec.setYRange(0, self.fmax, padding=0)
        self.p_spec.setMouseEnabled(x=False, y=True)

    def vider_file(self):
        blocs = []
        while True:
            try:
                blocs.append(self.q.get_nowait())
            except queue.Empty:
                return blocs

    def maj_direct(self):
        if self.src.error:
            err, n = self.src.error, self.src.n
            self.arreter()
            if n:
                txt = (f"L'acquisition s'est interrompue: {err}\n\n"
                       f"Les {fmt_duree(n / FS)} enregistrées jusque-là sont conservées "
                       "et affichées. Cliquer sur « Enregistrer .wav » pour les sauver.")
            elif not piccolo_present():
                txt = MSG_ABSENT
            else:
                txt = f"Impossible de démarrer l'acquisition: {err}"
            self.alerte("Acquisition interrompue", txt)
            return
        blocs = self.vider_file()
        if not blocs:
            return
        x = np.concatenate(blocs).astype(np.float32) / 32768
        cols = self.stft.push(x)
        if cols is not None:
            k = min(len(cols), self.ncols)
            self.img_data = np.roll(self.img_data, -k, axis=1)
            self.img_data[:, -k:] = cols[-k:].T
            self.img.setImage(self.img_data, autoLevels=False)
            self.c_inst.setData(self.freqs, cols[-1])
            f, lp = pic(self.freqs, cols[-1])
        else:
            f, lp = 0, float('nan')
        self.wave = np.concatenate([self.wave, x])[-len(self.wave):]
        self.c_wave.setData(self.t_wave, self.wave)
        n = self.src.n
        crete = np.abs(x).max()
        sat = ("  <span style='color:#c62828'><b>SATURATION</b></span>"
               if crete > 0.99 else "")
        pertes = (f"  <span style='color:#c62828'>débordements: {self.src.debordements}</span>"
                  if self.src.debordements else "")
        self.info(f"<b>Acquisition</b>: {fmt_duree(n / self.fs)} ({2 * n / 1e6:.1f} Mo)  |  "
                  f"niveau {niveau_db(x, self.cal):6.1f} dB  |  crête {100 * crete:5.1f} % FS  |  "
                  f"pic {f:7.1f} Hz à {lp:6.1f} dB{sat}{pertes}")

    def arreter(self):
        self.timer.stop()
        if self.src:
            self.src.stop()
            n = self.src.n
            self.src = None
        else:
            n = 0
        self.vider_file()
        if n:
            self.non_sauve = True
            self.charger(self.fichier_temp,
                         nom=f"Acquisition du {datetime.now():%Y-%m-%d %H:%M:%S}")
        else:
            self.liberer_temp()
            self.mode = 'vide'
            self.set_boutons(False)
            self.info("Acquisition arrêtée, aucune donnée reçue.")

    def changer_fenetre(self, txt):
        self.fenetre = int(txt.split()[0])
        if self.mode == 'direct':
            self.configurer_direct()

    # -- mode fichier ------------------------------------------------------
    def ouvrir(self):
        if not self.confirmer_perte():
            return
        chemin, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Ouvrir un enregistrement", str(DATA_DIR), "Fichiers WAV (*.wav)")
        if chemin:
            self.liberer_temp()
            self.charger(chemin)

    def sauver(self):
        defaut = DATA_DIR / f"{self.fichier_temp.stem}.wav"
        dest, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Enregistrer l'acquisition", str(defaut), "Fichiers WAV (*.wav)")
        if not dest:
            return
        vue = (self.p_spec.viewRange(), self.curseur.value())
        self.sig = None      # ferme le memmap pour pouvoir déplacer le fichier
        gc.collect()
        try:
            shutil.move(str(self.fichier_temp), dest)
        except OSError as e:
            self.charger(self.fichier_temp, nom=self.nom, vue=vue)
            self.erreur(f"Enregistrement impossible: {e}")
            return
        self.fichier_temp = None
        self.non_sauve = False
        self.charger(dest, vue=vue)

    def charger(self, chemin, nom=None, vue=None):
        try:
            sig = ouvrir_signal(chemin)
        except Exception as e:
            self.erreur(f"Lecture impossible: {e}")
            return
        self.sig, self.fs = sig, sig.fs
        self.nom = nom or Path(chemin).name
        self.mode = 'fichier'
        self.set_boutons(False)
        self.maj_resolution()
        self.p_wave.setXLink(self.p_spec)
        self.p_wave.enableAutoRange(axis='y')
        self.p_spec.setMouseEnabled(x=True, y=True)
        self.p_spec.setLimits(xMin=0, xMax=sig.duree, yMin=0, yMax=sig.fs / 2,
                              minXRange=min(0.01, sig.duree))
        self.p_inst.setXRange(0, self.fmax)
        self.curseur.setBounds([0, sig.duree])
        self.curseur.show()
        self.curseur_wave.show()
        self.setWindowTitle(f"Piccolo II  |  {self.nom}")
        if vue:
            (xr, yr), t = vue
            self.p_spec.setXRange(*xr, padding=0)
            self.p_spec.setYRange(*yr, padding=0)
            self.rendre_vue()
            self.curseur.setValue(t)
        else:
            self.tout_afficher()
            self.curseur.setValue(sig.duree / 2)
        self.maj_curseur()

    def tout_afficher(self):
        if self.mode != 'fichier':
            return
        self.p_spec.setXRange(0, self.sig.duree, padding=0)
        self.p_spec.setYRange(0, self.fmax, padding=0)
        self.rendre_vue()

    def rendre_vue(self):
        if self.mode != 'fichier' or self.sig is None:
            return
        (t0, t1), _ = self.p_spec.viewRange()
        i0 = max(0, int(t0 * self.fs))
        i1 = min(self.sig.n, int(math.ceil(t1 * self.fs)))
        if i1 - i0 < 2:
            return
        db, tc0, dt = self.sig.stft_zone(self.nfft, i0, i1, self.cal)
        df = self.fs / self.nfft
        self.img.setImage(db.T, autoLevels=False)
        self.img.setRect(QtCore.QRectF(tc0 - dt / 2, -df / 2, len(db) * dt, db.shape[1] * df))
        self.c_wave.setData(*self.sig.forme_onde(i0, i1))

    def maj_curseur(self):
        if self.mode != 'fichier' or self.sig is None:
            return
        t = self.curseur.value()
        trame = self.sig.trames(np.array([int(round(t * self.fs))]), self.nfft)[0]
        db = spectre_db(trame, np.hanning(self.nfft), self.cal)
        self.c_inst.setData(self.freqs, db)
        self.curseur_wave.setValue(t)
        f, lp = pic(self.freqs, db)
        etat = (" <span style='color:#c62828'>(non enregistrée)</span>"
                if self.non_sauve else "")
        self.info(f"<b>{self.nom}</b>{etat}  |  {fmt_duree(self.sig.duree)} à {self.fs} Hz  |  "
                  f"curseur t = {t:.3f} s: niveau trame {niveau_db(trame, self.cal):6.1f} dB, "
                  f"pic {f:7.1f} Hz à {lp:6.1f} dB")

    def changer_nfft(self, txt):
        self.nfft = int(txt)
        self.maj_resolution()
        if self.mode == 'fichier':
            self.rendre_vue()
            self.maj_curseur()

    def closeEvent(self, ev):
        if self.mode == 'direct':
            self.arreter()
        if not self.confirmer_perte():
            ev.ignore()
            return
        self.liberer_temp()
        ev.accept()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fichier', nargs='?', help=".wav à ouvrir au démarrage")
    ap.add_argument('--nfft', type=int, default=2048, choices=NFFT_CHOIX)
    ap.add_argument('--fmax', type=float, default=None, help="fréquence max affichée [Hz]")
    ap.add_argument('--fenetre', type=int, default=10, choices=FENETRES,
                    help="durée visible en direct [s]")
    ap.add_argument('--cal', type=float, default=0.0, help="décalage d'étalonnage [dB]")
    ap.add_argument('--min-db', type=float, default=-110)
    ap.add_argument('--max-db', type=float, default=-20)
    args = ap.parse_args()

    app = pg.mkQApp("Piccolo II")
    win = Fenetre(args)
    win.resize(1250, 900)
    win.show()
    if args.fichier:
        win.charger(args.fichier)
    sys.exit(app.exec() if hasattr(app, 'exec') else app.exec_())


if __name__ == '__main__':
    main()
