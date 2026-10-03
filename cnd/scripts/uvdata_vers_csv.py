#!/usr/bin/env python3
"""
uvdata_vers_csv.py : convertit un fichier .UVData (TOPAZ / UltraVision) en CSV, sans UltraVision et sans aucun
traitement des donnees. Les valeurs sont celles enregistrees dans le fichier.

Une ligne par A-scan :
    position, faisceau, e0, e1, e2, ...

    position   numero de la position, tel qu'enregistre par le Topaz (avec un encodeur : la case de position de
               l'encodeur; sans encodeur : le coup d'horloge)
    faisceau   numero du faisceau
    e0, e1...  valeur brute de chaque echantillon de l'A-scan, dans l'ordre du temps
               (cellules vides si le Topaz n'a enregistre aucun A-scan a cette position)

Avec --format long, une ligne par echantillon :
    position, faisceau, echantillon, valeur

Utilisation
    python uvdata_vers_csv.py                                   (fenetres pour choisir les fichiers)
    python uvdata_vers_csv.py mesure.UVData                     (ecrit mesure.csv a cote)
    python uvdata_vers_csv.py mesure.UVData -o sortie.csv --format long --sep ";"

Seul NumPy est necessaire (pip install numpy). Le fichier .UVData est lu en lecture seule. Son format a ete
deduit de fichiers UltraVision Touch 3.8R11 (TOPAZ16, balayage lineaire).
"""
import argparse
import os
import re
import struct
import sys

import numpy as np

BS = 0x8000      # taille d'un bloc du conteneur UltraVision
HDR = 0x20       # en-tete d'un bloc


# =====================================================================  lecture du conteneur
class UVContainer:
    """Mini systeme de fichiers d'UltraVision (blocs de 0x8000 octets)."""

    def __init__(self, path):
        with open(path, 'rb') as f:
            self.d = f.read()
        self.nb = len(self.d) // BS

    def _hdr(self, b):
        return struct.unpack('<8I', self.d[b * BS:b * BS + 32])

    def _type(self, b):
        return struct.unpack('<I', self.d[b * BS:b * BS + 4])[0]

    def _read(self, ino):
        rec = struct.unpack('<5I3dI', self.d[ino * BS + HDR:ino * BS + HDR + 48])
        recsize, size = rec[0], rec[8]
        tags = struct.unpack('<20I', self.d[ino * BS + HDR:ino * BS + HDR + 80])
        if tags[13] == 0x000C0002 and recsize in (96, 108):          # grand fichier : chaine de blocs
            b, parts, got = tags[15], [], 0
            while b and got < size:
                if self._type(b) != 0x22:
                    raise ValueError('chaine de blocs interrompue au bloc %d' % b)
                parts.append(self.d[b * BS + HDR:(b + 1) * BS])
                got += BS - HDR
                b = self._hdr(b)[4]
            return b''.join(parts)[:size]
        base = ino * BS + HDR                                          # petit fichier : dans l'inode
        return self.d[base + 96:base + 96 + size]

    def files(self):
        conts = {}
        for b in range(self.nb):
            if self._type(b) == 0x1f:
                s = self.d[b * BS + HDR:(b + 1) * BS].decode('utf-16le', 'ignore').split('\x00')[0]
                conts.setdefault(self._hdr(b)[5], []).append(s)
        names, parent, kinds = {}, {}, {}
        for db in range(self.nb):
            if self._type(db) != 0x21:
                continue
            owner = self._hdr(db)[3]
            off = db * BS + HDR
            while off + 64 <= (db + 1) * BS:
                nid, blk, kind, nlen, _ = struct.unpack('<5I', self.d[off:off + 20])
                if nid == 0 and blk == 0:
                    break
                nm = self.d[off + 20:off + 64].decode('utf-16le', 'ignore').split('\x00')[0]
                if nlen > len(nm):
                    for s in conts.get(nid, []):
                        nm += s
                names[blk], parent[blk], kinds[blk] = nm[:nlen], owner, kind
                off += 64

        def path(b):
            p = []
            while b in names:
                p.append(names[b])
                b = parent.get(b)
            return '/'.join(reversed(p))

        return {path(b): self._read(b) for b, k in kinds.items() if k == 1}


def _tag(xml, name, cast=str, default=None):
    m = re.search(r'<%s[^>]*>([^<]*)</%s>' % (name, name), xml)
    if not m:
        return default
    try:
        return cast(m.group(1))
    except ValueError:
        return default

def lire_uvdata(path):
    """Retourne (brut, present, horloge) : brut[position, faisceau, echantillon] = valeurs enregistrees,
    present[position, faisceau] vrai si un A-scan a ete enregistre, horloge vrai si la position vient de l'horloge."""
    files = UVContainer(path).files()
    txt = {k: v.decode('utf-8', 'ignore') for k, v in files.items() if k.endswith('.xml')}
    hw = txt.get('Setups/HardwareSetup.xml', '')
    blob_dir = next((k.rsplit('/', 1)[0] for k in files if k.endswith('/Blob.bin')), None)
    if blob_dir is None:
        raise ValueError("Ce fichier ne contient pas de données ultrasons (pas de Blob.bin).")
    blob = files[blob_dir + '/Blob.bin']
    table = files[blob_dir + '/Table.bin']
    lookup = re.findall(r'<Item>(.*?)</Item>', txt[blob_dir + '/Data Descriptor.xml'], re.S)

    for d in sorted({k.rsplit('/', 1)[0] for k in txt if k.startswith('UT Data/')}):
        desc = txt.get(d + '/Data Descriptor.xml', '')
        defi = txt.get(d + '/Data Definition.xml', '')
        if 'Zetec.Ultrasound.DAL.AscanData,' not in desc:
            continue
        short = _tag(defi, 'DataPath').split('/', 1)[1]
        item = next(i for i in lookup if '<Path>%s</Path>' % short in i)
        t_off, ncell, csize = (_tag(item, k, int) for k in ('Offset', 'CellQuantity', 'CellSize'))
        nx, ny = _tag(desc, 'LimitX', int), _tag(desc, 'LimitY', int)
        ns, ssize = _tag(desc, 'SampleQuantity', int), _tag(desc, 'SampleSize', int)
        signed = _tag(desc, 'IsSigned', str, 'false') == 'true'
        if ncell != nx * ny or csize != ns * ssize:
            raise ValueError('Format A-scan inattendu dans ' + d)
        offs = np.frombuffer(table[t_off:t_off + 8 * ncell], dtype='<u8').reshape(ny, nx)
        dtype = {1: 'i1' if signed else 'u1', 2: '<i2' if signed else '<u2'}[ssize]
        brut = np.zeros((nx, ny, ns), np.int64)
        present = np.zeros((nx, ny), bool)
        for y in range(ny):
            for x in range(nx):
                o = int(offs[y, x])
                if o > 0 and o + csize <= len(blob):                  # offset 0 : aucun A-scan a cette position
                    brut[x, y] = np.frombuffer(blob, dtype=dtype, count=ns, offset=o)
                    present[x, y] = True
        return brut, present, _tag(hw, 'UseInternalEncoder', str, '') == 'true'
    raise ValueError("Aucun A-scan enregistré : Record Condition était peut-être en C-Scan seulement.")


# =====================================================================  ecriture du CSV
def ecrire_csv(brut, present, sortie, fmt='large', sep=','):
    """Ecrit toutes les positions et tous les faisceaux du fichier, sans aucun traitement."""
    nx, ny, ns = brut.shape
    with open(sortie, 'w', encoding='utf-8', newline='') as f:
        if fmt == 'long':
            f.write(sep.join(['position', 'faisceau', 'echantillon', 'valeur']) + '\n')
            for p in range(nx):
                for b in range(ny):
                    if present[p, b]:
                        f.write(''.join('%d%s%d%s%d%s%d\n' % (p, sep, b, sep, k, sep, v)
                                        for k, v in enumerate(brut[p, b])))
        else:
            f.write(sep.join(['position', 'faisceau'] + ['e%d' % k for k in range(ns)]) + '\n')
            vide = sep * ns
            for p in range(nx):
                for b in range(ny):
                    corps = sep + sep.join(map(str, brut[p, b])) if present[p, b] else vide
                    f.write('%d%s%d%s\n' % (p, sep, b, corps))


def resume(source, sortie, brut, present, horloge):
    nx, ny, ns = brut.shape
    vides = int((~present.any(axis=1)).sum())
    return ('%s\n-> %s\n\n%d positions x %d faisceaux x %d échantillons (%d positions sans A-scan)\n'
            'Position : %s' % (os.path.basename(source), sortie, nx, ny, ns, vides,
                               "coup d'horloge (pas d'encodeur)" if horloge else 'encodeur'))


# =====================================================================  lancement
def fenetres():
    """Sans argument : fenetres pour choisir le .UVData et l'emplacement du CSV."""
    import tkinter as tk
    from tkinter import filedialog, messagebox
    root = tk.Tk()
    root.withdraw()
    source = filedialog.askopenfilename(title='Choisir le fichier de mesure',
                                        filetypes=[('Données UltraVision', '*.UVData'), ('Tous les fichiers', '*.*')])
    if not source:
        return
    sortie = filedialog.asksaveasfilename(title='Enregistrer le CSV', initialdir=os.path.dirname(source),
                                          initialfile=os.path.splitext(os.path.basename(source))[0] + '.csv',
                                          defaultextension='.csv', filetypes=[('CSV', '*.csv')])
    if not sortie:
        return
    try:
        root.config(cursor='watch')
        root.update()
        brut, present, horloge = lire_uvdata(source)
        ecrire_csv(brut, present, sortie)
    except Exception as e:
        messagebox.showerror('Conversion impossible', str(e))
        return
    messagebox.showinfo('Conversion terminée', resume(source, sortie, brut, present, horloge))


def main():
    if len(sys.argv) == 1:
        fenetres()
        return
    ap = argparse.ArgumentParser(description="Convertit un fichier .UVData (TOPAZ / UltraVision) en CSV, "
                                             "sans aucun traitement des données.")
    ap.add_argument('fichier', help='fichier .UVData à convertir')
    ap.add_argument('-o', '--sortie', help='fichier CSV à écrire (par défaut : même nom, extension .csv)')
    ap.add_argument('--format', choices=('large', 'long'), default='large',
                    help='large : une ligne par A-scan (défaut). long : une ligne par échantillon')
    ap.add_argument('--sep', default=',', help='séparateur de colonnes (défaut : ,)')
    a = ap.parse_args()
    sortie = a.sortie or os.path.splitext(a.fichier)[0] + '.csv'
    brut, present, horloge = lire_uvdata(a.fichier)
    ecrire_csv(brut, present, sortie, a.format, a.sep)
    print(resume(a.fichier, sortie, brut, present, horloge))


if __name__ == '__main__':
    main()
