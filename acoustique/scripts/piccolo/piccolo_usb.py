"""Accès direct au sonomètre Piccolo II (Soft dB) par USB, sans Piccolo Utility.

Protocole (obtenu par capture USBPcap du logiciel « Audio Recorder », PIC2_WRC.exe):

- Le Piccolo II s'annonce comme une carte « Signal_Ranger_mk3 » (VID 0x1612,
  PID 0x0108) avec le pilote WinUSB. L'interface n'a aucun endpoint: tout passe
  par des transferts de contrôle vendor sur l'endpoint 0.
- bRequest 0xA0: lecture (bmRequestType 0xC0) ou écriture (0x40) de la mémoire
  du STM32L476. Adresse 32 bits = (wIndex << 16) | wValue. Au plus 4096 octets
  par transfert.
- bRequest 0xB4 (0x40, wLength 0): appel de la fonction du firmware située à
  l'adresse donnée.
- Le firmware contient une table de symboles à 0x08008000 (nom de version,
  nombre de symboles, noms terminés par 0 et alignés sur 4 octets, puis une
  table d'adresses 32 bits). Les adresses des variables sont lues là.

Séquence de l'enregistrement audio:
  1. appel USB_Init_RecordTask, attendre USBAppEnd == 1, écrire Gbl_RecIn = 1
  2. boucle: lire FiFoSize (nombre d'échantillons disponibles). Dès qu'il y en
     a 4096, lire 8192 octets dans DAQBuffer (anneau de 32 ko) puis appeler
     UpdateOutIndex pour libérer le bloc. MissData compte les débordements.
  3. arrêt: appel USB_Init_RecordTask, attendre USBAppEnd == 1.

Les échantillons sont du PCM 16 bits signé little-endian, mono, à 44 117 Hz.
Validation: les données lues par USB sont identiques, échantillon par
échantillon, au fichier .wav sauvé par PIC2_WRC (sifflement_1.wav).

Utilisation en ligne de commande (enregistre un .wav de test):
    python piccolo_usb.py 5 test.wav
"""

import struct
import sys
import time
import wave

import numpy as np
import usb.core
import usb.util
import libusb_package

VID, PID = 0x1612, 0x0108
REQ_MEM = 0xA0
REQ_CALL = 0xB4
MAX_CTRL = 4096
SYMTAB_ADDR = 0x08008000

FS = 44117                 # fréquence d'échantillonnage réelle [Hz]
CHUNK_SAMPLES = 4096       # taille du bloc libéré par UpdateOutIndex
CHUNK_BYTES = 2 * CHUNK_SAMPLES
RING_BYTES = 32768         # taille de DAQBuffer


def piccolo_present():
    """Vrai si Windows voit un Piccolo II branché."""
    backend = libusb_package.get_libusb1_backend()
    return usb.core.find(idVendor=VID, idProduct=PID, backend=backend) is not None


class PiccoloII:
    """Connexion bas niveau: lecture/écriture mémoire, appels, symboles."""

    def __init__(self, timeout_ms=1000):
        backend = libusb_package.get_libusb1_backend()
        self.dev = usb.core.find(idVendor=VID, idProduct=PID, backend=backend)
        if self.dev is None:
            raise RuntimeError("Piccolo II introuvable (branché? Piccolo Utility fermé?)")
        self.timeout = timeout_ms
        self.version, self.sym = self._read_symbols()

    # -- accès mémoire -------------------------------------------------------
    def read(self, addr, n):
        out = bytearray()
        while n > 0:
            k = min(n, MAX_CTRL)
            out += self.dev.ctrl_transfer(0xC0, REQ_MEM, addr & 0xFFFF, addr >> 16, k,
                                          timeout=self.timeout)
            addr += k
            n -= k
        return bytes(out)

    def write(self, addr, data):
        self.dev.ctrl_transfer(0x40, REQ_MEM, addr & 0xFFFF, addr >> 16, data,
                               timeout=self.timeout)

    def call(self, addr):
        self.dev.ctrl_transfer(0x40, REQ_CALL, addr & 0xFFFF, addr >> 16, None,
                               timeout=self.timeout)

    def read_u8(self, addr):
        return self.read(addr, 1)[0]

    def read_u32(self, addr):
        return struct.unpack('<I', self.read(addr, 4))[0]

    # -- table de symboles ---------------------------------------------------
    def _read_symbols(self):
        head = self.read(SYMTAB_ADDR, 0x24)
        version = head[4:head.index(b'\0', 4)].decode()
        n = struct.unpack_from('<I', head, 0x20)[0]
        names, addr, buf = [], SYMTAB_ADDR + 0x24, b''
        while len(names) < n:
            buf += self.read(addr + len(buf), 512)
            off = 0
            names = []
            while len(names) < n:
                end = buf.find(b'\0', off)
                if end < 0:
                    break
                names.append(buf[off:end].decode('latin-1'))
                off = (end + 1 + 3) & ~3
        table = self.read(addr + off, 4 * n)
        addrs = struct.unpack('<%dI' % n, table)
        return version, dict(zip(names, addrs))

    def __getitem__(self, name):
        return self.sym[name]

    def close(self):
        usb.util.dispose_resources(self.dev)


class PiccoloStream:
    """Flux audio continu. Utiliser avec `with` pour toujours arrêter le flux."""

    def __init__(self, pic: PiccoloII):
        self.p = pic
        self.out_off = 0
        self.missed = 0
        self.overflows = 0     # nombre de débordements depuis start()
        self.running = False

    def _toggle_record_task(self):
        p = self.p
        p.call(p['USB_Init_RecordTask'])
        t0 = time.time()
        while p.read_u8(p['USBAppEnd']) != 1:
            if time.time() - t0 > 2:
                raise TimeoutError("USBAppEnd n'est jamais passé à 1")
            time.sleep(0.005)

    def start(self):
        self._toggle_record_task()
        self.p.write(self.p['Gbl_RecIn'], struct.pack('<H', 1))
        self.out_off = 0
        self.missed = self.p.read_u32(self.p['MissData'])
        self.running = True

    def stop(self):
        if self.running:
            self.running = False
            self._toggle_record_task()

    def read_chunk(self, poll_s=0.005, timeout_s=1.0):
        """Bloque jusqu'à avoir un bloc de 4096 échantillons (int16)."""
        p = self.p
        t0 = time.time()
        while p.read_u32(p['FiFoSize']) < CHUNK_SAMPLES:
            if time.time() - t0 > timeout_s:
                raise TimeoutError("aucune donnée audio reçue du Piccolo")
            time.sleep(poll_s)
        raw = p.read(p['DAQBuffer'] + self.out_off, CHUNK_BYTES)
        p.call(p['UpdateOutIndex'])
        self.out_off = (self.out_off + CHUNK_BYTES) % RING_BYTES
        miss = p.read_u32(p['MissData'])
        if miss != self.missed:
            print(f"[piccolo] débordement de l'anneau: MissData {self.missed} -> {miss}",
                  file=sys.stderr)
            self.missed = miss
            self.overflows += 1
        return np.frombuffer(raw, '<i2')

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *exc):
        self.stop()


def record_wav(seconds, path):
    pic = PiccoloII()
    print(f"Firmware: {pic.version} ({len(pic.sym)} symboles)")
    blocks = []
    with PiccoloStream(pic) as s:
        t0 = time.time()
        while len(blocks) * CHUNK_SAMPLES < seconds * FS:
            blocks.append(s.read_chunk())
        dt = time.time() - t0
    pic.close()
    x = np.concatenate(blocks)
    with wave.open(path, 'wb') as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(FS)
        w.writeframes(x.tobytes())
    print(f"{len(x)} échantillons en {dt:.2f} s (débit {len(x) / dt:.0f} éch/s), "
          f"crête {np.abs(x).max()}, MissData {s.missed}")


if __name__ == '__main__':
    record_wav(float(sys.argv[1]) if len(sys.argv) > 1 else 3,
               sys.argv[2] if len(sys.argv) > 2 else 'piccolo_test.wav')
