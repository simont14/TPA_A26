# Piccolo II en flux audio continu

Lit directement l'audio brut (PCM 16 bits, mono, 44 117 Hz) du sonomètre
Piccolo II par USB, sans Piccolo Utility, et affiche un spectrogramme en
temps réel sans moyenne.

Dépendances: `pip install pyusb libusb-package pyqtgraph PyQt5 numpy`

Fermer Piccolo Utility (`PIC2_*.exe`) avant de lancer un script en direct.

| Fichier | Rôle |
|---|---|
| `piccolo_usb.py` | Accès USB et flux audio. `python piccolo_usb.py 5 essai.wav` enregistre 5 s. |
| `piccolo_spectrogramme.py` | Interface: acquisition en direct (Démarrer/Arrêter) et visualisation de .wav. `--help` pour les options. |
| `sifflement_1.wav` | Enregistrement fait par l'Audio Recorder officiel, qui a servi à valider le décodage. |

```
python piccolo_spectrogramme.py
python piccolo_spectrogramme.py sifflement_1.wav       # ouvre un fichier au démarrage
python piccolo_spectrogramme.py --nfft 8192 --fmax 5000 --fenetre 60
```

- Démarrer/Arrêter: aucune limite de durée. L'audio est écrit sur le disque
  au fur et à mesure (`%TEMP%\piccolo_<date>.wav`, environ 320 Mo par heure),
  donc la mémoire ne grossit pas, et le fichier reste lisible si le
  programme plante. Le spectrogramme défile sur la « Fenêtre en direct »
  choisie (10 s à 10 min). À l'arrêt, tout l'enregistrement s'affiche et
  « Enregistrer .wav » déplace le fichier temporaire (dossier
  `acoustique/data/` par défaut). Une confirmation est demandée avant de
  perdre une acquisition non enregistrée.
- Ouvrir .wav: spectrogramme complet du fichier, lu directement sur le
  disque (un fichier de 30 min s'ouvre instantanément). Molette pour zoomer,
  glisser pour se déplacer (le calcul est refait sur la zone visible),
  « Tout afficher » pour revenir. Le curseur bleu donne le spectre de la
  trame à cet instant et la fréquence du pic.
- nfft: change la résolution (df = fs/nfft) sans recharger.
- Avec une longue fenêtre en direct ou une vue très dézoomée, une trame sur
  plusieurs est affichée (au plus 2500 ou 1500 colonnes). Zoomer pour voir
  toutes les trames. Les données enregistrées sont toujours complètes.

## Protocole (rétro-ingénierie, 5 novembre 2026)

Obtenu par capture USBPcap/tshark de `PIC2_WRC.exe` pendant un enregistrement.

- Périphérique Signal Ranger mk3 (VID `0x1612`, PID `0x0108`), pilote WinUSB,
  aucun endpoint: seulement des transferts de contrôle vendor sur EP0.
- `bRequest 0xA0`: lecture (`0xC0`) / écriture (`0x40`) de la mémoire du
  STM32L476, adresse = `wIndex << 16 | wValue`, 4096 octets max par transfert.
- `bRequest 0xB4`: appel de la fonction du firmware à cette adresse.
- Table de symboles du firmware à `0x08008000` (723 noms + adresses), lue au
  démarrage pour trouver les variables par leur nom.
- Démarrage/arrêt: appel `USB_Init_RecordTask`, attente `USBAppEnd == 1`,
  puis `Gbl_RecIn = 1` au démarrage.
- Boucle: `FiFoSize` (échantillons disponibles) ≥ 4096 → lire 8192 octets dans
  `DAQBuffer` (anneau de 32 ko) → appeler `UpdateOutIndex`. `MissData` compte
  les débordements.

Validation: le flux USB capturé est identique, échantillon par échantillon,
à `sifflement_1.wav`. Les jonctions entre blocs ne montrent aucune
discontinuité.

Les niveaux sont en dBFS. `--cal` ajoute un décalage en dB, à ajuster
contre la lecture LZeq de l'appareil (dépend de la plage Low/High).
