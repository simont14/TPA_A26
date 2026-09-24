#!/usr/bin/env python3
"""
uvdata_3d_gui.py : reconstruction 3D d'une acquisition TOPAZ / UltraVision (.UVData), sans UltraVision.

Deroulement :
  1. une fenetre te demande le fichier de mesure (.UVData)
  2. une fenetre te demande les parametres du balayage manuel
  3. une fenetre te demande ou enregistrer le fichier 3D :
       .html  vue 3D interactive, s'ouvre dans un navigateur   (demande scipy et scikit-image)
       .vtk   volume pour ParaView (logiciel gratuit)
       .npz   volume NumPy pour ton propre traitement

Installation :
    pip install numpy scipy scikit-image
Lancement :
    python uvdata_3d_gui.py

Le fichier .UVData est lu en lecture seule. Son format a ete deduit d'un fichier
UltraVision Touch 3.8R11 (TOPAZ16, balayage lineaire) ; aucune protection de licence n'est touchee.
"""
import base64
import html
import json
import math
import os
import re
import struct
import sys
import webbrowser

import numpy as np

BS = 0x8000
HDR = 0x20


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


def read_calibration(hw, comp):
    """Etat des calibrations enregistre dans le setup : etapes terminees, cible d'epaisseur,
    et corrections reellement appliquees aux faisceaux."""
    out = dict(cal_done={}, cal_thickness_mm=None, cal_applied=False)
    sec = re.search(r'<Calibration xsi:type="ultrasound:Calibration".*?</Calibration>', hw, re.S)
    if sec:
        sec = sec.group(0)
        for name, label in (('Velocity', 'vitesse'), ('WedgeDelay', 'délai du sabot'),
                            ('Sensitivity', 'sensibilité'), ('TGC', 'TCG'), ('ElementCheck', 'éléments')):
            m = re.search(r'<%s[\s>].*?</%s>' % (name, name), sec, re.S)
            if m:
                out['cal_done'][label] = _tag(m.group(0), 'IsCompleted', str, '') == 'true'
                if name in ('Velocity', 'WedgeDelay') and out['cal_done'][label] and not out['cal_thickness_mm']:
                    loc = _tag(m.group(0), 'Location', float)
                    if loc:
                        out['cal_thickness_mm'] = round(loc * 1000, 3)
    for t in ('BeamLongitudinalVelocityCorrection', 'WedgeDelayCorrection', 'CalibrationGain'):
        for v in re.findall(r'<%s>([^<]*)</%s>' % (t, t), hw + comp):
            try:
                if abs(float(v)) > 0:
                    out['cal_applied'] = True
            except ValueError:
                pass
    return out


def read_acquisition(path):
    """Retourne (info generale, [(volume % ecran, meta), ...])."""
    files = UVContainer(path).files()
    txt = {k: v.decode('utf-8', 'ignore') for k, v in files.items() if k.endswith('.xml')}
    hw = txt.get('Setups/HardwareSetup.xml', '')
    acq = txt.get('Setups/AcquisitionInformation.xml', '')
    comp = txt.get('Setups/ComponentSetup.xml', '')
    spec = txt.get('Setups/SpecimenSetup.xml', '')
    stop = _tag(hw, 'StopAtEndOfScanMode', str, None)
    info = dict(
        file=os.path.basename(path),
        setup=_tag(acq, 'SetupFileName', str, ''),
        date=(_tag(acq, 'AcquisitionDate', str, '') or '')[:16].replace('T', ' '),
        device=_tag(acq, 'FriendlyName', str, ''),
        rate_hz=_tag(hw, 'AcquisitionRate', float),
        rate_mode=_tag(hw, 'AcquisitionRate_AutoBehavior', str, ''),   # UserValue = valeur fixee par toi
        internal_clock=_tag(hw, 'UseInternalEncoder', str, '') == 'true',
        stop_at_end=None if stop is None else stop == 'true',      # false : le TOPAZ reecrit le debut (tampon)
        specimen_thickness_mm=(_tag(spec, 'Thickness', float) or 0) * 1000 or None,
    )
    info.update(read_calibration(hw, comp))

    blob_dir = next((k.rsplit('/', 1)[0] for k in files if k.endswith('/Blob.bin')), None)
    if blob_dir is None:
        raise ValueError("Ce fichier ne contient pas de données ultrasons (pas de Blob.bin).")
    blob = files[blob_dir + '/Blob.bin']
    table = files[blob_dir + '/Table.bin']
    lookup = re.findall(r'<Item>(.*?)</Item>', txt[blob_dir + '/Data Descriptor.xml'], re.S)

    vols = []
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
        vmax = float(np.iinfo(np.dtype(dtype)).max)
        vol = np.zeros((nx, ny, ns), np.float32)
        for y in range(ny):
            for x in range(nx):
                o = int(offs[y, x])
                if 0 < o + csize <= len(blob):
                    vol[x, y] = np.frombuffer(blob, dtype=dtype, count=ns, offset=o)
        vol *= 100.0 / vmax
        wave = _tag(defi, 'CurrentWaveType', str, 'Longitudinal')
        meta = dict(
            beam=_tag(defi, 'BeamName', str, ''),
            scan_res_mm=_tag(defi, 'ScanSamplingResolution', float) * 1000,
            index_res_mm=_tag(defi, 'IndexSamplingResolution', float) * 1000,
            sample_period_s=_tag(defi, 'UsoundSamplingResolution', float),
            velocity_m_s=_tag(defi, 'SpecimenLongitudinalSoundSpeed' if wave == 'Longitudinal'
                              else 'SpecimenTransversalSoundSpeed', float),
            wave_type=wave,
        )
        vols.append((vol, meta))
    if not vols:
        raise ValueError("Aucun A-scan enregistré : Record Condition était peut-être en C-Scan seulement.")
    return info, vols


# =====================================================================  traitement
def detect_wrap(vol):
    """Tampon circulaire : en mode horloge, si l'acquisition dure plus longtemps que Scan Start-Stop,
    le TOPAZ reecrit le debut. On cherche une coupure franche alors que la derniere et la premiere
    position se raccordent normalement."""
    nx = vol.shape[0]
    if nx < 10:
        return None
    f = vol.reshape(nx, -1)
    d = np.linalg.norm(np.diff(f, axis=0), axis=1)
    med = float(np.median(d)) or 1e-9
    k = int(np.argmax(d))
    jump, wrapd = float(d[k]), float(np.linalg.norm(f[0] - f[-1]))
    if jump > 8 * med and wrapd < 4 * med and jump > 5 * wrapd:
        return dict(index=k + 1, ratio=jump / med)
    return None


def _centroid_region(p, k):
    """Centroide de la zone contigue autour de k au-dessus de 50 % du pic (gere la saturation)."""
    lim = 0.5 * p[k]
    i0 = k
    while i0 > 0 and p[i0 - 1] >= lim:
        i0 -= 1
    i1 = k
    while i1 + 1 < len(p) and p[i1 + 1] >= lim:
        i1 += 1
    idx = np.arange(i0, i1 + 1)
    return float((idx * p[i0:i1 + 1]).sum() / p[i0:i1 + 1].sum())


def find_surface(vol):
    p = vol[:, :, : vol.shape[2] // 4].mean(axis=(0, 1))
    i0 = int(np.argmax(p >= 0.5 * p.max()))              # debut du premier echo fort
    i1 = i0
    while i1 + 1 < len(p) and p[i1 + 1] >= 0.5 * p.max():
        i1 += 1
    return _centroid_region(p, i0 + int(np.argmax(p[i0:i1 + 1])))


def find_backwall(vol, z0, dz, thickness):
    ns = vol.shape[2]
    e = z0 + thickness / dz
    lo, hi = int(z0 + 0.8 * (e - z0)), int(min(ns - 1, z0 + 1.2 * (e - z0)))
    if hi - lo < 3 or lo >= ns:
        return None
    p = np.median(vol.reshape(-1, ns), axis=0)                          # robuste aux zones masquees
    k = lo + int(np.argmax(p[lo:hi + 1]))
    return _centroid_region(p, k) if p[k] > 5 else None


def process(vol, meta, params, wrap):
    v = vol
    notes = []
    if params['unwrap'] and wrap:
        v = np.roll(v, -wrap['index'], axis=0)
        notes.append('Tampon circulaire remis dans l\'ordre : la coupure était à la position %d.' % wrap['index'])
    if params['reverse']:
        v = v[::-1]
        notes.append('Axe du scan retourné (balayage fait en sens inverse).')
    v = np.ascontiguousarray(v)
    dx = params['speed'] / params['rate']
    dy = meta['index_res_mm']
    vel = params['velocity']
    dz = meta['sample_period_s'] * vel / 2 * 1000
    z0 = find_surface(v)
    calib = None
    if params['thickness']:
        kb = find_backwall(v, z0, dz, params['thickness'])
        if kb is None:
            notes.append("Écho de fond introuvable près de l'épaisseur donnée : vitesse non recalée.")
        else:
            measured = (kb - z0) * dz
            vel = vel * params['thickness'] / measured
            dz = meta['sample_period_s'] * vel / 2 * 1000
            calib = dict(measured_mm=measured, velocity=vel)
    nx, ny, ns = v.shape
    return dict(v=v, dx=dx, dy=dy, dz=dz, z0=z0, vel=vel, calib=calib, notes=notes,
                speed=params['speed'], rate=params['rate'],
                X=(nx - 1) * dx, Y=(ny - 1) * dy, zmin=-z0 * dz, zmax=(ns - 1 - z0) * dz,
                thickness=params['thickness'])


# =====================================================================  exports
def write_vtk(path, R):
    v = np.clip(R['v'], 0, 100).astype('>f4')
    nx, ny, nz = v.shape
    with open(path, 'wb') as f:
        f.write(b'# vtk DataFile Version 3.0\nUltraVision A-scan volume (amplitude % ecran)\nBINARY\n')
        f.write(b'DATASET STRUCTURED_POINTS\n')
        f.write(('DIMENSIONS %d %d %d\n' % (nx, ny, nz)).encode())
        f.write(('SPACING %g %g %g\n' % (R['dx'], R['dy'], R['dz'])).encode())
        f.write(('ORIGIN 0 0 %g\n' % R['zmin']).encode())
        f.write(('POINT_DATA %d\nSCALARS amplitude_pct float 1\nLOOKUP_TABLE default\n' % v.size).encode())
        f.write(v.transpose(2, 1, 0).tobytes())


def write_npz(path, R, info):
    nx, ny, ns = R['v'].shape
    np.savez_compressed(
        path, volume=R['v'].astype(np.float32),
        scan_mm=np.arange(nx) * R['dx'], index_mm=np.arange(ny) * R['dy'],
        depth_mm=(np.arange(ns) - R['z0']) * R['dz'],
        info=json.dumps(dict(info, velocity_m_s=R['vel'], notes=R['notes'])))


def _b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode('ascii')


def write_html(path, R, info, texts=None):
    try:
        from scipy import ndimage
        from skimage import measure
    except ImportError:
        raise RuntimeError("La vue .html demande scipy et scikit-image :\n    pip install scipy scikit-image\n"
                           "Ou choisis plutôt .vtk (ParaView) ou .npz.")
    v = R['v']
    nx, ny, ns = v.shape
    sm = ndimage.gaussian_filter(v, sigma=(1.5, 0.7, 1.0))
    sx, sz = max(1, math.ceil(nx / 140)), max(1, math.ceil(ns / 130))
    g = sm[::sx, :, ::sz]
    q = lambda a, lo, hi: np.round((a - lo) / ((hi - lo) or 1) * 65535).clip(0, 65535).astype('<u2')
    meshes = {}
    for lvl in (15, 30, 50):
        if not (g.min() < lvl < g.max()):
            continue
        verts, faces, _, _ = measure.marching_cubes(g, level=lvl, spacing=(R['dx'] * sx, R['dy'], R['dz'] * sz))
        dep = verts[:, 2] - R['z0'] * R['dz']
        vb = np.stack([q(verts[:, 0], 0, R['X']), q(verts[:, 1], 0, R['Y']), q(dep, R['zmin'], R['zmax'])], 1)
        fw = 2 if len(verts) < 65536 else 4
        fb = faces.astype('<u2' if fw == 2 else '<u4')
        meshes[str(lvl)] = dict(v=_b64(vb), f=_b64(fb), nv=len(verts), nf=len(faces), fw=fw)
    if not meshes:
        raise RuntimeError('Aucune surface à 15, 30 ou 50 % : le signal est peut-être vide.')

    ss = max(1, math.ceil(nx / 250))
    s8 = np.round(v[::ss].clip(0, 100) / 100 * 255).astype(np.uint8)
    thick = R['thickness']
    M = dict(nx=s8.shape[0], ny=ny, ns=ns, dx=R['dx'] * ss, dy=R['dy'], dz=R['dz'], z0=R['z0'],
             X=R['X'], Y=R['Y'], zmin=R['zmin'], zmax=R['zmax'], cmax=thick or R['zmax'],
             zlo=3.0 if (thick or R['zmax']) > 8 else 0.0, zhi=R['zmax'],
             cx=round(R['X'] / 3, 1), cy=round(R['Y'] / 4, 1))
    data = json.dumps(dict(meta=M, meshes=meshes, slices=_b64(s8)))

    fr = lambda x: ('%.1f' % x).replace('.', ',')
    t = texts or {}
    title = t.get('title', 'Reconstruction 3D : ' + info['file'])
    fileline = t.get('fileline', ', '.join(x for x in (info['file'], info['device'], info['setup'], info['date']) if x))
    h1 = t.get('h1', 'Reconstruction 3D de ' + os.path.splitext(info['file'])[0])
    lede = t.get('lede', '%d positions × %d faisceaux × %d échantillons, lus directement dans le fichier. '
                         'Fais tourner la vue avec la souris ou le doigt, et resserre la fenêtre de profondeur '
                         'pour isoler une couche.' % (nx, ny, ns))
    notes = t.get('notes')
    if notes is None:
        notes = [('0 mm', "Surface : écho d'interface sabot/pièce.", False)]
        if thick:
            notes.append((fr(thick) + ' mm', 'Fond de la pièce (épaisseur donnée).', False))
        notes.append(('2×, 3×', "Un écho à 2 ou 3 fois la profondeur d'un autre est souvent un écho multiple, "
                                "pas un vrai défaut.", False))
    cave = t.get('caveats')
    if cave is None:
        cave = ['Pas du scan : %s mm par position (%s mm/s à %s Hz), longueur %s mm. Sans encodeur, '
                'les longueurs le long du scan dépendent de la régularité du déplacement.'
                % (('%.3f' % R['dx']).replace('.', ','), ('%g' % R['speed']).replace('.', ','),
                   ('%g' % R['rate']).replace('.', ','), fr(R['X']))]
        if R['calib']:
            cave.append("Profondeurs recalées sur l'épaisseur de %s mm : le fond sortait à %s mm, d'où une vitesse "
                        "d'environ %d m/s." % (fr(thick), fr(R['calib']['measured_mm']), round(R['calib']['velocity'])))
        else:
            cave.append('Vitesse du son utilisée : %d m/s. Le zéro est au centre de l\'écho de surface.' % round(R['vel']))
        cave += R['notes']
    dl = ''.join('<dt>%s</dt><dd%s>%s</dd>' % (html.escape(a), ' class="key"' if k else '', html.escape(b))
                 for a, b, k in notes)
    ps = ''.join('<p>%s</p>' % html.escape(p) for p in cave)
    page = (HTML_TEMPLATE.replace('__TITLE__', html.escape(title)).replace('__FILELINE__', html.escape(fileline))
            .replace('__H1__', html.escape(h1)).replace('__LEDE__', html.escape(lede))
            .replace('__NOTES__', dl).replace('__CAVEATS__', ps).replace('__DATA__', data))
    with open(path, 'w', encoding='utf-8') as f:
        f.write(page)


def export(path, R, info, texts=None):
    ext = os.path.splitext(path)[1].lower()
    if ext == '.vtk':
        write_vtk(path, R)
    elif ext == '.npz':
        write_npz(path, R, info)
    else:
        write_html(path, R, info, texts)


# =====================================================================  interface Tk
def _num(s):
    s = (s or '').strip().replace(',', '.')
    return float(s) if s else None


def _fr(x, fmt='%g'):
    return (fmt % x).replace('.', ',')


def file_defaults(info, vol, meta):
    """Tout ce qu'on peut tirer du fichier, sans rien demander a l'utilisateur."""
    wrap = None if info['stop_at_end'] else detect_wrap(vol)
    vv = np.roll(vol, -wrap['index'], axis=0) if wrap else vol
    thick, thick_src = None, None
    if not info['cal_applied']:
        if info['cal_thickness_mm']:
            thick, thick_src = info['cal_thickness_mm'], 'cible de ta calibration de vitesse'
        elif info['specimen_thickness_mm']:
            dz_f = meta['sample_period_s'] * meta['velocity_m_s'] / 2 * 1000
            z0_f = find_surface(vv)
            kb = find_backwall(vv, z0_f, dz_f, info['specimen_thickness_mm'])
            if kb and abs((kb - z0_f) * dz_f / info['specimen_thickness_mm'] - 1) < 0.05:
                thick, thick_src = info['specimen_thickness_mm'], 'épaisseur du spécimen du setup'
    rate_ok = bool(info['rate_hz']) and info['rate_mode'] == 'UserValue'
    return dict(wrap=wrap, thickness=thick, thickness_src=thick_src, rate_ok=rate_ok)


def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    try:                                                   # texte net sur les ecrans haute resolution
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass

    root = tk.Tk()
    root.title('Reconstruction 3D UltraVision')
    root.withdraw()

    src = filedialog.askopenfilename(title='Choisis le fichier de mesure',
                                     filetypes=[('Données UltraVision', '*.UVData'), ('Tous les fichiers', '*.*')])
    if not src:
        return
    try:
        info, vols = read_acquisition(src)
        vol, meta = vols[0]
        D = file_defaults(info, vol, meta)
        # apercu avec les valeurs du fichier : vitesse recalee, profondeur couverte
        rate_f = info['rate_hz'] or 20.0
        R0 = process(vol, meta, dict(speed=meta['scan_res_mm'] * rate_f, rate=rate_f,
                                     velocity=meta['velocity_m_s'], thickness=D['thickness'],
                                     unwrap=bool(D['wrap']), reverse=False), D['wrap'])
    except Exception as e:
        messagebox.showerror('Lecture impossible', str(e))
        return
    nx, ny, ns = vol.shape

    root.deiconify()
    root.resizable(False, False)
    frm = ttk.Frame(root, padding=16)
    frm.grid(sticky='nsew')
    ttk.Label(frm, text=info['file'], font=('TkDefaultFont', 12, 'bold')).grid(row=0, column=0, sticky='w')

    # ------------------------------------------------ 1. lu dans le fichier (lecture seule)
    box = ttk.LabelFrame(frm, text=' Lu dans le fichier ', padding=(12, 8))
    box.grid(row=1, column=0, sticky='ew', pady=(10, 0))
    rows = []

    def section(title):
        rows.append((title, None))

    def item(label, value):
        rows.append((label, value))

    section('Acquisition')
    item('Appareil, date', ', '.join(x for x in (info['device'], info['date']) if x) or '?')
    if info['setup']:
        item('Setup', info['setup'])
    item('Volume', '%d positions × %d faisceaux × %d échantillons' % (nx, ny, ns))
    item('Faisceaux', '%s, onde %s' % (meta['beam'] or '?', {'Longitudinal': 'longitudinale',
         'Transversal': 'transversale'}.get(meta['wave_type'], meta['wave_type'].lower())))
    item('Index', '%s mm entre faisceaux, %s mm de couverture' % (_fr(meta['index_res_mm']), _fr(R0['Y'], '%.1f')))
    item('Positions', ('horloge interne, sans encodeur' if info['internal_clock'] else 'encodeur')
         + ', pas nominal %s mm' % _fr(meta['scan_res_mm']))
    if info['rate_hz']:
        item("Fréquence d'acquisition", '%s Hz%s' % (_fr(info['rate_hz']),
             ' (valeur fixée)' if D['rate_ok'] else ' (mode %s : à vérifier)' % (info['rate_mode'] or '?')))
    item('Échantillonnage', '%s ns (%s MHz)' % (_fr(meta['sample_period_s'] * 1e9, '%.0f'),
                                                _fr(1e-6 / meta['sample_period_s'], '%.0f')))
    if info['stop_at_end'] is False:
        item('Tampon circulaire', ('coupure à la position %d (%s mm), remise en ordre automatique'
                                   % (D['wrap']['index'], _fr(D['wrap']['index'] * meta['scan_res_mm'], '%.2f')))
             if D['wrap'] else 'arrêt en fin de scan désactivé, aucune coupure détectée')
    section('Matériau et calibration')
    item('Vitesse du son (setup)', '%d m/s' % round(meta['velocity_m_s']))
    if info['cal_done']:
        item('Calibrations', ', '.join('%s %s' % (k, '✓' if v else '✗') for k, v in info['cal_done'].items()))
        item('Corrections', 'appliquées aux faisceaux (non relues par ce script)' if info['cal_applied']
             else 'aucune appliquée aux données')
    if D['thickness']:
        item('Épaisseur', '%s mm (%s)' % (_fr(D['thickness']), D['thickness_src']))
    if R0['calib']:
        item('Vitesse recalée', '%d m/s (le fond sortait à %s mm)' % (round(R0['vel']),
                                                                     _fr(R0['calib']['measured_mm'], '%.2f')))
    item('Profondeur couverte', '%s à %s mm (zéro au centre de l\'écho de surface)'
         % (_fr(R0['zmin'], '%.1f'), _fr(R0['zmax'], '%.1f')))

    r = 0
    for label, value in rows:
        if value is None:
            ttk.Label(box, text=label, font=('TkDefaultFont', 10, 'bold')).grid(
                row=r, column=0, columnspan=2, sticky='w', pady=(6 if r else 0, 2))
        else:
            ttk.Label(box, text=label, foreground='#555').grid(row=r, column=0, sticky='nw', padx=(0, 14))
            ttk.Label(box, text=value, wraplength=400).grid(row=r, column=1, sticky='w')
        r += 1

    # ------------------------------------------------ 2. a entrer (ce que le fichier ne sait pas)
    ask = ttk.LabelFrame(frm, text=' À entrer ', padding=(12, 8))
    ask.grid(row=2, column=0, sticky='ew', pady=(12, 0))
    a = 0

    def field(label, var, unit, hint):
        nonlocal a
        ttk.Label(ask, text=label).grid(row=a, column=0, sticky='w', padx=(0, 10), pady=2)
        ttk.Entry(ask, textvariable=var, width=9, justify='right').grid(row=a, column=1, sticky='w', pady=2)
        ttk.Label(ask, text=unit).grid(row=a, column=2, sticky='w', padx=(6, 0))
        a += 1
        ttk.Label(ask, text=hint, foreground='#555', wraplength=520).grid(row=a, column=0, columnspan=3,
                                                                          sticky='w', pady=(0, 6))
        a += 1

    v_speed = tk.StringVar(value=_fr(meta['scan_res_mm'] * info['rate_hz']) if info['rate_hz'] else '')
    field('Vitesse de déplacement de la sonde', v_speed, 'mm/s',
          "Le fichier ne mesure pas ta vitesse réelle. La valeur proposée est celle que suppose le TOPAZ "
          "(pas nominal × fréquence).")
    v_rate = tk.StringVar(value=_fr(info['rate_hz']) if info['rate_hz'] else '')
    if not D['rate_ok']:
        field("Fréquence d'acquisition", v_rate, 'Hz',
              "Pas fiable dans le fichier : entre la valeur réglée sur le TOPAZ (Digitizer > Acquisition Rate).")
    v_thick = tk.StringVar(value='')
    if not D['thickness'] and not info['cal_applied']:
        field('Épaisseur connue de la pièce (optionnel)', v_thick, 'mm',
              "Absente du fichier. Si tu la connais, la vitesse du son sera recalée sur l'écho de fond.")
    lbl_len = ttk.Label(ask, foreground='#2F5D7C')
    lbl_len.grid(row=a, column=0, columnspan=3, sticky='w', pady=(0, 6))
    a += 1
    v_rev = tk.BooleanVar(value=False)
    ttk.Checkbutton(ask, text="Balayage fait dans le sens inverse (le fichier ne connaît pas le sens de ta main)",
                    variable=v_rev).grid(row=a, column=0, columnspan=3, sticky='w')
    a += 1

    def upd(*_):
        try:
            step = _num(v_speed.get()) / _num(v_rate.get())
            lbl_len.config(text='Pas : %s mm par position, longueur totale : %s mm'
                           % (_fr(step, '%.3f'), _fr(step * (nx - 1), '%.1f')))
        except Exception:
            lbl_len.config(text='Entre une vitesse valide.')
    v_speed.trace_add('write', upd)
    v_rate.trace_add('write', upd)
    upd()

    # ------------------------------------------------ 3. export
    status = ttk.Label(frm, text='', foreground='#2F5D7C')
    status.grid(row=3, column=0, sticky='w', pady=(10, 0))
    btns = ttk.Frame(frm)
    btns.grid(row=4, column=0, sticky='e', pady=(8, 0))

    def go():
        try:
            params = dict(speed=_num(v_speed.get()), rate=_num(v_rate.get()), velocity=meta['velocity_m_s'],
                          thickness=D['thickness'] or _num(v_thick.get()), unwrap=bool(D['wrap']),
                          reverse=v_rev.get())
            if not params['speed'] or not params['rate'] or params['speed'] <= 0 or params['rate'] <= 0:
                raise ValueError
        except (ValueError, TypeError):
            messagebox.showerror('Paramètre invalide', 'Vérifie la vitesse de déplacement'
                                 + ('' if D['rate_ok'] else ' et la fréquence') + '.')
            return
        base = os.path.splitext(os.path.basename(src))[0]
        out = filedialog.asksaveasfilename(
            title='Enregistrer le fichier 3D', initialdir=os.path.dirname(src), initialfile=base + '_3D.html',
            defaultextension='.html',
            filetypes=[('Vue 3D interactive', '*.html'), ('Volume ParaView', '*.vtk'), ('Volume NumPy', '*.npz')])
        if not out:
            return
        status.config(text='Traitement en cours…')
        root.config(cursor='watch')
        root.update()
        try:
            R = process(vol, meta, params, D['wrap'])
            export(out, R, info)
        except Exception as e:
            root.config(cursor='')
            status.config(text='')
            messagebox.showerror('Export impossible', str(e))
            return
        root.config(cursor='')
        status.config(text='Enregistré : ' + os.path.basename(out))
        msg = 'Fichier enregistré :\n%s\n\nPas du scan : %s mm, longueur %s mm' % (
            out, _fr(R['dx'], '%.3f'), _fr(R['X'], '%.1f'))
        if out.lower().endswith('.html'):
            if messagebox.askyesno('Terminé', msg + '\n\nOuvrir la vue 3D dans le navigateur ?'):
                webbrowser.open('file://' + os.path.abspath(out))
        else:
            messagebox.showinfo('Terminé', msg)

    ttk.Button(btns, text='Fermer', command=root.destroy).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(btns, text='Choisir où enregistrer…', command=go).grid(row=0, column=1)
    root.mainloop()


# =====================================================================  page HTML
HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="fr">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>__TITLE__</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Barlow:wght@400;500;600&family=Barlow+Condensed:wght@500;600&display=swap" rel="stylesheet">
<style>
:root { --paper:#EDF1F4; --panel:#F8FAFB; --ink:#17222E; --ink-soft:#4A5A69; --steel:#2F5D7C; --rule:#C9D3DC;
  --signal:#C98A00; --signal-bg:#FFF4D6; box-sizing:border-box;
  padding-top:env(safe-area-inset-top,0px); padding-bottom:env(safe-area-inset-bottom,0px); }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --paper:#121A22; --panel:#19232D;
  --ink:#E3E9EF; --ink-soft:#9DAEBD; --steel:#86B5D4; --rule:#2C3A47; --signal:#F2C14E; --signal-bg:#2E2714; } }
:root[data-theme="dark"] { --paper:#121A22; --panel:#19232D; --ink:#E3E9EF; --ink-soft:#9DAEBD; --steel:#86B5D4;
  --rule:#2C3A47; --signal:#F2C14E; --signal-bg:#2E2714; }
*,*::before,*::after { box-sizing:inherit; }
html { scroll-padding-top:env(safe-area-inset-top,0px); }
body { margin:0; background:var(--paper); color:var(--ink); font-family:"Barlow","Segoe UI",system-ui,-apple-system,sans-serif;
  font-size:16px; line-height:1.5; }
.wrap { max-width:1320px; margin:0 auto; padding:28px 24px 48px; }
header { max-width:78ch; margin-bottom:22px; }
.file { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:500; font-size:15px; color:var(--ink-soft); margin:0 0 6px; }
h1 { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:600; font-size:clamp(30px,4.6vw,52px);
  line-height:1.04; margin:0 0 12px; }
.lede { margin:0; color:var(--ink-soft); max-width:72ch; }
.stage { display:grid; grid-template-columns:minmax(0,1fr) 300px; gap:20px; align-items:start; }
.panel { background:var(--panel); border:1px solid var(--rule); border-radius:6px; }
#view3d { height:min(64vh,620px); min-height:380px; }
.loading { display:flex; align-items:center; justify-content:center; height:100%; color:var(--ink-soft); padding:20px; text-align:center; }
aside { padding:18px 18px 8px; }
aside h2, .cut h2 { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:600; font-size:19px; margin:0 0 10px; }
.control { margin-bottom:18px; }
.control label, .control .lbl { display:flex; justify-content:space-between; align-items:baseline; font-weight:500; font-size:15px; margin-bottom:6px; }
.control output { font-variant-numeric:tabular-nums; color:var(--steel); font-weight:600; }
.hint { font-size:13.5px; color:var(--ink-soft); margin:4px 0 0; }
input[type=range] { width:100%; accent-color:var(--steel); }
.seg { display:flex; border:1px solid var(--rule); border-radius:5px; overflow:hidden; }
.seg button { flex:1; padding:7px 0; border:0; background:transparent; color:var(--ink); font:inherit; font-size:15px; cursor:pointer; }
.seg button + button { border-left:1px solid var(--rule); }
.seg button[aria-pressed="true"] { background:var(--steel); color:var(--panel); font-weight:600; }
.seg button:disabled { opacity:.4; cursor:default; }
button:focus-visible, input:focus-visible { outline:2px solid var(--signal); outline-offset:2px; }
.reading { border-top:1px solid var(--rule); padding-top:14px; margin-top:4px; }
.reading dl { margin:0; display:grid; grid-template-columns:74px 1fr; gap:8px 10px; font-size:14.5px; }
.reading dt { font-variant-numeric:tabular-nums; font-weight:600; color:var(--steel); text-align:right; }
.reading dd { margin:0; }
.reading dd.key { background:var(--signal-bg); border-left:3px solid var(--signal); padding:2px 6px; margin-left:-9px; }
.cuts { margin-top:20px; display:grid; grid-template-columns:1.15fr 1.15fr 0.7fr; gap:20px; }
.cut { padding:14px 14px 6px; min-width:0; }
.cut .plot { height:290px; }
.cut p { font-size:13.5px; color:var(--ink-soft); margin:0 0 8px; }
.caveats { margin-top:26px; max-width:78ch; color:var(--ink-soft); font-size:14.5px; }
.caveats p { margin:0 0 8px; }
@media (max-width:1080px) { .stage { grid-template-columns:1fr; } .cuts { grid-template-columns:1fr 1fr; } .cuts .cut:last-child { grid-column:span 2; } }
@media (max-width:680px) { .wrap { padding:20px 14px 36px; } .cuts { grid-template-columns:1fr; } .cuts .cut:last-child { grid-column:auto; } #view3d { height:440px; min-height:0; } }
</style>
</head>
<body>
<div class="wrap">
  <header>
    <p class="file">__FILELINE__</p>
    <h1>__H1__</h1>
    <p class="lede">__LEDE__</p>
  </header>
  <section class="stage">
    <div class="panel"><div id="view3d"><div class="loading">Chargement du volume…</div></div></div>
    <aside class="panel" aria-label="Réglages de la vue 3D">
      <h2>Réglages</h2>
      <div class="control">
        <span class="lbl" id="thr-lbl">Seuil de la surface <output id="thr-out"></output></span>
        <div class="seg" role="group" aria-labelledby="thr-lbl">
          <button type="button" data-thr="15" aria-pressed="false">15 %</button>
          <button type="button" data-thr="30" aria-pressed="true">30 %</button>
          <button type="button" data-thr="50" aria-pressed="false">50 %</button>
        </div>
        <p class="hint">Amplitude d'écran à partir de laquelle un écho devient une surface.</p>
      </div>
      <div class="control">
        <label for="zlo">Profondeur minimale <output id="zlo-out"></output></label>
        <input type="range" id="zlo" step="0.1">
      </div>
      <div class="control">
        <label for="zhi">Profondeur maximale <output id="zhi-out"></output></label>
        <input type="range" id="zhi" step="0.1">
        <p class="hint">Monte la profondeur minimale pour cacher l'écho de surface, qui masque le reste.</p>
      </div>
      <div class="reading"><h2>Lire l'image</h2><dl>__NOTES__</dl></div>
    </aside>
  </section>
  <section class="cuts">
    <div class="panel cut">
      <h2>Vue de dessus</h2>
      <p>Amplitude max dans la fenêtre de profondeur. Clique pour placer les coupes.</p>
      <div class="plot" id="top"></div>
    </div>
    <div class="panel cut">
      <h2>Coupe le long du scan</h2>
      <div class="control"><label for="cy">Index <output id="cy-out"></output></label><input type="range" id="cy"></div>
      <div class="plot" id="side"></div>
    </div>
    <div class="panel cut">
      <h2>Coupe transversale</h2>
      <div class="control"><label for="cx">Scan <output id="cx-out"></output></label><input type="range" id="cx"></div>
      <div class="plot" id="end"></div>
    </div>
  </section>
  <section class="caveats">__CAVEATS__</section>
</div>
<script id="uvdata" type="application/json">__DATA__</script>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-strict-dist-min@2.35.2/plotly-strict.min.js"></script>
<script>
(function () {
  const raw = JSON.parse(document.getElementById('uvdata').textContent);
  const M = raw.meta;
  function b64(s) { const bin = atob(s); const out = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i); return out; }
  const meshes = {};
  for (const key of Object.keys(raw.meshes)) {
    const m = raw.meshes[key];
    const vq = new Uint16Array(b64(m.v).buffer);
    const fb = b64(m.f).buffer;
    const f = m.fw === 4 ? new Uint32Array(fb) : new Uint16Array(fb);
    const n = m.nv, x = new Float32Array(n), y = new Float32Array(n), z = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = vq[3 * i] / 65535 * M.X; y[i] = vq[3 * i + 1] / 65535 * M.Y;
      z[i] = M.zmin + vq[3 * i + 2] / 65535 * (M.zmax - M.zmin);
    }
    meshes[key] = { x, y, z, f, nf: m.nf };
  }
  const S = b64(raw.slices), NX = M.nx, NY = M.ny, NS = M.ns;
  const xs = Array.from({ length: NX }, (_, i) => +(i * M.dx).toFixed(2));
  const ys = Array.from({ length: NY }, (_, i) => +(i * M.dy).toFixed(2));
  const depth = Array.from({ length: NS }, (_, k) => +((k - M.z0) * M.dz).toFixed(3));
  const at = (ix, iy, k) => S[(ix * NY + iy) * NS + k] * (100 / 255);
  const kOf = d => Math.min(NS - 1, Math.max(0, Math.round(d / M.dz + M.z0)));
  const first = ['30', '15', '50'].find(k => meshes[k]);
  const state = { thr: first, zlo: M.zlo, zhi: M.zhi, cx: M.cx, cy: M.cy };
  document.querySelectorAll('.seg button').forEach(b => { if (!meshes[b.dataset.thr]) b.disabled = true; b.setAttribute('aria-pressed', String(b.dataset.thr === first)); });

  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const depthScale = [[0, '#9CC3DC'], [0.26, '#F2B705'], [0.58, '#C2410C'], [1, '#4B2356']];
  const cfg = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['toImage', 'resetCameraLastSave3d'] };
  function baseLayout() {
    const ink = css('--ink'), rule = css('--rule');
    return { paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { family: 'Barlow, Segoe UI, sans-serif', color: ink, size: 13 },
      margin: { l: 52, r: 12, t: 8, b: 44 }, xaxis: { gridcolor: rule, zerolinecolor: rule, linecolor: rule }, yaxis: { gridcolor: rule, zerolinecolor: rule, linecolor: rule } };
  }
  let camera = { eye: { x: 0.95, y: -2.15, z: 1.2 }, center: { x: 0, y: 0, z: -0.08 }, up: { x: 0, y: 0, z: 1 } };
  function draw3d() {
    const m = meshes[state.thr], I = [], J = [], K = [], f = m.f, z = m.z;
    for (let t = 0; t < m.nf; t++) {
      const a = f[3 * t], b = f[3 * t + 1], c = f[3 * t + 2];
      if (z[a] >= state.zlo && z[b] >= state.zlo && z[c] >= state.zlo && z[a] <= state.zhi && z[b] <= state.zhi && z[c] <= state.zhi) { I.push(a); J.push(b); K.push(c); }
    }
    const ink = css('--ink'), rule = css('--rule'), el = document.getElementById('view3d'), narrow = el.clientWidth < 600;
    const ax = t => ({ title: { text: t }, gridcolor: rule, zerolinecolor: rule, showbackground: false, color: ink });
    const trace = { type: 'mesh3d', x: m.x, y: m.y, z: m.z, i: I, j: J, k: K, intensity: m.z, cmin: 0, cmax: M.cmax,
      colorscale: depthScale, showscale: !narrow,
      colorbar: { title: { text: 'Profondeur (mm)', side: 'right' }, thickness: 12, len: 0.7, outlinewidth: 0, tickfont: { color: ink } },
      lighting: { ambient: 0.62, diffuse: 0.4, specular: 0.05, roughness: 0.85, fresnel: 0.05 }, lightposition: { x: 100, y: -200, z: 300 },
      hovertemplate: 'scan %{x:.1f} mm<br>index %{y:.1f} mm<br>profondeur %{z:.1f} mm<extra></extra>' };
    const layout = { paper_bgcolor: 'rgba(0,0,0,0)', font: { family: 'Barlow, Segoe UI, sans-serif', color: ink, size: 13 }, margin: { l: 0, r: 0, t: 0, b: 0 }, uirevision: 'keep',
      scene: { xaxis: Object.assign(ax('Scan (mm)'), { range: [0, M.X] }), yaxis: Object.assign(ax('Index (mm)'), { range: [0, M.Y] }),
        zaxis: Object.assign(ax('Profondeur (mm)'), { range: [M.zmax, M.zmin] }), aspectmode: 'data', camera } };
    if (el.querySelector('.loading')) el.innerHTML = '';
    Plotly.react(el, [trace], layout, cfg);
  }
  function drawTop() {
    const k0 = kOf(state.zlo), k1 = kOf(state.zhi), zz = [];
    for (let iy = 0; iy < NY; iy++) { const row = new Array(NX);
      for (let ix = 0; ix < NX; ix++) { let mx = 0; const base = (ix * NY + iy) * NS; for (let k = k0; k <= k1; k++) { const v = S[base + k]; if (v > mx) mx = v; } row[ix] = mx * (100 / 255); }
      zz.push(row); }
    const sig = css('--signal'), L = baseLayout();
    L.xaxis.title = { text: 'Scan (mm)' }; L.yaxis.title = { text: 'Index (mm)' };
    L.shapes = [{ type: 'line', x0: state.cx, x1: state.cx, y0: 0, y1: M.Y, line: { color: sig, width: 2 } }, { type: 'line', x0: 0, x1: M.X, y0: state.cy, y1: state.cy, line: { color: sig, width: 2 } }];
    Plotly.react('top', [{ type: 'heatmap', x: xs, y: ys, z: zz, zmin: 0, zmax: 100, colorscale: 'Jet', colorbar: { title: { text: '%' }, thickness: 10, outlinewidth: 0 },
      hovertemplate: 'scan %{x} mm<br>index %{y} mm<br>%{z:.0f} %<extra></extra>' }], L, cfg);
  }
  const lines = xmax => [state.zlo, state.zhi].map(d => ({ type: 'line', x0: 0, x1: xmax, y0: d, y1: d, line: { color: css('--signal'), width: 1.5, dash: 'dot' } }));
  function drawSide() {
    const iy = Math.min(NY - 1, Math.round(state.cy / M.dy)), zz = [];
    for (let k = 0; k < NS; k++) { const row = new Array(NX); for (let ix = 0; ix < NX; ix++) row[ix] = at(ix, iy, k); zz.push(row); }
    const L = baseLayout(); L.xaxis.title = { text: 'Scan (mm)' }; L.yaxis.title = { text: 'Profondeur (mm)' }; L.yaxis.autorange = 'reversed'; L.shapes = lines(M.X);
    Plotly.react('side', [{ type: 'heatmap', x: xs, y: depth, z: zz, zmin: 0, zmax: 100, colorscale: 'Jet', showscale: false, hovertemplate: 'scan %{x} mm<br>profondeur %{y:.1f} mm<br>%{z:.0f} %<extra></extra>' }], L, cfg);
  }
  function drawEnd() {
    const ix = Math.min(NX - 1, Math.round(state.cx / M.dx)), zz = [];
    for (let k = 0; k < NS; k++) { const row = new Array(NY); for (let iy = 0; iy < NY; iy++) row[iy] = at(ix, iy, k); zz.push(row); }
    const L = baseLayout(); L.xaxis.title = { text: 'Index (mm)' }; L.yaxis.title = { text: 'Profondeur (mm)' }; L.yaxis.autorange = 'reversed'; L.shapes = lines(M.Y);
    Plotly.react('end', [{ type: 'heatmap', x: ys, y: depth, z: zz, zmin: 0, zmax: 100, colorscale: 'Jet', showscale: false, hovertemplate: 'index %{x} mm<br>profondeur %{y:.1f} mm<br>%{z:.0f} %<extra></extra>' }], L, cfg);
  }
  const fmt = v => (Math.round(v * 10) / 10).toFixed(1).replace('.', ',') + ' mm';
  const zlo = document.getElementById('zlo'), zhi = document.getElementById('zhi'), cx = document.getElementById('cx'), cy = document.getElementById('cy');
  zlo.min = zhi.min = M.zmin.toFixed(1); zlo.max = zhi.max = M.zmax.toFixed(1); zlo.value = state.zlo; zhi.value = state.zhi;
  cx.min = 0; cx.max = M.X.toFixed(2); cx.step = M.dx; cx.value = state.cx; cy.min = 0; cy.max = M.Y.toFixed(2); cy.step = M.dy; cy.value = state.cy;
  function sync() { document.getElementById('zlo-out').textContent = fmt(state.zlo); document.getElementById('zhi-out').textContent = fmt(state.zhi);
    document.getElementById('cx-out').textContent = fmt(state.cx); document.getElementById('cy-out').textContent = fmt(state.cy); document.getElementById('thr-out').textContent = state.thr + ' %'; }
  let pending = null; const later = fn => { if (pending) cancelAnimationFrame(pending); pending = requestAnimationFrame(() => { pending = null; fn(); }); };
  zlo.addEventListener('input', () => { state.zlo = Math.min(+zlo.value, state.zhi - 0.5); zlo.value = state.zlo; sync(); later(all2); });
  zhi.addEventListener('input', () => { state.zhi = Math.max(+zhi.value, state.zlo + 0.5); zhi.value = state.zhi; sync(); later(all2); });
  cx.addEventListener('input', () => { state.cx = +cx.value; sync(); later(() => { drawTop(); drawEnd(); }); });
  cy.addEventListener('input', () => { state.cy = +cy.value; sync(); later(() => { drawTop(); drawSide(); }); });
  document.querySelectorAll('.seg button').forEach(btn => btn.addEventListener('click', () => { state.thr = btn.dataset.thr;
    document.querySelectorAll('.seg button').forEach(b => b.setAttribute('aria-pressed', String(b === btn))); sync(); draw3d(); }));
  function all2() { draw3d(); drawTop(); drawSide(); drawEnd(); }
  function all() { sync(); all2(); }
  if (!window.Plotly) { document.querySelector('#view3d .loading').textContent = 'La bibliothèque de graphiques n\u2019a pas pu se charger. Vérifie ta connexion Internet et recharge la page.'; return; }
  if (document.getElementById('view3d').clientWidth < 600) camera.eye = { x: 1.35, y: -3.0, z: 1.7 };
  all();
  document.getElementById('view3d').on('plotly_relayout', ev => { if (ev && ev['scene.camera']) camera = ev['scene.camera']; });
  document.getElementById('top').on('plotly_click', ev => { const p = ev.points && ev.points[0]; if (!p) return;
    state.cx = p.x; state.cy = p.y; cx.value = state.cx; cy.value = state.cy; sync(); drawTop(); drawSide(); drawEnd(); });
  const mq = window.matchMedia('(prefers-color-scheme: dark)'); if (mq.addEventListener) mq.addEventListener('change', all);
  new MutationObserver(all).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
})();
</script>
</body>
</html>
"""


if __name__ == '__main__':
    run_gui()
