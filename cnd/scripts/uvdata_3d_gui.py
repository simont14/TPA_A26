#!/usr/bin/env python3
"""
uvdata_3d_gui.py : reconstruction 3D d'une acquisition TOPAZ / UltraVision (.UVData), sans UltraVision.

Le script detecte tout seul comment le balayage a ete fait (Encoder ID du setup) :

  Horloge interne (sonde deplacee a la main, sans encodeur)
    La position de balayage avance avec l'horloge du TOPAZ. Le fichier ne contient aucune distance : il faut
    entrer la vitesse de deplacement. Le tampon circulaire est remis en ordre s'il a ete reecrit.

  Encodeur de position
    La distance n'est pas enregistree avec chaque mesure : elle est dans l'indice de la position
    (distance = indice x ScanSamplingResolution, lu dans le fichier). Rien a entrer. Si la sonde avance de plus
    d'un pas entre deux tirs, des cases restent vides (trous) : elles sont comblees par interpolation lineaire.
    Pas de tampon circulaire. Les calibrations sont deja dans le setup.

Deroulement :
  1. une fenetre te demande le fichier de mesure (.UVData)
  2. une fenetre montre ce qui a ete lu dans le fichier (et ce qu'il reste a entrer, en mode horloge)
  3. une fenetre te demande ou enregistrer le fichier 3D :
       .html  vue 3D interactive, s'ouvre dans un navigateur   (demande scipy et scikit-image)
       .vtk   volume pour ParaView (logiciel gratuit)
       .npz   volume NumPy pour ton propre traitement

Installation :
    pip install numpy scipy scikit-image
Lancement :
    python uvdata_3d_gui.py                       (fenetres)
    python uvdata_3d_gui.py fichier.UVData        (sans fenetre, ecrit fichier_3D.html)
        options : -o sortie.(html|vtk|npz)   --inverser   --sans-interpolation (encodeur)
                  --vitesse mm/s   --frequence Hz   --epaisseur mm (horloge)

Le fichier .UVData est lu en lecture seule. Son format a ete deduit de fichiers
UltraVision Touch 3.8R11 (TOPAZ16, balayage lineaire) ; aucune protection de licence n'est touchee.
"""
import argparse
import base64
import html
import json
import math
import os
import re
import struct
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
    enc = re.search(r'<Encoder [^>]*mechanical:Encoder[^>]*>\s*<Name>Encoder 1</Name>\s*<Inverted>(\w*)</Inverted>'
                    r'.*?<Resolution>([^<]*)</Resolution>', hw, re.S)
    info['encoder_inverted'] = bool(enc) and enc.group(1) == 'true'          # sens de comptage de l'encodeur 1
    info['encoder_res_mm'] = float(enc.group(2)) * 1000 if enc else None     # mm par coup d'encodeur
    info['encoder_divider'] = _tag(hw, 'Divider', int, None)
    info['velocity_setup_m_s'] = _tag(spec, 'SoundVelocityLongitudinal', float)   # vitesse du materiau avant calibration
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
                if o > 0 and o + csize <= len(blob):                  # offset 0 : position sans mesure
                    vol[x, y] = np.frombuffer(blob, dtype=dtype, count=ns, offset=o)
        filled = (offs > 0).any(axis=0)                                # positions du scan qui ont une mesure
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
            filled=filled,
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
    dx = params['step'] if params.get('step') else params['speed'] / params['rate']   # 'step' : pas lu d'un encodeur
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
                speed=params.get('speed'), rate=params.get('rate'),
                X=(nx - 1) * dx, Y=(ny - 1) * dy, zmin=-z0 * dz, zmax=(ns - 1 - z0) * dz,
                thickness=params['thickness'])


# =====================================================================  mode encodeur
def gap_runs(filled):
    """Longueurs des series de positions vides entre la premiere et la derniere position mesuree."""
    idx = np.nonzero(filled)[0]
    if len(idx) < 2:
        return []
    return [int(g) for g in np.diff(idx) - 1 if g > 0]


def prepare_positions(vol, filled, interpolate=True):
    """Garde la zone balayee (premiere a derniere position mesuree) et comble les trous par interpolation
    lineaire entre les deux positions mesurees voisines. Sans interpolation, les trous restent a zero."""
    idx = np.nonzero(filled)[0]
    if len(idx) < 2:
        raise ValueError("Moins de deux positions contiennent une mesure : rien à reconstruire.")
    first, last = int(idx[0]), int(idx[-1])
    v = vol[first:last + 1].copy()
    f = filled[first:last + 1]
    holes = np.nonzero(~f)[0]
    if interpolate and len(holes):
        known = np.nonzero(f)[0]
        nxt = np.searchsorted(known, holes)
        lo, hi = known[nxt - 1], known[nxt]
        w = ((holes - lo) / (hi - lo)).astype(np.float32)[:, None, None]
        v[holes] = v[lo] * (1 - w) + v[hi] * w
    return dict(v=v, first=first, last=last, n_pos=len(f), n_holes=int(len(holes)), gaps=gap_runs(filled),
                interpolated=bool(interpolate and len(holes)))


def window_depth_mm(meta, ns):
    return ns * meta['sample_period_s'] * meta['velocity_m_s'] / 2 * 1000


def calibration_in_setup(info, meta):
    """Vrai si la calibration est deja dans le setup : la vitesse du fichier n'est plus celle du materiau."""
    nominal = info.get('velocity_setup_m_s')
    return bool(info['cal_applied']) or (nominal is not None and abs(meta['velocity_m_s'] - nominal) > 1.0)


def backwall_target(info, meta, v):
    """Epaisseur pour recaler la vitesse, seulement si la calibration n'est pas deja dans le setup ET si
    l'echo de fond est dans la fenetre enregistree. Retourne (epaisseur ou None, raison)."""
    if calibration_in_setup(info, meta):
        return None, 'calibration déjà dans le setup'
    thick = info['cal_thickness_mm'] or info['specimen_thickness_mm']
    if not thick:
        return None, 'épaisseur inconnue'
    dz = meta['sample_period_s'] * meta['velocity_m_s'] / 2 * 1000
    z0 = find_surface(v)
    if z0 + 1.1 * thick / dz >= v.shape[2]:
        return None, "écho de fond hors de la fenêtre enregistrée"
    kb = find_backwall(v, z0, dz, thick)
    if kb is None or abs((kb - z0) * dz / thick - 1) > 0.05:
        return None, "écho de fond introuvable près de l'épaisseur"
    return thick, 'écho de fond trouvé'


def build_encoder(info, vol, meta, interpolate=True, reverse=False):
    """Tout ce qu'il faut pour afficher et exporter une acquisition faite avec encodeur."""
    pos = prepare_positions(vol, meta['filled'], interpolate)
    thick, why = backwall_target(info, meta, pos['v'])
    params = dict(step=meta['scan_res_mm'], velocity=meta['velocity_m_s'], thickness=thick,
                  unwrap=False, reverse=reverse)
    R = process(pos['v'], meta, params, None)
    return dict(info=info, meta=meta, pos=pos, thick=thick, thick_why=why, R=R, shape=vol.shape)


def encoder_texts(B):
    """Textes de la page HTML pour une acquisition avec encodeur."""
    info, meta, pos, R = B['info'], B['meta'], B['pos'], B['R']
    nx, ny, ns = R['v'].shape
    start = pos['first'] * meta['scan_res_mm']
    lede = ("%d positions × %d faisceaux × %d échantillons, lus directement dans le fichier. Les distances le long "
            "du scan viennent de l'encodeur. Fais tourner la vue avec la souris ou le doigt, et resserre la fenêtre "
            "de profondeur pour isoler une couche." % (nx, ny, ns))
    cave = ["Pas du scan : %s mm par position, lu dans le fichier (encodeur), longueur %s mm à partir de %s mm. "
            "Les longueurs le long du scan sont mesurées par l'encodeur, pas estimées à partir d'une vitesse."
            % (_fr(R['dx'], '%.3f'), _fr(R['X'], '%.1f'), _fr(start, '%.1f'))]
    if pos['n_holes']:
        longest = max(pos['gaps'])
        cave.append("%d positions sur %d (%s %%) n'avaient aucune mesure, le plus long trou faisant %d positions "
                    "(%s mm). %s" % (pos['n_holes'], pos['n_pos'], _fr(100.0 * pos['n_holes'] / pos['n_pos'], '%.0f'),
                                     longest, _fr(longest * R['dx'], '%.2f'),
                                     "Elles ont été comblées par interpolation linéaire entre les positions voisines. "
                                     "Ce ne sont pas des mesures." if pos['interpolated']
                                     else "Elles sont laissées vides (valeur nulle) : elles apparaissent comme des "
                                          "creux dans la surface 3D."))
    if R['calib']:
        cave.append("Profondeurs recalées sur l'épaisseur de %s mm : le fond sortait à %s mm, d'où une vitesse "
                    "d'environ %d m/s." % (_fr(B['thick']), _fr(R['calib']['measured_mm']), round(R['calib']['velocity'])))
    else:
        cave.append("Vitesse du son utilisée : %d m/s%s. Le zéro est au centre de l'écho de surface."
                    % (round(R['vel']), ", valeur du setup après calibration" if calibration_in_setup(info, meta)
                       else ''))
    ep = info['cal_thickness_mm'] or info['specimen_thickness_mm']
    if ep and B['thick_why'] == "écho de fond hors de la fenêtre enregistrée":
        cave.append("La fenêtre enregistrée couvre %s mm de profondeur : l'écho de fond du bloc (%s mm) n'y est pas."
                    % (_fr(window_depth_mm(meta, ns), '%.1f'), _fr(ep)))
    cave += R['notes']
    return dict(lede=lede, caveats=cave)


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


HTML_LEVELS = (10, 15, 20, 30, 40, 50, 60, 70, 80)      # seuils d'amplitude (%) proposes dans la vue 3D


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
    for lvl in HTML_LEVELS:
        if not (g.min() < lvl < g.max()):
            continue
        verts, faces, _, _ = measure.marching_cubes(g, level=lvl, spacing=(R['dx'] * sx, R['dy'], R['dz'] * sz))
        dep = verts[:, 2] - R['z0'] * R['dz']
        vb = np.stack([q(verts[:, 0], 0, R['X']), q(verts[:, 1], 0, R['Y']), q(dep, R['zmin'], R['zmax'])], 1)
        fw = 2 if len(verts) < 65536 else 4
        fb = faces.astype('<u2' if fw == 2 else '<u4')
        meshes[str(lvl)] = dict(v=_b64(vb), f=_b64(fb), nv=len(verts), nf=len(faces), fw=fw)
    if not meshes:
        raise RuntimeError('Aucune surface aux seuils proposés : le signal est peut-être vide.')

    n_s = min(nx, 250)                                    # coupes : au plus 250 positions, etendue du scan conservee
    sel = np.round(np.linspace(0, nx - 1, n_s)).astype(int)
    s8 = np.round(v[sel].clip(0, 100) / 100 * 255).astype(np.uint8)
    dx_s = R['X'] / (n_s - 1) if n_s > 1 else R['dx']
    thick = R['thickness']
    M = dict(nx=n_s, ny=ny, ns=ns, dx=dx_s, dy=R['dy'], dz=R['dz'], z0=R['z0'],
             X=R['X'], Y=R['Y'], zmin=R['zmin'], zmax=R['zmax'], cmax=thick or R['zmax'],
             zlo=3.0 if (thick or R['zmax']) > 8 else 0.0, zhi=R['zmax'],
             cx=round(R['X'] / 3, 1), cy=round(R['Y'] / 4, 1), name=os.path.splitext(info['file'])[0])
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
        if not info['internal_clock']:                     # encodeur : autre fenetre, rien a entrer
            run_gui_encoder(root, src, info, vol, meta)
            return
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


def run_gui_encoder(root, src, info, vol, meta):
    """Fenetre pour une acquisition faite avec encodeur : tout est lu dans le fichier, deux options a cocher."""
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    B0 = build_encoder(info, vol, meta, True, False)
    pos, R0 = B0['pos'], B0['R']
    nx, ny, ns = vol.shape

    root.deiconify()
    root.resizable(False, False)
    frm = ttk.Frame(root, padding=16)
    frm.grid(sticky='nsew')
    ttk.Label(frm, text=info['file'], font=('TkDefaultFont', 12, 'bold')).grid(row=0, column=0, sticky='w')
    ttk.Label(frm, text="Balayage avec encodeur détecté", foreground='#2F5D7C').grid(row=1, column=0, sticky='w')

    box = ttk.LabelFrame(frm, text=' Lu dans le fichier ', padding=(12, 8))
    box.grid(row=2, column=0, sticky='ew', pady=(8, 0))
    rows = []
    section = lambda t: rows.append((t, None))
    item = lambda a, b: rows.append((a, b))

    section('Acquisition')
    item('Appareil, date', ', '.join(x for x in (info['device'], info['date']) if x) or '?')
    if info['setup']:
        item('Setup', info['setup'])
    item('Grille enregistrée', '%d positions × %d faisceaux × %d échantillons' % (nx, ny, ns))
    item('Faisceaux', '%s, onde %s' % (meta['beam'] or '?', {'Longitudinal': 'longitudinale',
         'Transversal': 'transversale'}.get(meta['wave_type'], meta['wave_type'].lower())))
    item('Index', '%s mm entre faisceaux, %s mm de couverture' % (_fr(meta['index_res_mm']), _fr(R0['Y'], '%.1f')))
    item('Échantillonnage', '%s ns (%s MHz)' % (_fr(meta['sample_period_s'] * 1e9, '%.0f'),
                                                _fr(1e-6 / meta['sample_period_s'], '%.0f')))
    section('Positions (encodeur)')
    enc = 'encodeur'
    if info['encoder_res_mm']:
        enc += ', %s µm par coup' % _fr(info['encoder_res_mm'] * 1000, '%.1f')
        if info['encoder_divider']:
            enc += ' (diviseur %d)' % info['encoder_divider']
    item('Source', enc)
    item('Pas du scan', '%s mm par position (lu dans le fichier)' % _fr(meta['scan_res_mm']))
    item('Zone balayée', 'de %s à %s mm, longueur %s mm (%d positions)'
         % (_fr(pos['first'] * meta['scan_res_mm'], '%.1f'), _fr(pos['last'] * meta['scan_res_mm'], '%.1f'),
            _fr(R0['X'], '%.1f'), pos['n_pos']))
    if pos['n_holes']:
        longest = max(pos['gaps'])
        item('Trous', '%d positions sur %d sans mesure (%s %%), plus long trou : %d (%s mm)' % (
            pos['n_holes'], pos['n_pos'], _fr(100.0 * pos['n_holes'] / pos['n_pos'], '%.0f'), longest,
            _fr(longest * meta['scan_res_mm'], '%.2f')))
        if info['rate_hz']:
            item('Vitesse sans trou', 'au plus %s mm/s (pas × %s Hz)'
                 % (_fr(meta['scan_res_mm'] * info['rate_hz'], '%.1f'), _fr(info['rate_hz'])))
    else:
        item('Trous', 'aucun, toutes les positions ont une mesure')
    section('Matériau et calibration')
    in_setup = calibration_in_setup(info, meta)
    item('Vitesse du son', '%d m/s%s' % (round(meta['velocity_m_s']), (' (calibrée dans le setup, matériau %d m/s)'
         % round(info['velocity_setup_m_s'])) if in_setup and info['velocity_setup_m_s'] else ''))
    if info['cal_done']:
        item('Calibrations', ', '.join('%s %s' % (k, '✓' if v else '✗') for k, v in info['cal_done'].items()))
    ep = info['cal_thickness_mm'] or info['specimen_thickness_mm']
    if B0['thick']:
        item('Épaisseur', '%s mm (écho de fond trouvé, vitesse recalée)' % _fr(B0['thick']))
    elif ep and B0['thick_why'] == "écho de fond hors de la fenêtre enregistrée":
        item('Fenêtre', "%s mm de profondeur, sans l'écho de fond (%s mm)"
             % (_fr(window_depth_mm(meta, ns), '%.1f'), _fr(ep)))
    item('Profondeur couverte', "%s à %s mm (zéro au centre de l'écho de surface)"
         % (_fr(R0['zmin'], '%.1f'), _fr(R0['zmax'], '%.1f')))

    r = 0
    for label, value in rows:
        if value is None:
            ttk.Label(box, text=label, font=('TkDefaultFont', 10, 'bold')).grid(
                row=r, column=0, columnspan=2, sticky='w', pady=(6 if r else 0, 2))
        else:
            ttk.Label(box, text=label, foreground='#555').grid(row=r, column=0, sticky='nw', padx=(0, 14))
            ttk.Label(box, text=value, wraplength=420).grid(row=r, column=1, sticky='w')
        r += 1

    ask = ttk.LabelFrame(frm, text=' Options ', padding=(12, 8))
    ask.grid(row=3, column=0, sticky='ew', pady=(12, 0))
    ttk.Label(ask, text="Rien à entrer : les distances viennent de l'encodeur.", foreground='#555').grid(
        row=0, column=0, sticky='w', pady=(0, 6))
    v_interp = tk.BooleanVar(value=True)
    cb = ttk.Checkbutton(ask, text="Combler les trous par interpolation linéaire (sinon ils restent vides)",
                         variable=v_interp)
    cb.grid(row=1, column=0, sticky='w')
    if not pos['n_holes']:
        cb.state(['disabled'])
    v_rev = tk.BooleanVar(value=False)
    ttk.Checkbutton(ask, text="Inverser le sens du scan (si la vue est à l'envers par rapport au bloc)",
                    variable=v_rev).grid(row=2, column=0, sticky='w')

    status = ttk.Label(frm, text='', foreground='#2F5D7C')
    status.grid(row=4, column=0, sticky='w', pady=(10, 0))
    btns = ttk.Frame(frm)
    btns.grid(row=5, column=0, sticky='e', pady=(8, 0))

    def go():
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
            B = build_encoder(info, vol, meta, v_interp.get(), v_rev.get())
            export(out, B['R'], info, encoder_texts(B))
        except Exception as e:
            root.config(cursor='')
            status.config(text='')
            messagebox.showerror('Export impossible', str(e))
            return
        root.config(cursor='')
        status.config(text='Enregistré : ' + os.path.basename(out))
        msg = 'Fichier enregistré :\n%s\n\nPas du scan : %s mm, longueur %s mm' % (
            out, _fr(B['R']['dx'], '%.3f'), _fr(B['R']['X'], '%.1f'))
        if out.lower().endswith('.html'):
            if messagebox.askyesno('Terminé', msg + '\n\nOuvrir la vue 3D dans le navigateur ?'):
                webbrowser.open('file://' + os.path.abspath(out))
        else:
            messagebox.showinfo('Terminé', msg)

    ttk.Button(btns, text='Fermer', command=root.destroy).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(btns, text='Choisir où enregistrer…', command=go).grid(row=0, column=1)
    root.mainloop()


# =====================================================================  ligne de commande
def run_cli(a):
    info, vols = read_acquisition(a.fichier)
    vol, meta = vols[0]
    out = a.o or os.path.splitext(a.fichier)[0] + '_3D.html'
    if not info['internal_clock']:
        B = build_encoder(info, vol, meta, not a.sans_interpolation, a.inverser)
        export(out, B['R'], info, encoder_texts(B))
        R, pos = B['R'], B['pos']
        print('Mode encodeur : pas %.3f mm lu dans le fichier, longueur %.1f mm' % (R['dx'], R['X']))
        print('%d positions dont %d trous%s' % (pos['n_pos'], pos['n_holes'],
                                               ' (comblés)' if pos['interpolated'] else ''))
    else:
        D = file_defaults(info, vol, meta)
        rate = a.frequence or info['rate_hz']
        if not rate:
            raise SystemExit("Fréquence d'acquisition introuvable dans le fichier : ajoute --frequence Hz.")
        speed = a.vitesse or meta['scan_res_mm'] * rate
        R = process(vol, meta, dict(speed=speed, rate=rate, velocity=meta['velocity_m_s'],
                                    thickness=a.epaisseur or D['thickness'], unwrap=bool(D['wrap']),
                                    reverse=a.inverser), D['wrap'])
        export(out, R, info)
        print('Mode horloge interne : vitesse %g mm/s à %g Hz, pas %.3f mm, longueur %.1f mm'
              % (speed, rate, R['dx'], R['X']))
    print('Vitesse du son %d m/s, profondeur %.1f à %.1f mm' % (round(R['vel']), R['zmin'], R['zmax']))
    for n in R['notes']:
        print('Note :', n)
    print('Enregistré :', out)


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
  --signal:#C98A00; --signal-bg:#FFF4D6; --amber:#d97706; --teal:#0d9488; --rose:#e11d48; --blue:#2563eb;
  box-sizing:border-box; padding-top:env(safe-area-inset-top,0px); padding-bottom:env(safe-area-inset-bottom,0px); }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --paper:#121A22; --panel:#19232D;
  --ink:#E3E9EF; --ink-soft:#9DAEBD; --steel:#86B5D4; --rule:#2C3A47; --signal:#F2C14E; --signal-bg:#2E2714;
  --amber:#f59e0b; --teal:#2dd4bf; --rose:#fb7185; --blue:#60a5fa; } }
:root[data-theme="dark"] { --paper:#121A22; --panel:#19232D; --ink:#E3E9EF; --ink-soft:#9DAEBD; --steel:#86B5D4;
  --rule:#2C3A47; --signal:#F2C14E; --signal-bg:#2E2714; --amber:#f59e0b; --teal:#2dd4bf; --rose:#fb7185; --blue:#60a5fa; }
*,*::before,*::after { box-sizing:inherit; }
html { scroll-padding-top:env(safe-area-inset-top,0px); }
body { margin:0; background:var(--paper); color:var(--ink); font-family:"Barlow","Segoe UI",system-ui,-apple-system,sans-serif;
  font-size:16px; line-height:1.5; }
.wrap { max-width:1500px; margin:0 auto; padding:24px 20px 48px; }
header { max-width:78ch; margin-bottom:18px; }
.file { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:500; font-size:15px; color:var(--ink-soft); margin:0 0 6px; }
h1 { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:600; font-size:clamp(28px,4.2vw,46px); line-height:1.04; margin:0 0 10px; }
.lede { margin:0; color:var(--ink-soft); max-width:72ch; }
h2 { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:600; font-size:19px; margin:0; }
.panel { background:var(--panel); border:1px solid var(--rule); border-radius:6px; }
.stage { display:grid; grid-template-columns:minmax(0,1fr) 350px; gap:16px; align-items:start; }
.card-head { display:flex; align-items:center; flex-wrap:wrap; gap:8px 12px; padding:9px 14px; border-bottom:1px solid var(--rule); }
.card-head .sp { flex:1; }
.hint { font-size:13.5px; color:var(--ink-soft); margin:4px 0 0; }
button, select, input { font:inherit; color:inherit; }
.btn { border:1px solid var(--rule); background:transparent; border-radius:5px; padding:4px 10px; font-size:14px; cursor:pointer; color:var(--ink); }
.btn:hover { border-color:var(--steel); }
.btn.pri { background:var(--steel); color:var(--panel); border-color:var(--steel); font-weight:600; }
.btn:disabled { opacity:.45; cursor:default; }
button:focus-visible, input:focus-visible, select:focus-visible, .th:focus-visible { outline:2px solid var(--signal); outline-offset:2px; }
.seg { display:flex; border:1px solid var(--rule); border-radius:5px; overflow:hidden; }
.seg button { flex:1; padding:5px 8px; border:0; background:transparent; font-size:14px; cursor:pointer; white-space:nowrap; }
.seg button + button { border-left:1px solid var(--rule); }
.seg button[aria-pressed="true"] { background:var(--steel); color:var(--panel); font-weight:600; }
#view3d { height:min(68vh,660px); min-height:400px; }
.loading { display:flex; align-items:center; justify-content:center; height:100%; color:var(--ink-soft); padding:20px; text-align:center; }
aside { padding:6px 16px 10px; max-height:calc(68vh + 52px); overflow-y:auto; }
details.grp { border-top:1px solid var(--rule); padding:8px 0; }
details.grp:first-child { border-top:0; }
details.grp > summary { cursor:pointer; font-weight:600; font-size:15.5px; padding:2px 0; }
.row { display:flex; align-items:center; justify-content:space-between; gap:10px; margin:7px 0; font-size:14.5px; }
.row > label, .row > .lbl { flex:0 0 auto; }
.row input[type=range] { flex:1; min-width:70px; accent-color:var(--steel); }
.row select, .row input.num, .row input[type=color] { border:1px solid var(--rule); background:var(--paper); border-radius:5px; padding:3px 6px; }
.row input.num { width:84px; text-align:right; }
.pair { display:flex; align-items:center; gap:6px; } .pair input.num { width:66px; }
.rb { position:absolute; border:1.5px dashed var(--steel); background:rgba(47,93,124,.16); pointer-events:none; z-index:5; }
#plans.crop .cell, #plans.crop .nsewdrag { cursor:crosshair !important; }
.row input[type=color] { padding:0; width:40px; height:26px; }
.row output { font-variant-numeric:tabular-nums; color:var(--steel); font-weight:600; min-width:44px; text-align:right; }
.row.chk label { display:flex; align-items:center; gap:8px; cursor:pointer; }
.btns { display:flex; flex-wrap:wrap; gap:6px; margin:6px 0; }
.btns .btn { flex:1 1 auto; }
.reading { border-top:1px solid var(--rule); padding-top:12px; margin-top:6px; }
.reading dl { margin:8px 0 0; display:grid; grid-template-columns:70px 1fr; gap:6px 10px; font-size:14px; }
.reading dt { font-variant-numeric:tabular-nums; font-weight:600; color:var(--steel); text-align:right; }
.reading dd { margin:0; }
.reading dd.key { background:var(--signal-bg); border-left:3px solid var(--signal); padding:2px 6px; margin-left:-9px; }
.plans { margin-top:16px; }
#plans { padding:8px 10px 12px; overflow-x:auto; }
.pgrid { display:grid; margin:0 auto; width:max-content; }
.cell { position:relative; min-width:0; min-height:0; }
.cell .cap { position:absolute; left:58px; top:4px; font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:600; font-size:15px; pointer-events:none; }
.cell .cap i { font-style:normal; font-weight:500; color:var(--ink-soft); margin-left:8px; font-size:13px; }
.cell.top .cap { color:var(--ink); }
.slbl { position:absolute; font-size:11.5px; line-height:1.15; color:var(--ink-soft); text-align:center; font-variant-numeric:tabular-nums; pointer-events:none; }
.slbl b { display:block; color:var(--ink); font-weight:600; }
.corner { padding:8px 10px; font-size:14px; }
.corner .big { font-family:"Barlow Condensed","Arial Narrow",sans-serif; font-weight:600; font-size:30px; line-height:1.1; color:var(--steel); }
.corner .sub { color:var(--ink-soft); font-size:13px; }
.legend { margin-top:10px; font-size:12.5px; color:var(--ink-soft); }
.legend .bar { height:10px; border-radius:3px; margin:4px 0 2px; border:1px solid var(--rule); }
.legend .ends { display:flex; justify-content:space-between; font-variant-numeric:tabular-nums; }
.sl { position:absolute; touch-action:none; user-select:none; -webkit-user-select:none; cursor:pointer; }
.sl.h { height:26px; } .sl.v { width:26px; }
.sl .tr { position:absolute; background:var(--rule); border-radius:3px; }
.sl.h .tr { left:0; right:0; top:11px; height:4px; }
.sl.v .tr { top:0; bottom:0; left:11px; width:4px; }
.sl .fill { position:absolute; background:var(--c); opacity:.55; border-radius:3px; }
.sl .th { position:absolute; width:16px; height:16px; border-radius:50%; background:var(--c); border:2px solid var(--panel);
  box-shadow:0 0 0 1px var(--c); transform:translate(-50%,-50%); cursor:grab; }
.mrow { display:grid; grid-template-columns:22px 1fr; gap:2px 8px; font-size:14px; margin:3px 0; font-variant-numeric:tabular-nums; }
.dotA, .dotB { font-weight:600; } .dotA { color:var(--rose); } .dotB { color:var(--blue); }
.caveats { margin-top:24px; max-width:78ch; color:var(--ink-soft); font-size:14.5px; }
.caveats p { margin:0 0 8px; }
.toast { position:fixed; left:50%; bottom:22px; transform:translateX(-50%); background:var(--ink); color:var(--paper); padding:8px 14px;
  border-radius:6px; font-size:14px; opacity:0; transition:opacity .2s; pointer-events:none; z-index:50; }
.toast.on { opacity:.94; }
@media (max-width:1100px) { .stage { grid-template-columns:1fr; } aside { max-height:none; } }
@media (max-width:680px) { .wrap { padding:18px 12px 36px; } #view3d { height:440px; min-height:0; } }
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
    <div class="panel">
      <div class="card-head">
        <h2>Vue 3D</h2>
        <span class="hint" style="margin:0">Glisser pour tourner, molette pour zoomer, clic pour placer les coupes ou un point de mesure.</span>
        <span class="sp"></span>
        <button type="button" class="btn" data-save="three">Enregistrer l'image</button>
      </div>
      <div id="view3d"><div class="loading">Chargement du volume…</div></div>
    </div>

    <aside class="panel" id="controls" aria-label="Réglages">
      <details class="grp" open>
        <summary>Orientation</summary>
        <div class="row chk"><label><input type="checkbox" id="flipX"> Inverser le sens du scan</label></div>
        <div class="row chk"><label><input type="checkbox" id="flipY"> Inverser l'index (faisceaux)</label></div>
        <div class="btns" id="camBtns">
          <button type="button" class="btn" data-cam="iso">Isométrique</button>
          <button type="button" class="btn" data-cam="top">Dessus</button>
          <button type="button" class="btn" data-cam="front">Face</button>
          <button type="button" class="btn" data-cam="side">Côté</button>
        </div>
        <div class="row"><label for="proj">Projection</label>
          <select id="proj"><option value="perspective">Perspective</option><option value="orthographic">Orthographique</option></select></div>
      </details>

      <details class="grp" open>
        <summary>Région d'intérêt (rognage)</summary>
        <p class="hint">Limite le modèle 3D à une boîte, et les coupes si la case est cochée. En mode Rogner, trace un rectangle dans une des coupes.</p>
        <div class="row"><label for="rX0">Scan (mm)</label><span class="pair"><input type="text" inputmode="decimal" class="num" id="rX0" aria-label="Scan min"> à <input type="text" inputmode="decimal" class="num" id="rX1" aria-label="Scan max"></span></div>
        <div class="row"><label for="rY0">Index (mm)</label><span class="pair"><input type="text" inputmode="decimal" class="num" id="rY0" aria-label="Index min"> à <input type="text" inputmode="decimal" class="num" id="rY1" aria-label="Index max"></span></div>
        <div class="row"><label for="nZlo">Profondeur (mm)</label><span class="pair"><input type="text" inputmode="decimal" class="num" id="nZlo" aria-label="Profondeur min"> à <input type="text" inputmode="decimal" class="num" id="nZhi" aria-label="Profondeur max"></span></div>
        <p class="hint">Monte la profondeur minimale pour cacher l'écho de surface, qui masque le reste.</p>
        <div class="row chk"><label><input type="checkbox" id="cropPlans"> Recadrer les coupes sur la région</label></div>
        <div class="btns"><button type="button" class="btn" id="roiFit">Ajuster à la surface visible</button><button type="button" class="btn" id="roiAll">Tout afficher</button></div>
      </details>

      <details class="grp" open>
        <summary>Surface 3D</summary>
        <div class="row"><label for="thr">Seuil de la surface</label><input type="range" id="thr" min="0" max="0" step="1"><output id="thr-out"></output></div>
        <p class="hint">Amplitude d'écran à partir de laquelle un écho devient une surface.</p>
        <div class="row"><label for="colorMode">Couleur</label>
          <select id="colorMode"><option value="depth">Selon la profondeur</option><option value="solid">Couleur unie</option></select></div>
        <div class="row" id="rowPal"><label for="palette">Palette</label><select id="palette"></select></div>
        <div class="row" id="rowSolid"><label for="solid">Couleur de la surface</label><input type="color" id="solid" value="#c2410c"></div>
        <div class="row"><label for="opacity">Opacité</label><input type="range" id="opacity" min="0.1" max="1" step="0.05"><output id="opacity-out"></output></div>
        <div class="row"><label for="light">Éclairage</label>
          <select id="light"><option value="mat">Mat</option><option value="brillant">Brillant</option><option value="plat">Plat (sans relief)</option></select></div>
        <div class="row chk"><label><input type="checkbox" id="flat"> Facettes visibles</label></div>
        <div class="row chk"><label><input type="checkbox" id="colorbar"> Barre de couleur</label></div>
      </details>

      <details class="grp" open>
        <summary>Cadre et axes de la vue 3D</summary>
        <div class="row"><label for="bg">Fond</label>
          <select id="bg"><option value="theme">Selon le thème</option><option value="blanc">Blanc</option><option value="gris">Gris clair</option><option value="noir">Noir</option></select></div>
        <div class="row chk"><label><input type="checkbox" id="axes"> Axes et grille</label></div>
        <div class="row chk"><label><input type="checkbox" id="box"> Cadre du volume</label></div>
        <div class="row chk"><label><input type="checkbox" id="cuts3d"> Plans de coupe dans la vue 3D</label></div>
        <div class="row"><label for="zstretch">Étirement de la profondeur</label><input type="range" id="zstretch" min="0.5" max="6" step="0.1"><output id="zstretch-out"></output></div>
      </details>

      <details class="grp">
        <summary>Coupes</summary>
        <div class="row"><label for="ampPalette">Palette d'amplitude</label><select id="ampPalette"></select></div>
        <div class="row"><label for="ampMin">Amplitude min</label><input type="range" id="ampMin" min="0" max="100" step="1"><output id="ampMin-out"></output></div>
        <div class="row"><label for="ampMax">Amplitude max</label><input type="range" id="ampMax" min="0" max="100" step="1"><output id="ampMax-out"></output></div>
        <div class="row chk"><label><input type="checkbox" id="smooth"> Lissage de l'image</label></div>
        <div class="row chk"><label><input type="checkbox" id="grid"> Grille</label></div>
        <div class="row chk"><label><input type="checkbox" id="cursors"> Repères de coupe</label></div>
        <div class="row chk"><label><input type="checkbox" id="trueScale"> Échelle réelle (mm égaux sur les axes)</label></div>
        <div class="row"><label for="nScan">Position en scan (mm)</label><input type="text" inputmode="decimal" class="num" id="nScan"></div>
        <div class="row"><label for="nIndex">Position en index (mm)</label><input type="text" inputmode="decimal" class="num" id="nIndex"></div>
      </details>

      <details class="grp" open>
        <summary>Mesure de distance</summary>
        <div class="seg" role="group" aria-label="Action du clic">
          <button type="button" data-mode="cursor" aria-pressed="true">Placer les coupes</button>
          <button type="button" data-mode="measure" aria-pressed="false">Mesurer</button>
          <button type="button" data-mode="crop" aria-pressed="false">Rogner</button>
        </div>
        <p class="hint">En mode Mesurer, clique deux points dans n'importe quelle vue. Un troisième clic recommence. En mode Rogner, trace un rectangle dans une coupe.</p>
        <div class="mrow"><span class="dotA">A</span><span class="rd-a">pas encore placé</span></div>
        <div class="mrow"><span class="dotB">B</span><span class="rd-b">pas encore placé</span></div>
        <div class="reading" style="margin-top:6px;padding-top:8px"><div class="rd-main" style="font-weight:600;font-size:17px;color:var(--steel)"></div><div class="rd-sub hint"></div></div>
        <div class="btns"><button type="button" class="btn" id="mClear">Effacer</button><button type="button" class="btn" id="mCopy">Copier le résultat</button></div>
      </details>

      <details class="grp" open>
        <summary>Enregistrer les figures</summary>
        <div class="row"><label for="fmt">Format</label>
          <select id="fmt"><option value="png">PNG</option><option value="jpeg">JPEG</option><option value="webp">WebP</option><option value="svg">SVG (coupes seulement)</option></select></div>
        <div class="row"><label for="scale">Résolution</label>
          <select id="scale"><option value="1">×1 (écran)</option><option value="2">×2</option><option value="3">×3</option><option value="4">×4 (impression)</option></select></div>
        <div class="row"><label for="exBg">Fond de l'image</label>
          <select id="exBg"><option value="white">Blanc, texte noir</option><option value="transparent">Transparent, texte noir</option></select></div>
        <div class="row chk"><label><input type="checkbox" id="exCursors"> Garder les repères de coupe</label></div>
        <div class="row chk"><label><input type="checkbox" id="exMeasure"> Garder la mesure</label></div>
        <div class="btns">
          <button type="button" class="btn pri" data-save="sheet">Les 3 coupes en une planche</button>
          <button type="button" class="btn" data-save="three">Vue 3D</button>
          <button type="button" class="btn" data-save="top">Vue de dessus</button>
          <button type="button" class="btn" data-save="side">Coupe longitudinale</button>
          <button type="button" class="btn" data-save="end">Coupe transversale</button>
          <button type="button" class="btn" data-save="all">Tout enregistrer</button>
        </div>
        <p class="hint">Les images sont sans titre, avec les axes en mm. La légende se met dans le rapport.</p>
      </details>

      <details class="grp">
        <summary>Lire l'image</summary>
        <div class="reading" style="border-top:0;margin-top:0;padding-top:0"><dl>__NOTES__</dl></div>
      </details>
      <div class="btns" style="margin-top:8px"><button type="button" class="btn" id="reset">Réinitialiser les réglages</button></div>
    </aside>
  </section>

  <section class="panel plans">
    <div class="card-head">
      <h2>Plans de coupe</h2>
      <span class="hint" style="margin:0">Dépliés comme un dessin technique : dessus en haut, coupes longitudinale et transversale en dessous.</span>
      <span class="sp"></span>
      <div class="seg" role="group" aria-label="Action du clic" id="modeSeg2">
        <button type="button" data-mode="cursor" aria-pressed="true">Placer les coupes</button>
        <button type="button" data-mode="measure" aria-pressed="false">Mesurer</button>
        <button type="button" data-mode="crop" aria-pressed="false">Rogner</button>
      </div>
      <button type="button" class="btn" data-save="sheet">Enregistrer la planche</button>
    </div>
    <div id="plans"></div>
  </section>
  <section class="caveats">__CAVEATS__</section>
</div>
<div class="toast" id="toast" role="status"></div>
<script id="uvdata" type="application/json">__DATA__</script>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
<script>
(function () {
  'use strict';
  const raw = JSON.parse(document.getElementById('uvdata').textContent);
  const M = raw.meta;
  const $ = id => document.getElementById(id);
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
  const fmt = (v, d) => (+v).toFixed(d === undefined ? 1 : d).replace('.', ',');
  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  function b64(s) { const bin = atob(s); const out = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i); return out; }

  // ---------------------------------------------------------------- données
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
    meshes[key] = { x, y, z, f, nf: m.nf, nv: n };
  }
  const levelKeys = Object.keys(meshes).sort((a, b) => a - b);
  const S = b64(raw.slices), NX = M.nx, NY = M.ny, NS = M.ns;
  const X = M.X, Y = M.Y, Z0 = M.zmin, Z1 = M.zmax;
  const xs = Array.from({ length: NX }, (_, i) => +(i * M.dx).toFixed(4));
  const ys = Array.from({ length: NY }, (_, i) => +(i * M.dy).toFixed(4));
  const depth = Array.from({ length: NS }, (_, k) => +((k - M.z0) * M.dz).toFixed(4));
  const amp = (ix, iy, k) => S[(ix * NY + iy) * NS + k] * (100 / 255);
  const kOf = d => clamp(Math.round(d / M.dz + M.z0), 0, NS - 1);

  // ---------------------------------------------------------------- palettes
  const PAL = {
    defaut: [[0, '#9CC3DC'], [0.26, '#F2B705'], [0.58, '#C2410C'], [1, '#4B2356']],
    viridis: [[0, '#440154'], [0.25, '#3b528b'], [0.5, '#21918c'], [0.75, '#5ec962'], [1, '#fde725']],
    plasma: [[0, '#0d0887'], [0.25, '#7e03a8'], [0.5, '#cc4778'], [0.75, '#f89540'], [1, '#f0f921']],
    cividis: [[0, '#00204d'], [0.25, '#414d6b'], [0.5, '#7c7b78'], [0.75, '#bcaf6f'], [1, '#ffe945']],
    chaud: [[0, '#000000'], [0.35, '#cc0000'], [0.7, '#ffaa00'], [1, '#ffffff']],
    arc: [[0, '#00007f'], [0.125, '#0000ff'], [0.375, '#00ffff'], [0.625, '#ffff00'], [0.875, '#ff0000'], [1, '#7f0000']],
    gris: [[0, '#f2f2f2'], [1, '#222222']],
    corrosion: [[0, '#ff0000'], [0.0118, '#ff4300'], [0.0235, '#ff6f00'], [0.0392, '#ff9200'], [0.0745, '#ffbf00'], [0.1255, '#fee501'],
      [0.1765, '#f8fa08'], [0.2627, '#bffa21'], [0.3922, '#63f15d'], [0.4941, '#26dc8d'], [0.5843, '#05c3b1'], [0.651, '#00adc6'],
      [0.8314, '#006cf1'], [0.9529, '#0036fd'], [0.9843, '#0019ff'], [1, '#0000ff']]          // EVDNT_Corrosion_5.pal (Omniscan)
  };
  const PAL_NAMES = { defaut: 'Défaut (bleu à violet)', viridis: 'Viridis', plasma: 'Plasma', cividis: 'Cividis', chaud: 'Chaud', arc: 'Arc-en-ciel', gris: 'Gris', corrosion: 'Corrosion EVDNT (Omniscan)' };
  const AMP = {
    Jet: PAL.arc,
    Viridis: PAL.viridis,
    Cividis: PAL.cividis,
    Chaud: PAL.chaud,
    Gris: [[0, '#000000'], [1, '#ffffff']],
    Gris_inv: [[0, '#ffffff'], [1, '#000000']],
    Corrosion: PAL.corrosion,
    Corrosion_inv: PAL.corrosion.map(s => [+(1 - s[0]).toFixed(4), s[1]]).reverse()
  };
  const AMP_NAMES = { Jet: 'Arc-en-ciel (Jet)', Viridis: 'Viridis', Cividis: 'Cividis', Chaud: 'Chaud', Gris: 'Gris (noir à blanc)', Gris_inv: 'Gris inversé', Corrosion: 'Corrosion EVDNT (rouge à bleu)', Corrosion_inv: 'Corrosion EVDNT inversée (bleu à rouge)' };
  const LIGHT = {
    mat: { ambient: 0.62, diffuse: 0.4, specular: 0.05, roughness: 0.85, fresnel: 0.05 },
    brillant: { ambient: 0.4, diffuse: 0.85, specular: 0.9, roughness: 0.2, fresnel: 0.3 },
    plat: { ambient: 0.97, diffuse: 0.05, specular: 0, roughness: 1, fresnel: 0 }
  };
  const CAMS = {
    iso: { eye: { x: 0.5, y: -1.2, z: 0.65 }, up: { x: 0, y: 0, z: 1 } },
    top: { eye: { x: 0, y: 0.0001, z: 1.7 }, up: { x: 0, y: 1, z: 0 } },
    front: { eye: { x: 0, y: -1.7, z: 0.0001 }, up: { x: 0, y: 0, z: 1 } },
    side: { eye: { x: 1.7, y: 0, z: 0.0001 }, up: { x: 0, y: 0, z: 1 } }
  };
  const C_SIDE = '#d97706', C_END = '#0d9488', C_A = '#e11d48', C_B = '#2563eb', C_WIN = '#8b5cf6';
  const FONT = 'Barlow, Segoe UI, sans-serif';

  // ---------------------------------------------------------------- état
  const first = ['30', '15', '50'].find(k => meshes[k]) || levelKeys[0];
  const DEF = {
    thr: first, zlo: M.zlo, zhi: M.zhi, cx: M.cx, cy: M.cy, flipX: false, flipY: false, mode: 'cursor',
    colorMode: 'depth', palette: 'defaut', solid: '#c2410c', opacity: 1, light: 'mat', flat: false, colorbar: true,
    bg: 'theme', axes: true, box: true, cuts3d: true, zstretch: 1, projection: 'perspective',
    ampPalette: 'Jet', ampMin: 0, ampMax: 100, smooth: false, grid: true, cursors: true, trueScale: true,
    fmt: 'png', scale: 2, exBg: 'white', exCursors: false, exMeasure: true, cropPlans: true,
    roi: { x0: 0, x1: X, y0: 0, y1: Y }                     // région d'intérêt en coordonnées des données (profondeur : zlo, zhi)
  };
  const st = Object.assign({}, DEF, { pts: [], camera: null, roi: Object.assign({}, DEF.roi) });
  st.camera = Object.assign({ center: { x: 0, y: 0, z: -0.08 } }, CAMS.iso);

  // miroir : x_affiche = X - x_donnees (et inversement)
  const fx = x => st.flipX ? X - x : x;
  const fy = y => st.flipY ? Y - y : y;
  const ixData = i => st.flipX ? NX - 1 - i : i;
  const iyData = j => st.flipY ? NY - 1 - j : j;
  const toDisp = p => ({ x: fx(p.x), y: fy(p.y), z: p.z });
  // région d'intérêt : bornes affichées (après miroir) et plages des axes des coupes
  const rx = () => { const a = fx(st.roi.x0), b = fx(st.roi.x1); return [Math.min(a, b), Math.max(a, b)]; };
  const ry = () => { const a = fy(st.roi.y0), b = fy(st.roi.y1); return [Math.min(a, b), Math.max(a, b)]; };
  const roiFull = () => st.roi.x0 <= 0 && st.roi.x1 >= X && st.roi.y0 <= 0 && st.roi.y1 >= Y;
  const rangeX = () => st.cropPlans && !roiFull() ? rx() : RX;
  const rangeY = () => st.cropPlans && !roiFull() ? ry() : RY;
  function setROI(r) {
    const q = Object.assign({}, st.roi, r);
    const fix = (a, b, L, step) => { let lo = clamp(Math.min(a, b), 0, L), hi = clamp(Math.max(a, b), 0, L);
      if (hi - lo < 2 * step) { const c = (lo + hi) / 2; lo = clamp(c - step, 0, L); hi = clamp(c + step, 0, L); } return [lo, hi]; };
    const [x0, x1] = fix(q.x0, q.x1, X, M.dx), [y0, y1] = fix(q.y0, q.y1, Y, M.dy);
    st.roi = { x0, x1, y0, y1 };
    st.cx = clamp(st.cx, x0, x1); st.cy = clamp(st.cy, y0, y1);
    invalidate('layout', 'three', 'ui');
  }
  function setZ(a, b) {
    st.zlo = clamp(Math.min(a, b), Z0, Z1); st.zhi = clamp(Math.max(a, b), Z0, Z1);
    if (st.zhi - st.zlo < 0.3) st.zhi = Math.min(Z1, st.zlo + 0.3);
  }

  // ---------------------------------------------------------------- calculs mis en cache
  let topCache = {}, sideCache = {}, endCache = {}, meshCache = {};
  function topData() {
    const k0 = kOf(st.zlo), k1 = kOf(st.zhi), key = [k0, k1, st.flipX, st.flipY].join();
    if (topCache.key === key) return topCache;
    const z = [], arg = new Int16Array(NX * NY);
    for (let j = 0; j < NY; j++) {
      const row = new Array(NX), iy = iyData(j);
      for (let i = 0; i < NX; i++) {
        const base = (ixData(i) * NY + iy) * NS; let mx = 0, ka = k0;
        for (let k = k0; k <= k1; k++) { const v = S[base + k]; if (v > mx) { mx = v; ka = k; } }
        row[i] = mx * (100 / 255); arg[j * NX + i] = ka;
      }
      z.push(row);
    }
    topCache = { key, z, arg };
    return topCache;
  }
  const idxX = () => clamp(Math.round(st.cx / M.dx), 0, NX - 1);
  const idxY = () => clamp(Math.round(st.cy / M.dy), 0, NY - 1);
  function sideData() {
    const iy = idxY(), key = [iy, st.flipX].join();
    if (sideCache.key === key) return sideCache;
    const z = [];
    for (let k = 0; k < NS; k++) { const row = new Array(NX); for (let i = 0; i < NX; i++) row[i] = amp(ixData(i), iy, k); z.push(row); }
    sideCache = { key, z };
    return sideCache;
  }
  function endData() {
    const ix = idxX(), key = [ix, st.flipY].join();
    if (endCache.key === key) return endCache;
    const z = [];
    for (let k = 0; k < NS; k++) { const row = new Array(NY); for (let j = 0; j < NY; j++) row[j] = amp(ix, iyData(j), k); z.push(row); }
    endCache = { key, z };
    return endCache;
  }
  function getMesh() {
    const R = st.roi, key = [st.thr, st.zlo, st.zhi, st.flipX, st.flipY, R.x0, R.x1, R.y0, R.y1].join();
    if (meshCache.key === key) return meshCache;
    const m = meshes[st.thr], n = m.nv;
    const x = new Float32Array(n), y = new Float32Array(n);
    for (let i = 0; i < n; i++) { x[i] = st.flipX ? X - m.x[i] : m.x[i]; y[i] = st.flipY ? Y - m.y[i] : m.y[i]; }
    const f = m.f, z = m.z, mxd = m.x, myd = m.y, I = new Uint32Array(m.nf), J = new Uint32Array(m.nf), K = new Uint32Array(m.nf);
    const inBox = v => z[v] >= st.zlo && z[v] <= st.zhi && mxd[v] >= R.x0 && mxd[v] <= R.x1 && myd[v] >= R.y0 && myd[v] <= R.y1;
    const swap = st.flipX !== st.flipY;                      // un miroir inverse l'orientation des triangles
    let c = 0;
    for (let t = 0; t < m.nf; t++) {
      const a = f[3 * t], b = f[3 * t + 1], d = f[3 * t + 2];
      if (inBox(a) && inBox(b) && inBox(d)) {
        I[c] = a; J[c] = swap ? d : b; K[c] = swap ? b : d; c++;
      }
    }
    if (!c) { I[0] = J[0] = K[0] = 0; c = 1; }                  // aucun triangle dans la boîte : triangle dégénéré
    meshCache = { key, x, y, z, I: I.slice(0, c), J: J.slice(0, c), K: K.slice(0, c) };
    return meshCache;
  }

  // ---------------------------------------------------------------- mesure
  const dist3 = (a, b) => Math.hypot(b.x - a.x, b.y - a.y, b.z - a.z);
  function measureText() {
    const [a, b] = st.pts;
    const pt = p => 'scan ' + fmt(p.x, 2) + ', index ' + fmt(p.y, 2) + ', prof. ' + fmt(p.z, 2) + ' mm';
    const r = { a: a ? pt(a) : 'pas encore placé', b: b ? pt(b) : 'pas encore placé', main: '', sub: '', copy: '' };
    if (!a) { r.main = '—'; r.sub = st.mode === 'measure' ? 'Clique le point A.' : 'Passe en mode Mesurer, puis clique deux points.'; }
    else if (!b) { r.main = 'A placé'; r.sub = 'Clique le point B.'; }
    else {
      const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z, d = dist3(a, b), dp = Math.hypot(dx, dy);
      r.main = fmt(d, 2) + ' mm';
      r.sub = 'Δscan ' + fmt(Math.abs(dx), 2) + ' · Δindex ' + fmt(Math.abs(dy), 2) + ' · Δprof. ' + fmt(Math.abs(dz), 2) +
        ' mm. Dans le plan de la surface : ' + fmt(dp, 2) + ' mm.';
      r.copy = 'Distance A-B : ' + fmt(d, 2) + ' mm (Δscan ' + fmt(Math.abs(dx), 2) + ' mm, Δindex ' + fmt(Math.abs(dy), 2) +
        ' mm, Δprofondeur ' + fmt(Math.abs(dz), 2) + ' mm). A : ' + pt(a) + '. B : ' + pt(b) + '.';
    }
    return r;
  }
  function updateReadout() {
    const r = measureText();
    document.querySelectorAll('.rd-a').forEach(e => e.textContent = r.a);
    document.querySelectorAll('.rd-b').forEach(e => e.textContent = r.b);
    document.querySelectorAll('.rd-main').forEach(e => e.textContent = r.main);
    document.querySelectorAll('.rd-sub').forEach(e => e.textContent = r.sub);
  }
  function addPoint(p) {
    if (st.pts.length >= 2) st.pts = [];
    st.pts.push({ x: p.x, y: p.y, z: p.z });
    invalidate('three', 'top', 'side', 'end', 'ui');
  }

  // ---------------------------------------------------------------- figures (données + mise en page)
  function palette3(name) { return PAL[name] || PAL.defaut; }
  function bgSpec(o, is3d) {
    if (o.light) return st.exBg === 'transparent' ? { bg: 'rgba(0,0,0,0)', ink: '#111111', rule: '#cfd6dd', soft: '#555555' }
                                                   : { bg: '#ffffff', ink: '#111111', rule: '#cfd6dd', soft: '#555555' };
    if (is3d) {
      const t = { blanc: { bg: '#ffffff', ink: '#111111', rule: '#cfd6dd', soft: '#555555' },
                  gris: { bg: '#e6ebf0', ink: '#111111', rule: '#c3ccd5', soft: '#555555' },
                  noir: { bg: '#0b0f14', ink: '#e8edf2', rule: '#2b3641', soft: '#9aa8b5' } }[st.bg];
      if (t) return t;
    }
    return { bg: 'rgba(0,0,0,0)', ink: css('--ink'), rule: css('--rule'), soft: css('--ink-soft') };
  }
  const MARG = { l: 58, r: 16, t: 26, b: 44 };
  function axis2d(c, title, range, extra) {
    return Object.assign({ title: { text: title, standoff: 4 }, range, gridcolor: c.rule, zerolinecolor: c.rule, linecolor: c.rule,
      showgrid: st.grid, ticks: 'outside', tickcolor: c.rule, fixedrange: true, showline: true, mirror: true, tickfont: { size: 12 } }, extra || {});
  }
  function layout2d(c, w, h, o, mr) {
    const m = Object.assign({}, MARG); if (mr) m.r += mr;
    return { width: w + (mr || 0), height: h, margin: m, paper_bgcolor: c.bg, plot_bgcolor: 'rgba(0,0,0,0)', separators: ', ',
      font: { family: FONT, color: c.ink, size: 13 }, showlegend: false, hovermode: 'closest', dragmode: false, shapes: [] };
  }
  const line = (x0, y0, x1, y1, color, dash, w) => ({ type: 'line', xref: 'x', yref: 'y', x0, y0, x1, y1, line: { color, width: w || 1.6, dash: dash || 'solid' } });
  function heat(z, x, y, hover, o) {
    return { type: 'heatmap', x, y, z, zmin: st.ampMin, zmax: Math.max(st.ampMax, st.ampMin + 1), colorscale: AMP[st.ampPalette], zsmooth: st.smooth ? 'best' : false,
      showscale: !!(o && o.cbar), colorbar: { title: { text: 'Amplitude (%)', side: 'right' }, thickness: 12, len: 1, outlinewidth: 0, x: 1.02, xpad: 0 },
      hovertemplate: hover + '<extra></extra>' };
  }
  function ptsTrace(proj, fade, o) {
    const P = st.pts.map(toDisp); if (!P.length || o.measure === false) return [];
    const cols = [C_A, C_B].slice(0, P.length), out = [];
    if (P.length === 2) out.push({ type: 'scatter', mode: 'lines', x: P.map(p => proj(p)[0]), y: P.map(p => proj(p)[1]), line: { color: '#ffffff', width: 3 }, hoverinfo: 'skip' },
      { type: 'scatter', mode: 'lines', x: P.map(p => proj(p)[0]), y: P.map(p => proj(p)[1]), line: { color: '#111111', width: 1.4, dash: 'dot' }, hoverinfo: 'skip' });
    out.push({ type: 'scatter', mode: 'markers+text', x: P.map(p => proj(p)[0]), y: P.map(p => proj(p)[1]), text: ['A', 'B'].slice(0, P.length),
      textposition: 'top right', textfont: { size: 14, color: cols }, marker: { size: 10, color: cols, line: { color: '#ffffff', width: 1.5 }, opacity: P.map(p => fade(p)) }, hoverinfo: 'skip' });
    return out;
  }
  function winLines(x0, x1) { return [line(x0, st.zlo, x1, st.zlo, C_WIN, 'dot', 1.3), line(x0, st.zhi, x1, st.zhi, C_WIN, 'dot', 1.3)]; }
  const RX = [-M.dx / 2, X + M.dx / 2], RY = [-M.dy / 2, Y + M.dy / 2], RZ = [Z1 + M.dz / 2, Z0 - M.dz / 2];

  const SHADE = 'rgba(110,120,130,0.55)';
  const shadeRect = (x0, x1, y0, y1) => ({ type: 'rect', xref: 'x', yref: 'y', x0, x1, y0, y1, fillcolor: SHADE, line: { width: 0 }, layer: 'above' });
  function shadeX(yr, full, roiR) {               // bandes hors de la région le long de l'axe horizontal (coupes non recadrées)
    if (st.cropPlans || roiFull()) return [];
    return [shadeRect(full[0], roiR[0], yr[0], yr[1]), shadeRect(roiR[1], full[1], yr[0], yr[1])].filter(r => r.x1 > r.x0);
  }
  function shadeZ(xr) { return [shadeRect(xr[0], xr[1], RZ[1], st.zlo), shadeRect(xr[0], xr[1], st.zhi, RZ[0])]; }
  function figTop(o) {
    o = o || {}; const c = bgSpec(o, false), sz = SZ.top, t = topData();
    const L = layout2d(c, sz.w, sz.h, o, o.cbar ? 70 : 0);
    const D = [heat(t.z, xs, ys, 'scan %{x:.1f} mm<br>index %{y:.1f} mm<br>%{z:.0f} %', o)];
    const gx = rangeX(), gy = rangeY();
    L.xaxis = axis2d(c, 'Scan (mm)', gx); L.yaxis = axis2d(c, 'Index (mm)', gy);
    if (st.cursors && o.cursors !== false) {
      const r1 = rx(), r2 = ry();
      L.shapes.push(...shadeX(gy, RX, r1), ...(st.cropPlans || roiFull() ? [] : [shadeRect(r1[0], r1[1], RY[0], r2[0]), shadeRect(r1[0], r1[1], r2[1], RY[1])]),
        line(fx(st.cx), gy[0], fx(st.cx), gy[1], C_END, 'solid', 2), line(gx[0], fy(st.cy), gx[1], fy(st.cy), C_SIDE, 'solid', 2));
    }
    return { data: D.concat(ptsTrace(p => [p.x, p.y], () => 1, o)), layout: L };
  }
  function figSide(o) {
    o = o || {}; const c = bgSpec(o, false), sz = SZ.side, t = sideData();
    const L = layout2d(c, sz.w, sz.h, o);
    const D = [heat(t.z, xs, depth, 'scan %{x:.1f} mm<br>profondeur %{y:.2f} mm<br>%{z:.0f} %', o)];
    const gx = rangeX();
    L.xaxis = axis2d(c, 'Scan (mm)', gx); L.yaxis = axis2d(c, 'Profondeur (mm)', RZ);
    if (st.cursors && o.cursors !== false) L.shapes.push(...shadeZ(gx), ...shadeX([RZ[1], RZ[0]], RX, rx()),
      line(fx(st.cx), RZ[0], fx(st.cx), RZ[1], C_END, 'solid', 2), ...winLines(gx[0], gx[1]));
    return { data: D.concat(ptsTrace(p => [p.x, p.z], p => Math.abs(p.y - fy(st.cy)) <= M.dy * 0.75 ? 1 : 0.3, o)), layout: L };
  }
  function figEnd(o) {
    o = o || {}; const c = bgSpec(o, false), sz = SZ.end, t = endData();
    const L = layout2d(c, sz.w, sz.h, o);
    const D = [heat(t.z, ys, depth, 'index %{x:.1f} mm<br>profondeur %{y:.2f} mm<br>%{z:.0f} %', o)];
    const gy = rangeY();
    L.xaxis = axis2d(c, 'Index (mm)', gy); L.yaxis = axis2d(c, 'Profondeur (mm)', RZ);
    if (st.cursors && o.cursors !== false) L.shapes.push(...shadeZ(gy), ...shadeX([RZ[1], RZ[0]], RY, ry()),
      line(fy(st.cy), RZ[0], fy(st.cy), RZ[1], C_SIDE, 'solid', 2), ...winLines(gy[0], gy[1]));
    return { data: D.concat(ptsTrace(p => [p.y, p.z], p => Math.abs(p.x - fx(st.cx)) <= M.dx * 0.75 ? 1 : 0.3, o)), layout: L };
  }
  function rect3(fixed, val, c, w) {
    const [a, b] = rx(), [p, q] = ry(), u = st.zlo, v = st.zhi;
    const A = fixed === 'y' ? [[a, val, u], [b, val, u], [b, val, v], [a, val, v], [a, val, u]]
                            : [[val, p, u], [val, q, u], [val, q, v], [val, p, v], [val, p, u]];
    return { type: 'scatter3d', mode: 'lines', x: A.map(p => p[0]), y: A.map(p => p[1]), z: A.map(p => p[2]), line: { color: c, width: w || 4 }, hoverinfo: 'skip', showlegend: false };
  }
  function fig3d(o) {
    o = o || {}; const c = bgSpec(o, true), mc = getMesh();
    const t = { type: 'mesh3d', x: mc.x, y: mc.y, z: mc.z, i: mc.I, j: mc.J, k: mc.K, opacity: st.opacity, flatshading: st.flat,
      lighting: LIGHT[st.light], lightposition: { x: 100, y: -200, z: 300 }, name: 'surface',
      hovertemplate: 'scan %{x:.1f} mm<br>index %{y:.1f} mm<br>profondeur %{z:.2f} mm<extra></extra>' };
    if (st.colorMode === 'depth') {
      Object.assign(t, { intensity: mc.z, cmin: st.zlo, cmax: st.zhi, colorscale: palette3(st.palette), showscale: st.colorbar,
        colorbar: { title: { text: 'Profondeur (mm)', side: 'right' }, thickness: 12, len: 0.7, outlinewidth: 0 } });
    } else { t.color = st.solid; t.showscale = false; }
    const D = [t];
    if (st.box) {
      const [a, b] = rx(), [p, q] = ry(), u = st.zlo, v = st.zhi;
      const P = [[a, p, u], [b, p, u], [b, q, u], [a, q, u], [a, p, v], [b, p, v], [b, q, v], [a, q, v]];
      const E = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]], bx = [], by = [], bz = [];
      E.forEach(e => { [0, 1].forEach(q => { bx.push(P[e[q]][0]); by.push(P[e[q]][1]); bz.push(P[e[q]][2]); }); bx.push(null); by.push(null); bz.push(null); });
      D.push({ type: 'scatter3d', mode: 'lines', x: bx, y: by, z: bz, line: { color: c.soft, width: 2 }, hoverinfo: 'skip', showlegend: false });
    }
    if (st.cuts3d && st.cursors && o.cursors !== false) D.push(rect3('y', fy(st.cy), C_SIDE), rect3('x', fx(st.cx), C_END));
    if (st.pts.length && o.measure !== false) {
      const P = st.pts.map(toDisp), cols = [C_A, C_B].slice(0, P.length);
      D.push({ type: 'scatter3d', mode: 'markers+text', x: P.map(p => p.x), y: P.map(p => p.y), z: P.map(p => p.z), text: ['A', 'B'].slice(0, P.length),
        textposition: 'top center', textfont: { size: 15, color: cols }, marker: { size: 5, color: cols, line: { color: '#ffffff', width: 1 } }, hoverinfo: 'skip', showlegend: false });
      if (P.length === 2) {
        const mid = { x: (P[0].x + P[1].x) / 2, y: (P[0].y + P[1].y) / 2, z: (P[0].z + P[1].z) / 2 };
        D.push({ type: 'scatter3d', mode: 'lines+text', x: [P[0].x, mid.x, P[1].x], y: [P[0].y, mid.y, P[1].y], z: [P[0].z, mid.z, P[1].z],
          text: ['', fmt(dist3(st.pts[0], st.pts[1]), 2) + ' mm', ''], textposition: 'top center', textfont: { size: 14, color: c.ink },
          line: { color: c.ink, width: 5 }, hoverinfo: 'skip', showlegend: false });
      }
    }
    const [ra, rb] = rx(), [rp, rq] = ry(), sx = rb - ra, sy = rq - rp, zs = (st.zhi - st.zlo) * st.zstretch, mx = Math.max(sx, sy, zs);
    const axs = (title, range) => ({ title: { text: title }, range, visible: st.axes, gridcolor: c.rule, zerolinecolor: c.rule, showbackground: false, color: c.ink, showspikes: false, tickfont: { size: 11 }, nticks: 6 });
    const L = { paper_bgcolor: c.bg, font: { family: FONT, color: c.ink, size: 13 }, margin: { l: 0, r: 0, t: 0, b: 0 }, uirevision: 'keep', separators: ', ', showlegend: false,
      scene: { xaxis: axs('Scan (mm)', [ra, rb]), yaxis: axs('Index (mm)', [rp, rq]), zaxis: axs('Profondeur (mm)', [st.zhi, st.zlo]),
        bgcolor: c.bg === 'rgba(0,0,0,0)' ? 'rgba(0,0,0,0)' : c.bg, aspectmode: 'manual', aspectratio: { x: sx / mx, y: sy / mx, z: zs / mx },
        camera: Object.assign({}, st.camera, { projection: { type: st.projection } }) } };
    if (o.w) { L.width = o.w; L.height = o.h; }
    return { data: D, layout: L };
  }
  const FIG = { three: fig3d, top: figTop, side: figSide, end: figEnd };
  const IDS = { three: 'view3d', top: 'pTop', side: 'pSide', end: 'pEnd' };
  const CFG2 = { displaylogo: false, displayModeBar: false, responsive: false, scrollZoom: false, doubleClick: false };
  const CFG3 = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['toImage', 'resetCameraLastSave3d', 'tableRotation', 'orbitRotation', 'hoverClosest3d'] };

  // ---------------------------------------------------------------- curseurs à glissière
  function makeSlider(o) {
    const el = document.createElement('div'); el.className = 'sl ' + (o.vertical ? 'v' : 'h'); el.style.setProperty('--c', o.color);
    const tr = document.createElement('div'); tr.className = 'tr'; el.appendChild(tr);
    const n = o.values.length, vals = o.values.slice(), ths = [];
    const fill = document.createElement('div'); fill.className = 'fill'; if (n === 2) tr.appendChild(fill);
    for (let t = 0; t < n; t++) {
      const th = document.createElement('div'); th.className = 'th'; th.tabIndex = 0; th.setAttribute('role', 'slider');
      th.setAttribute('aria-label', o.labels[t]); th.setAttribute('aria-valuemin', o.min); th.setAttribute('aria-valuemax', o.max);
      el.appendChild(th); ths.push(th);
    }
    const frac = v => (v - o.min) / (o.max - o.min);
    const along = f => (o.vertical && !o.fromTop) ? 1 - f : f;
    function paint() {
      ths.forEach((th, t) => {
        const p = along(clamp(frac(vals[t]), 0, 1)) * 100;
        if (o.vertical) { th.style.top = p + '%'; th.style.left = '13px'; } else { th.style.left = p + '%'; th.style.top = '13px'; }
        th.setAttribute('aria-valuenow', (+vals[t]).toFixed(2));
      });
      if (n === 2) {
        const a = along(clamp(frac(vals[0]), 0, 1)) * 100, b = along(clamp(frac(vals[1]), 0, 1)) * 100, lo = Math.min(a, b), hi = Math.max(a, b);
        if (o.vertical) { fill.style.top = lo + '%'; fill.style.height = (hi - lo) + '%'; fill.style.left = 0; fill.style.width = '100%'; }
        else { fill.style.left = lo + '%'; fill.style.width = (hi - lo) + '%'; fill.style.top = 0; fill.style.height = '100%'; }
      }
    }
    function snap(v) { v = clamp(v, o.min, o.max); return o.step ? clamp(Math.round((v - o.min) / o.step) * o.step + o.min, o.min, o.max) : v; }
    function valueAt(ev) {
      const r = tr.getBoundingClientRect(); let f = o.vertical ? (ev.clientY - r.top) / r.height : (ev.clientX - r.left) / r.width;
      f = clamp(f, 0, 1); if (o.vertical && !o.fromTop) f = 1 - f; return snap(o.min + f * (o.max - o.min));
    }
    let active = -1;
    function move(v) {
      if (n === 2) { const g = o.minGap || 0; v = active === 0 ? Math.min(v, vals[1] - g) : Math.max(v, vals[0] + g); }
      if (v === vals[active]) return; vals[active] = v; paint(); o.onInput(vals.slice());
    }
    el.addEventListener('pointerdown', ev => {
      ev.preventDefault(); const v = valueAt(ev);
      active = n === 1 ? 0 : (Math.abs(v - vals[0]) <= Math.abs(v - vals[1]) ? 0 : 1);
      try { el.setPointerCapture(ev.pointerId); } catch (e) { /* ignore */ } move(v);
    });
    el.addEventListener('pointermove', ev => { if (active >= 0) move(valueAt(ev)); });
    const end = () => { active = -1; }; el.addEventListener('pointerup', end); el.addEventListener('pointercancel', end);
    ths.forEach((th, t) => th.addEventListener('keydown', ev => {
      const step = (o.step || (o.max - o.min) / 100) * (ev.shiftKey ? 10 : 1); let d = 0;
      if (ev.key === 'ArrowLeft' || ev.key === 'ArrowDown') d = -step; else if (ev.key === 'ArrowRight' || ev.key === 'ArrowUp') d = step; else return;
      ev.preventDefault(); if (o.vertical && o.fromTop) d = -d; active = t; move(snap(vals[t] + d)); active = -1;
    }));
    paint();
    return { el, set(v) { v.forEach((x, i) => { vals[i] = x; }); paint(); } };
  }

  // ---------------------------------------------------------------- disposition des plans
  let SZ = null, SL = {}, layoutKey = '';
  function computeSizes(avail) {
    const m = MARG, zr = Z1 - Z0, narrow = avail < 760;
    if (narrow) {
      const w = Math.max(280, avail - 8), mk = h => ({ w, h, iw: w - m.l - m.r, ih: h - m.t - m.b });
      return { narrow: true, top: mk(300), side: mk(260), end: mk(260), colB: 0, rowS: 0 };
    }
    const colB = 38, rowS = 40, zs = zr * st.zstretch, gx = rangeX(), gy = rangeY(), SX = gx[1] - gx[0], SY = gy[1] - gy[0];
    let s = (avail - 12 - colB - 2 * (m.l + m.r)) / (SX + SY);
    s = Math.min(s, (800 - 2 * (m.t + m.b) - rowS) / (SY + zs));
    let iw1, iw2, ih1, ih2;
    if (st.trueScale) { iw1 = SX * s; iw2 = Math.max(SY * s, 190); ih1 = Math.max(SY * s, 150); ih2 = Math.max(zs * s, 150); }
    else { iw1 = (avail - 12 - colB - 2 * (m.l + m.r)) * 0.72; iw2 = (avail - 12 - colB - 2 * (m.l + m.r)) * 0.28; ih1 = 230; ih2 = 230; }
    const mk = (iw, ih) => ({ w: Math.round(iw + m.l + m.r), h: Math.round(ih + m.t + m.b), iw: Math.round(iw), ih: Math.round(ih) });
    const sz = { narrow: false, top: mk(iw1, ih1), side: mk(iw1, ih2), end: mk(iw2, ih2), colB, rowS };
    sz.scale = s; sz.depthK = st.trueScale ? ih2 / (zr * s) : null;
    return sz;
  }
  function el(tag, cls, parent) { const e = document.createElement(tag); if (cls) e.className = cls; if (parent) parent.appendChild(e); return e; }
  function plotCell(parent, id, sz, cap, sub, clsx) {
    const c = el('div', 'cell ' + (clsx || ''), parent); const d = el('div', '', c); d.id = id; d.style.width = sz.w + 'px'; d.style.height = sz.h + 'px';
    const k = el('div', 'cap', c); k.textContent = cap; if (sub) { const i = el('i', '', k); i.textContent = sub; }
    return c;
  }
  function label(parent, css, name, valId) {
    const l = el('div', 'slbl', parent); Object.assign(l.style, css); const b = el('b', '', l); b.textContent = name; const v = el('span', '', l); v.id = valId; return l;
  }
  function readoutCard(parent) {
    const c = el('div', 'corner', parent);
    c.innerHTML = '<div class="sub">Distance A-B</div><div class="big rd-main"></div><div class="sub rd-sub"></div>' +
      '<div class="legend"><div>Amplitude</div><div class="bar" id="legBar"></div><div class="ends"><span id="legMin"></span><span id="legMax"></span></div>' +
      '<div id="legScale" style="margin-top:6px"></div><div id="roiTxt" style="margin-top:4px"></div></div>';
    return c;
  }
  function buildPlans() {
    const host = $('plans'); ['pTop', 'pSide', 'pEnd'].forEach(id => { const d = $(id); if (d) { try { Plotly.purge(d); } catch (e) { /* ignore */ } } });
    host.innerHTML = ''; SL = {};
    const avail = Math.floor(host.clientWidth || host.parentElement.clientWidth || 900);
    SZ = computeSizes(avail); const m = MARG;
    const g = el('div', 'pgrid', host);
    const gx = rangeX(), gy = rangeY();                                  // glissières sur la même étendue que les axes
    const cxSl = { min: gx[0], max: gx[1], color: C_END, labels: ['Position en scan'] };
    const cySl = { min: gy[0], max: gy[1], color: C_SIDE, labels: ['Position en index'] };
    const zSl = { min: RZ[1], max: RZ[0], step: 0.05, color: C_WIN, labels: ['Profondeur min', 'Profondeur max'], minGap: 0.3, fromTop: true };
    const onCx = v => { st.cx = clamp(fx(v[0]), st.roi.x0, st.roi.x1); invalidate('end', 'top', 'side', 'three', 'ui'); };
    const onCy = v => { st.cy = clamp(fy(v[0]), st.roi.y0, st.roi.y1); invalidate('side', 'top', 'end', 'three', 'ui'); };
    const onZ = v => { setZ(v[0], v[1]); invalidate('top', 'side', 'end', 'three', 'ui'); };
    if (SZ.narrow) {
      g.style.gridTemplateColumns = SZ.top.w + 'px';
      plotCell(g, 'pTop', SZ.top, 'Vue de dessus', 'scan × index', 'top');
      const sx = el('div', 'cell', g); sx.style.height = '46px'; label(sx, { left: '0', top: '0', width: '100%' }, 'Position en scan (coupe transversale)', 'lblCx');
      SL.cx = makeSlider(Object.assign({ vertical: false, values: [fx(st.cx)], onInput: onCx }, cxSl)); sx.appendChild(SL.cx.el);
      Object.assign(SL.cx.el.style, { left: m.l + 'px', width: SZ.top.iw + 'px', top: '18px' });
      const sy = el('div', 'cell', g); sy.style.height = '46px'; label(sy, { left: '0', top: '0', width: '100%' }, 'Position en index (coupe longitudinale)', 'lblCy');
      SL.cy = makeSlider(Object.assign({ vertical: false, values: [fy(st.cy)], onInput: onCy }, cySl)); sy.appendChild(SL.cy.el);
      Object.assign(SL.cy.el.style, { left: m.l + 'px', width: SZ.top.iw + 'px', top: '18px' });
      plotCell(g, 'pSide', SZ.side, 'Coupe longitudinale', 'scan × profondeur', 'side');
      plotCell(g, 'pEnd', SZ.end, 'Coupe transversale', 'index × profondeur', 'end');
      const sz = el('div', 'cell', g); sz.style.height = '46px'; label(sz, { left: '0', top: '0', width: '100%' }, 'Fenêtre de profondeur', 'lblZ');
      SL.z = makeSlider(Object.assign({ vertical: false, values: [st.zlo, st.zhi], onInput: onZ }, zSl, { fromTop: false })); sz.appendChild(SL.z.el);
      Object.assign(SL.z.el.style, { left: m.l + 'px', width: SZ.top.iw + 'px', top: '18px' });
      readoutCard(g);
    } else {
      g.style.gridTemplateColumns = SZ.top.w + 'px ' + SZ.colB + 'px ' + SZ.end.w + 'px';
      g.style.gridTemplateRows = SZ.top.h + 'px ' + SZ.rowS + 'px ' + SZ.side.h + 'px';
      plotCell(g, 'pTop', SZ.top, 'Vue de dessus', 'scan × index', 'top');
      const cyHost = el('div', 'cell', g);                                 // glissière d'index, le long de l'axe vertical de la vue de dessus
      label(cyHost, { left: '-6px', top: '2px', width: SZ.colB + 12 + 'px' }, 'Index', 'lblCy');
      SL.cy = makeSlider(Object.assign({ vertical: true, values: [fy(st.cy)], onInput: onCy }, cySl)); cyHost.appendChild(SL.cy.el);
      Object.assign(SL.cy.el.style, { top: m.t + 'px', height: SZ.top.ih + 'px', left: (SZ.colB - 26) / 2 + 'px' });
      readoutCard(el('div', 'cell', g));                                   // case en haut à droite : mesure et légende
      const cxHost = el('div', 'cell', g);                                 // glissière de scan, entre la vue de dessus et la coupe longitudinale
      label(cxHost, { left: '0', top: '2px', width: m.l - 6 + 'px' }, 'Scan', 'lblCx');
      SL.cx = makeSlider(Object.assign({ vertical: false, values: [fx(st.cx)], onInput: onCx }, cxSl)); cxHost.appendChild(SL.cx.el);
      Object.assign(SL.cx.el.style, { left: m.l + 'px', width: SZ.top.iw + 'px', top: '7px' });
      el('div', 'cell', g); el('div', 'cell', g);
      plotCell(g, 'pSide', SZ.side, 'Coupe longitudinale', 'scan × profondeur', 'side');
      const zHost = el('div', 'cell', g);                                  // fenêtre de profondeur, entre les deux coupes qui partagent l'axe de profondeur
      label(zHost, { left: '-6px', top: '2px', width: SZ.colB + 12 + 'px' }, 'Prof.', 'lblZ');
      SL.z = makeSlider(Object.assign({ vertical: true, values: [st.zlo, st.zhi], onInput: onZ }, zSl)); zHost.appendChild(SL.z.el);
      Object.assign(SL.z.el.style, { top: m.t + 'px', height: SZ.side.ih + 'px', left: (SZ.colB - 26) / 2 + 'px' });
      plotCell(g, 'pEnd', SZ.end, 'Coupe transversale', 'index × profondeur', 'end');
    }
    layoutKey = key();
    ['top', 'side', 'end'].forEach(k => { $(IDS[k]).__hooked = false; cropDrag($(IDS[k]).parentNode, k); });
  }
  const key = () => [Math.floor($('plans').clientWidth), st.trueScale, st.zstretch, st.cropPlans, st.roi.x0, st.roi.x1, st.roi.y0, st.roi.y1].join();

  // rognage à la souris : rectangle tracé dans une coupe (mode Rogner)
  function cropDrag(cell, k) {
    let p0 = null, rb = null;
    const rel = ev => {
      const d = $(IDS[k]), fl = d && d._fullLayout; if (!fl || !fl.xaxis) return null;
      const r = d.getBoundingClientRect(), xa = fl.xaxis, ya = fl.yaxis;
      const px = clamp(ev.clientX - r.left - xa._offset, 0, xa._length), py = clamp(ev.clientY - r.top - ya._offset, 0, ya._length);
      return { px: px + xa._offset, py: py + ya._offset, x: xa.p2l(px), y: ya.p2l(py) };
    };
    const paint = q => { if (!q || !rb) return; Object.assign(rb.style, { left: Math.min(p0.px, q.px) + 'px', top: Math.min(p0.py, q.py) + 'px',
      width: Math.abs(q.px - p0.px) + 'px', height: Math.abs(q.py - p0.py) + 'px' }); };
    const stop = () => { p0 = null; if (rb) { rb.remove(); rb = null; } };
    cell.addEventListener('pointerdown', ev => {
      if (st.mode !== 'crop' || ev.button !== 0) return; const q = rel(ev); if (!q) return;
      ev.preventDefault(); ev.stopPropagation(); p0 = q; rb = el('div', 'rb', cell);
      try { cell.setPointerCapture(ev.pointerId); } catch (e) { /* ignore */ } paint(q);
    }, true);
    cell.addEventListener('pointermove', ev => { if (p0) paint(rel(ev)); }, true);
    cell.addEventListener('pointerup', ev => {
      if (!p0) return; const a = p0, q = rel(ev); stop();
      if (!q || Math.abs(q.px - a.px) < 6 || Math.abs(q.py - a.py) < 6) return;   // trop petit : pas de rognage
      if (k === 'top') setROI({ x0: fx(a.x), x1: fx(q.x), y0: fy(a.y), y1: fy(q.y) });
      else if (k === 'side') { setZ(a.y, q.y); setROI({ x0: fx(a.x), x1: fx(q.x) }); }
      else { setZ(a.y, q.y); setROI({ y0: fy(a.x), y1: fy(q.x) }); }
    }, true);
    cell.addEventListener('pointercancel', stop, true);
  }
  function fitROI() {                                                     // boîte englobante de la surface visible
    const m = meshes[st.thr], f = m.f; let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity, z0 = Infinity, z1 = -Infinity;
    const ok = v => m.z[v] >= st.zlo && m.z[v] <= st.zhi;
    for (let t = 0; t < m.nf; t++) {
      const a = f[3 * t], b = f[3 * t + 1], d = f[3 * t + 2]; if (!(ok(a) && ok(b) && ok(d))) continue;
      for (const v of [a, b, d]) { x0 = Math.min(x0, m.x[v]); x1 = Math.max(x1, m.x[v]); y0 = Math.min(y0, m.y[v]); y1 = Math.max(y1, m.y[v]); z0 = Math.min(z0, m.z[v]); z1 = Math.max(z1, m.z[v]); }
    }
    if (!isFinite(x0)) { toast('Aucune surface dans la fenêtre de profondeur.'); return; }
    setZ(Math.max(st.zlo, z0 - 0.3), Math.min(st.zhi, z1 + 0.3));
    setROI({ x0: x0 - 1, x1: x1 + 1, y0: y0 - 1, y1: y1 + 1 });
  }

  // ---------------------------------------------------------------- dessin
  const dirty = new Set(); let raf = null;
  function invalidate() { for (const p of arguments) dirty.add(p); if (!raf) raf = requestAnimationFrame(flush); }
  function flush() {
    raf = null; const d = new Set(dirty); dirty.clear();
    if (!window.Plotly) return;
    if (d.has('layout') || !SZ) { buildPlans(); ['top', 'side', 'end', 'ui'].forEach(x => d.add(x)); }
    if (d.has('ui')) syncUI();
    ['top', 'side', 'end'].forEach(k => { if (d.has(k)) draw(k); });
    if (d.has('three')) draw('three');
  }
  function hook(k, div) {
    if (div.__hooked) return; div.__hooked = true;
    div.on('plotly_click', ev => { if (ev.points && ev.points.length) onClick(k, ev.points[0]); });
    if (k === 'three') div.on('plotly_relayout', ev => { if (ev && ev['scene.camera']) st.camera = ev['scene.camera']; });
  }
  function draw(k) {
    const div = $(IDS[k]); if (!div) return;
    const f = FIG[k]();
    Plotly.react(div, f.data, f.layout, k === 'three' ? CFG3 : CFG2).then(() => hook(k, div));
  }
  function onClick(k, p) {
    if (st.mode === 'crop') return;
    let x, y, z;
    if (k === 'top') {
      const i = clamp(Math.round(p.x / M.dx), 0, NX - 1), j = clamp(Math.round(p.y / M.dy), 0, NY - 1);
      x = fx(p.x); y = fy(p.y); z = depth[topData().arg[j * NX + i]];
    } else if (k === 'side') { x = fx(p.x); y = st.cy; z = p.y; }
    else if (k === 'end') { x = st.cx; y = fy(p.x); z = p.y; }
    else { x = fx(p.x); y = fy(p.y); z = p.z; }
    if (st.mode === 'measure') { addPoint({ x, y, z }); return; }
    const R = st.roi;
    if (k === 'top' || k === 'three') { st.cx = clamp(x, R.x0, R.x1); st.cy = clamp(y, R.y0, R.y1); }
    else if (k === 'side') st.cx = clamp(x, R.x0, R.x1); else st.cy = clamp(y, R.y0, R.y1);
    invalidate('top', 'side', 'end', 'three', 'ui');
  }

  // ---------------------------------------------------------------- interface
  const setOut = (id, txt) => { const e = $(id); if (e) e.textContent = txt; };
  function syncUI() {
    if (SL.cx) SL.cx.set([fx(st.cx)]); if (SL.cy) SL.cy.set([fy(st.cy)]); if (SL.z) SL.z.set([st.zlo, st.zhi]);
    setOut('lblCx', fmt(fx(st.cx)) + ' mm'); setOut('lblCy', fmt(fy(st.cy)) + ' mm'); setOut('lblZ', fmt(st.zlo) + ' à ' + fmt(st.zhi));
    const set = (id, v) => { const e = $(id); if (e && document.activeElement !== e) e.value = fmt(v, 2); };
    set('nScan', fx(st.cx)); set('nIndex', fy(st.cy)); set('nZlo', st.zlo); set('nZhi', st.zhi);
    set('rX0', rx()[0]); set('rX1', rx()[1]); set('rY0', ry()[0]); set('rY1', ry()[1]);
    $('plans').classList.toggle('crop', st.mode === 'crop');
    setOut('roiTxt', 'Région : scan ' + fmt(rx()[0]) + ' à ' + fmt(rx()[1]) + ', index ' + fmt(ry()[0]) + ' à ' + fmt(ry()[1]) +
      ', prof. ' + fmt(st.zlo) + ' à ' + fmt(st.zhi) + ' mm');
    setOut('thr-out', st.thr + ' %'); setOut('opacity-out', fmt(st.opacity, 2)); setOut('zstretch-out', '×' + fmt(st.zstretch));
    setOut('ampMin-out', st.ampMin + ' %'); setOut('ampMax-out', st.ampMax + ' %');
    document.querySelectorAll('[data-mode]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === st.mode)));
    $('rowPal').style.display = st.colorMode === 'depth' ? '' : 'none'; $('rowSolid').style.display = st.colorMode === 'solid' ? '' : 'none';
    const bar = $('legBar');
    if (bar) { bar.style.background = 'linear-gradient(to right,' + AMP[st.ampPalette].map(s => s[1] + ' ' + (s[0] * 100) + '%').join(',') + ')';
      setOut('legMin', st.ampMin + ' %'); setOut('legMax', st.ampMax + ' %');
      setOut('legScale', SZ && !SZ.narrow && SZ.scale ? 'Échelle : ' + fmt(SZ.scale, 1) + ' px/mm' + (SZ.depthK && Math.abs(SZ.depthK - 1) > 0.02 ? ', profondeur ×' + fmt(SZ.depthK) : '') : ''); }
    updateReadout();
  }
  function toast(t) { const e = $('toast'); e.textContent = t; e.classList.add('on'); clearTimeout(toast.t); toast.t = setTimeout(() => e.classList.remove('on'), 2200); }
  function bindCheck(id, key, then) { const e = $(id); e.checked = !!st[key]; e.addEventListener('change', () => { st[key] = e.checked; (then || (() => invalidate('three')))(); }); }
  function bindSel(id, key, then, num) { const e = $(id); e.value = String(st[key]); e.addEventListener('change', () => { st[key] = num ? +e.value : e.value; (then || (() => invalidate('three')))(); }); }
  function bindRange(id, key, then, num) { const e = $(id); e.value = st[key]; e.addEventListener('input', () => { st[key] = +e.value; (then || (() => invalidate('three')))(); }); }
  function fillSelect(id, names) { const s = $(id); Object.keys(names).forEach(k => { const o = document.createElement('option'); o.value = k; o.textContent = names[k]; s.appendChild(o); }); }

  function initUI() {
    fillSelect('palette', PAL_NAMES); fillSelect('ampPalette', AMP_NAMES);
    const thr = $('thr'); thr.max = levelKeys.length - 1; thr.value = levelKeys.indexOf(st.thr);
    thr.addEventListener('input', () => { st.thr = levelKeys[+thr.value]; invalidate('three', 'ui'); });
    bindCheck('flipX', 'flipX', () => invalidate('top', 'side', 'end', 'three', 'ui'));
    bindCheck('flipY', 'flipY', () => invalidate('top', 'side', 'end', 'three', 'ui'));
    bindSel('proj', 'projection', () => Plotly.relayout($('view3d'), { 'scene.camera.projection.type': st.projection })); bindSel('colorMode', 'colorMode', () => invalidate('three', 'ui')); bindSel('palette', 'palette');
    $('solid').value = st.solid; $('solid').addEventListener('input', () => { st.solid = $('solid').value; invalidate('three'); });
    bindRange('opacity', 'opacity', () => invalidate('three', 'ui'));
    bindSel('light', 'light'); bindCheck('flat', 'flat'); bindCheck('colorbar', 'colorbar');
    bindSel('bg', 'bg'); bindCheck('axes', 'axes'); bindCheck('box', 'box'); bindCheck('cuts3d', 'cuts3d');
    bindRange('zstretch', 'zstretch', () => invalidate('layout', 'three', 'ui'));
    const allPlans = () => invalidate('top', 'side', 'end', 'ui');
    bindSel('ampPalette', 'ampPalette', allPlans);
    bindRange('ampMin', 'ampMin', () => { st.ampMax = Math.max(st.ampMax, st.ampMin + 1); $('ampMax').value = st.ampMax; allPlans(); });
    bindRange('ampMax', 'ampMax', () => { st.ampMin = Math.min(st.ampMin, st.ampMax - 1); $('ampMin').value = st.ampMin; allPlans(); });
    bindCheck('smooth', 'smooth', allPlans); bindCheck('grid', 'grid', allPlans);
    bindCheck('cursors', 'cursors', () => invalidate('top', 'side', 'end', 'three'));
    bindCheck('trueScale', 'trueScale', () => invalidate('layout', 'ui'));
    const num = (id, fn) => $(id).addEventListener('change', () => { const v = parseFloat(String($(id).value).replace(',', '.')); if (isFinite(v)) { fn(v); invalidate('top', 'side', 'end', 'three', 'ui'); } });
    num('nScan', v => { st.cx = fx(clamp(v, 0, X)); }); num('nIndex', v => { st.cy = fy(clamp(v, 0, Y)); });
    num('nZlo', v => { st.zlo = clamp(Math.min(v, st.zhi - 0.3), Z0, Z1); }); num('nZhi', v => { st.zhi = clamp(Math.max(v, st.zlo + 0.3), Z0, Z1); });
    const roiNum = (id, ax, i) => $(id).addEventListener('change', () => {
      const v = parseFloat(String($(id).value).replace(',', '.')); if (!isFinite(v)) { invalidate('ui'); return; }
      const r = ax === 'x' ? rx() : ry(); r[i] = v;
      setROI(ax === 'x' ? { x0: fx(r[0]), x1: fx(r[1]) } : { y0: fy(r[0]), y1: fy(r[1]) });
    });
    roiNum('rX0', 'x', 0); roiNum('rX1', 'x', 1); roiNum('rY0', 'y', 0); roiNum('rY1', 'y', 1);
    bindCheck('cropPlans', 'cropPlans', () => invalidate('layout', 'ui'));
    $('roiFit').addEventListener('click', fitROI);
    $('roiAll').addEventListener('click', () => { setZ(DEF.zlo, DEF.zhi); setROI({ x0: 0, x1: X, y0: 0, y1: Y }); });
    document.querySelectorAll('[data-mode]').forEach(b => b.addEventListener('click', () => { st.mode = b.dataset.mode; invalidate('ui'); }));
    document.querySelectorAll('[data-cam]').forEach(b => b.addEventListener('click', () => {
      st.camera = Object.assign({ center: { x: 0, y: 0, z: -0.08 } }, JSON.parse(JSON.stringify(CAMS[b.dataset.cam])));
      Plotly.relayout($('view3d'), { 'scene.camera': Object.assign({}, st.camera, { projection: { type: st.projection } }) });
    }));
    $('mClear').addEventListener('click', () => { st.pts = []; invalidate('three', 'top', 'side', 'end', 'ui'); });
    $('mCopy').addEventListener('click', () => {
      const t = measureText().copy; if (!t) { toast('Place deux points d\'abord.'); return; }
      try { navigator.clipboard.writeText(t).then(() => toast('Résultat copié.'), () => toast('Copie impossible.')); } catch (e) { toast('Copie impossible.'); }
    });
    bindSel('fmt', 'fmt'); bindSel('scale', 'scale', null, true); bindSel('exBg', 'exBg'); bindCheck('exCursors', 'exCursors', () => {}); bindCheck('exMeasure', 'exMeasure', () => {});
    document.querySelectorAll('[data-save]').forEach(b => b.addEventListener('click', () => saveFigs(b.dataset.save)));
    $('reset').addEventListener('click', () => {
      const pts = st.pts; Object.assign(st, DEF, { pts: pts, roi: Object.assign({}, DEF.roi) }); st.camera = Object.assign({ center: { x: 0, y: 0, z: -0.08 } }, JSON.parse(JSON.stringify(CAMS.iso)));
      ['flipX', 'flipY', 'flat', 'colorbar', 'axes', 'box', 'cuts3d', 'smooth', 'grid', 'cursors', 'trueScale', 'exCursors', 'exMeasure', 'cropPlans'].forEach(k => { $(k).checked = !!st[k]; });
      ['proj:projection', 'colorMode', 'palette', 'light', 'bg', 'ampPalette', 'fmt', 'exBg'].forEach(s => { const [id, k] = s.split(':'); $(id).value = st[k || id]; });
      ['opacity', 'zstretch', 'ampMin', 'ampMax'].forEach(k => { $(k).value = st[k]; }); $('scale').value = String(st.scale); $('solid').value = st.solid; $('thr').value = levelKeys.indexOf(st.thr);
      meshCache = {}; topCache = {}; sideCache = {}; endCache = {};
      Plotly.relayout($('view3d'), { 'scene.camera': st.camera }); invalidate('layout', 'three', 'ui');
    });
  }

  // ---------------------------------------------------------------- enregistrement des figures
  function download(url, name) { const a = document.createElement('a'); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove(); }
  function figOpts() { return { light: true, cursors: st.exCursors, measure: st.exMeasure }; }
  function exportSize(k) {
    const d = $(IDS[k]); const w = Math.round(d.clientWidth || (SZ[k] && SZ[k].w) || 800), h = Math.round(d.clientHeight || (SZ[k] && SZ[k].h) || 600); return { w, h };
  }
  async function renderFig(k, format) {
    const o = figOpts(), sz = exportSize(k);
    if (k === 'top') o.cbar = true;
    const f = FIG[k](o);
    let w = sz.w, h = sz.h;
    if (k === 'top') { w += 70; }
    if (k === 'three') { f.layout.width = w; f.layout.height = h; }
    if (format === 'jpeg') f.layout.paper_bgcolor = '#ffffff';
    return Plotly.toImage(f, { format, width: w, height: h, scale: st.scale });
  }
  const NAMES = { three: 'vue3D', top: 'dessus', side: 'coupe_longitudinale', end: 'coupe_transversale', sheet: 'planche' };
  async function saveOne(k) {
    let format = st.fmt; if (format === 'svg' && k === 'three') { format = 'png'; toast('La vue 3D n\'existe pas en SVG : PNG utilisé.'); }
    const url = await renderFig(k, format);
    download(url, M.name + '_' + NAMES[k] + '.' + (format === 'jpeg' ? 'jpg' : format));
  }
  function loadImg(url) { return new Promise((res, rej) => { const im = new Image(); im.onload = () => res(im); im.onerror = rej; im.src = url; }); }
  async function sheetURL(format) {
    const g = document.querySelector('.pgrid'); if (!g) return null; const gr = g.getBoundingClientRect(), s = st.scale;
    const parts = ['top', 'side', 'end'], rects = {}, imgs = {};
    for (const k of parts) { const r = $(IDS[k]).getBoundingClientRect(); rects[k] = { x: r.left - gr.left, y: r.top - gr.top, w: r.width, h: r.height }; imgs[k] = await loadImg(await renderFig(k, 'png')); }
    const W = Math.round(Math.max(...parts.map(k => rects[k].x + rects[k].w)) + 70), H = Math.round(Math.max(...parts.map(k => rects[k].y + rects[k].h)));
    const cv = document.createElement('canvas'); cv.width = W * s; cv.height = H * s; const ctx = cv.getContext('2d');
    if (st.exBg !== 'transparent' || format === 'jpeg') { ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, cv.width, cv.height); }
    parts.forEach(k => { const r = rects[k]; ctx.drawImage(imgs[k], r.x * s, r.y * s); });
    const mime = format === 'jpeg' ? 'image/jpeg' : format === 'webp' ? 'image/webp' : 'image/png';
    return cv.toDataURL(mime, 0.95);
  }
  async function saveSheet() {
    let format = st.fmt; if (format === 'svg') { format = 'png'; toast('La planche est enregistrée en PNG.'); }
    const url = await sheetURL(format); if (url) download(url, M.name + '_' + NAMES.sheet + '.' + (format === 'jpeg' ? 'jpg' : format));
  }
  async function saveFigs(what) {
    if (!window.Plotly) return;
    try {
      toast('Préparation de l\'image…');
      if (what === 'sheet') await saveSheet();
      else if (what === 'all') { for (const k of ['three', 'top', 'side', 'end']) await saveOne(k); await saveSheet(); }
      else await saveOne(what);
      toast('Image enregistrée.');
    } catch (e) { toast('Enregistrement impossible : ' + (e && e.message ? e.message : e)); }
  }

  // ---------------------------------------------------------------- démarrage
  if (!window.Plotly) { document.querySelector('#view3d .loading').textContent = 'La bibliothèque de graphiques n’a pas pu se charger. Vérifie ta connexion Internet et recharge la page.'; return; }
  document.getElementById('view3d').innerHTML = '';
  initUI(); invalidate('layout', 'three', 'ui');
  let rt = null; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => { if (key() !== layoutKey) invalidate('layout', 'ui'); }, 150); });
  const mq = window.matchMedia('(prefers-color-scheme: dark)'); if (mq.addEventListener) mq.addEventListener('change', () => invalidate('top', 'side', 'end', 'three'));
  new MutationObserver(() => invalidate('top', 'side', 'end', 'three')).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  window.__uv = { st, addPoint, invalidate, flush, FIG, SZ: () => SZ, measureText, onClick, saveFigs, topData, renderFig, sheetURL, setROI, fitROI, rx, ry };
})();
</script>
</body>
</html>
"""


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Reconstruction 3D d\'un .UVData (horloge ou encodeur, détecté seul)')
    ap.add_argument('fichier', nargs='?', help='fichier .UVData (sans argument : ouvre les fenêtres)')
    ap.add_argument('-o', help='fichier de sortie (.html, .vtk ou .npz)')
    ap.add_argument('--inverser', action='store_true', help='inverser le sens du scan')
    ap.add_argument('--sans-interpolation', action='store_true', help='encodeur : laisser les trous vides')
    ap.add_argument('--vitesse', type=float, help='horloge : vitesse de déplacement (mm/s)')
    ap.add_argument('--frequence', type=float, help="horloge : fréquence d'acquisition (Hz)")
    ap.add_argument('--epaisseur', type=float, help='horloge : épaisseur connue (mm) pour recaler la vitesse')
    args = ap.parse_args()
    run_cli(args) if args.fichier else run_gui()
