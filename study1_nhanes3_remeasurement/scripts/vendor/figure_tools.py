from __future__ import annotations
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json
import platform
import tempfile
import warnings
import xml.etree.ElementTree as ET
import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np
COLORS = {'before': '#0072B2', 'after': '#D55E00', 'neutral': '#657381', 'accent': '#009E73'}

@contextmanager
def publication_style(font_size=8):
    if not np.isfinite(font_size) or font_size <= 0:
        raise ValueError('font_size must be positive and finite')
    with mpl.rc_context({'font.family': 'DejaVu Sans', 'font.size': font_size, 'axes.labelsize': font_size, 'axes.titlesize': font_size + 0.5, 'xtick.labelsize': font_size, 'ytick.labelsize': font_size, 'legend.fontsize': font_size, 'axes.linewidth': 0.65, 'lines.linewidth': 1.1, 'lines.markersize': 3.2, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'xtick.direction': 'out', 'ytick.direction': 'out', 'text.usetex': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.bbox': None, 'savefig.facecolor': 'white'}):
        yield

def paired_values(records, *, subject='subject', condition='condition', value='value', conditions=('before', 'after')):
    if len(conditions) != 2 or conditions[0] == conditions[1]:
        raise ValueError('conditions must contain two distinct names')
    grouped = {}
    for row_number, row in enumerate(records, 1):
        try:
            raw_id = row[subject]
            if raw_id is None:
                raise ValueError('empty subject ID')
            key = str(raw_id).strip()
            cond = row[condition]
            number = float(row[value])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'invalid record {row_number}: {exc}') from exc
        if not key or cond not in conditions or (not np.isfinite(number)):
            raise ValueError(f'invalid subject, condition or nonfinite value at row {row_number}')
        pair = grouped.setdefault(key, {})
        if cond in pair:
            raise ValueError(f'duplicate subject-condition: {key}, {cond}; aggregate only with scientific justification')
        pair[cond] = number
    if not grouped:
        raise ValueError('no paired records')
    incomplete = [key for key, pair in grouped.items() if len(pair) != 2]
    if incomplete:
        raise ValueError(f'incomplete pairs; no automatic exclusion: {incomplete[:8]}')
    ids = sorted(grouped)
    return (ids, np.array([grouped[key][conditions[0]] for key in ids]), np.array([grouped[key][conditions[1]] for key in ids]))

def paired_mean_interval(before, after, *, subject_ids, confidence=0.95, n_boot=5000, seed=17):
    before, after = (np.asarray(before, dtype=float), np.asarray(after, dtype=float))
    if before.ndim != 1 or after.shape != before.shape or len(before) < 2:
        raise ValueError('need matching one-dimensional paired vectors with at least two subjects')
    if not np.isfinite(before).all() or not np.isfinite(after).all():
        raise ValueError('paired values must be finite; resolve missingness explicitly')
    ids = [str(x).strip() if x is not None else '' for x in subject_ids]
    if len(ids) != len(before) or len(set(ids)) != len(ids) or (not all(ids)):
        raise ValueError('subject IDs must be distinct, nonempty and match vector length')
    if not 0 < confidence < 1 or not isinstance(n_boot, int) or n_boot < 100:
        raise ValueError('confidence must be between 0 and 1, n_boot an integer >= 100')
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError('provide a nonnegative integer seed')
    with np.errstate(over='ignore', invalid='ignore'):
        delta = after - before
    if not np.isfinite(delta).all():
        raise ValueError('paired differences overflow; rescale units explicitly')
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot)
    batch_size = max(1, min(256, 1000000 // len(delta)))
    for start in range(0, n_boot, batch_size):
        count = min(batch_size, n_boot - start)
        means[start:start + count] = delta[rng.integers(0, len(delta), (count, len(delta)))].mean(axis=1)
    alpha = (1 - confidence) / 2
    if not np.isfinite(means).all():
        raise ValueError('bootstrap means overflow; rescale units explicitly')
    low, high = np.quantile(means, [alpha, 1 - alpha])
    return {'estimate': float(delta.mean()), 'low': float(low), 'high': float(high), 'n': len(delta), 'direction': 'after minus before', 'confidence': confidence, 'interval': 'pointwise percentile bootstrap of paired subject mean differences', 'n_boot': n_boot, 'seed': int(seed), 'interpretation_note': ('Fewer than 10 independent subjects: interval stability/coverage require particular caution. ' if len(delta) < 10 else '') + 'Unadjusted; assumes independent subjects and complete pairs. No multiplicity correction or causal identification.'}

def inspect_canvas(fig, *, min_font_pt=6.5):
    if not np.isfinite(min_font_pt) or min_font_pt <= 0:
        raise ValueError('min_font_pt must be positive and finite')
    with warnings.catch_warnings(record=True) as emitted:
        warnings.simplefilter('always')
        fig.canvas.draw()
        fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = (fig.bbox.width, fig.bbox.height)
    outside, small, fonts = ([], [], set())
    excluded_text = set()
    for ax in fig.axes:
        if not ax.get_visible():
            excluded_text.update((id(text) for text in ax.findobj(match=Text)))
            continue
        for axis in (ax.xaxis, ax.yaxis):
            if not ax.axison or not axis.get_visible():
                excluded_text.update((id(text) for text in axis.findobj(match=Text)))
                continue
            low, high = sorted(axis.get_view_interval())
            tolerance = max(abs(high - low), 1) * 1e-10
            for tick in list(axis.get_major_ticks()) + list(axis.get_minor_ticks()):
                if not tick.get_visible() or tick.get_loc() < low - tolerance or tick.get_loc() > high + tolerance:
                    excluded_text.update((id(tick.label1), id(tick.label2)))
    for artist in fig.findobj(match=Text):
        if id(artist) in excluded_text or not artist.get_visible() or (not artist.get_text().strip()):
            continue
        axes = artist.axes
        if axes is not None and (not axes.get_visible()):
            continue
        box = artist.get_window_extent(renderer)
        if box.width == 0 or box.height == 0:
            continue
        if box.x0 < -1 or box.y0 < -1 or box.x1 > width + 1 or (box.y1 > height + 1):
            outside.append(artist.get_text())
        if artist.get_fontsize() < min_font_pt:
            small.append({'text': artist.get_text(), 'size_pt': artist.get_fontsize()})
        fonts.add(font_manager.findfont(artist.get_fontproperties()))
    return {'outside_canvas': outside, 'undersized_text': small, 'min_font_pt': min_font_pt, 'font_files': sorted(fonts), 'render_warnings': sorted(set((str(w.message) for w in emitted))), 'limitations': 'Does not verify all overlaps, non-text clipping, color accessibility or scientific meaning.'}

def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def export_figure(fig, stem, *, sources, panels, title, synthetic=False, dpi=300, formats=('png', 'svg'), min_font_pt=6.5, strict=True, overwrite=False):
    stem = Path(stem)
    if stem.suffix:
        raise ValueError('stem must not have a filename extension')
    if not isinstance(dpi, int) or dpi < 72:
        raise ValueError('dpi must be an integer >= 72')
    formats = tuple(formats)
    if not formats or len(set(formats)) != len(formats) or (not set(formats) <= {'png', 'svg', 'pdf'}):
        raise ValueError('formats must be distinct members of png, svg, pdf')
    if not str(title).strip() or not isinstance(synthetic, bool) or (not isinstance(overwrite, bool)):
        raise ValueError('provide a title and boolean synthetic/overwrite flags')
    if not panels:
        raise ValueError('provide panel-level scientific metadata')
    ids = []
    for panel in panels:
        for field in ('id', 'question', 'units', 'sample_unit', 'transform', 'estimate', 'interval'):
            if not isinstance(panel.get(field), str) or not panel[field].strip():
                raise ValueError(f'panel requires nonempty {field}')
        if isinstance(panel.get('n'), bool) or not isinstance(panel.get('n'), int) or panel['n'] < 1:
            raise ValueError('panel n must be a positive integer independent-unit count')
        ids.append(panel['id'])
    if len(set(ids)) != len(ids):
        raise ValueError('panel IDs must be unique')
    source_records = []
    for source in sources:
        path = Path(source).resolve(strict=True)
        if not path.is_file():
            raise ValueError(f'source is not a file: {path}')
        source_records.append({'path': str(path), 'sha256': _sha256(path), 'bytes': path.stat().st_size})
    if not source_records:
        raise ValueError('provide actual input/script files for provenance')
    targets = [stem.with_suffix('.' + fmt) for fmt in formats] + [stem.with_suffix('.manifest.json')]
    source_paths = {Path(item['path']) for item in source_records}
    for target in targets:
        if target.resolve() in source_paths:
            raise ValueError(f'output must never replace a listed source: {target}')
        if target.exists() and (not overwrite):
            raise FileExistsError(f'refusing to replace existing output: {target}')
    qa = inspect_canvas(fig, min_font_pt=min_font_pt)
    if strict and (qa['outside_canvas'] or qa['undersized_text'] or qa['render_warnings']):
        raise ValueError('figure QA failed: ' + json.dumps(qa, ensure_ascii=False))
    stem.parent.mkdir(parents=True, exist_ok=True)
    dimensions_mm = (fig.get_size_inches() * 25.4).tolist()
    manifest = {'schema_version': 1, 'title': title, 'synthetic': synthetic, 'size_mm': dimensions_mm, 'dpi': dpi, 'panels': panels, 'sources': source_records, 'runtime': {'python': platform.python_version(), 'matplotlib': mpl.__version__, 'numpy': np.__version__}, 'qa': qa, 'outputs': [], 'replacement_policy': 'All outputs render and validate before replacement. Each file replacement is atomic; a process/OS failure between files can leave a mixed export set.', 'semantic_validation': 'Panel metadata is supplied by caller; not verified against plotted data.'}
    with tempfile.TemporaryDirectory(prefix='.figure-export-', dir=stem.parent) as temporary:
        temp = Path(temporary)
        with mpl.rc_context({'svg.fonttype': 'none', 'savefig.bbox': None, 'pdf.fonttype': 42}):
            for fmt in formats:
                path = temp / (stem.name + '.' + fmt)
                fig.savefig(path, format=fmt, dpi=dpi, bbox_inches=None, facecolor='white')
                detail = {'path': str(stem.with_suffix('.' + fmt).resolve()), 'sha256': _sha256(path), 'bytes': path.stat().st_size}
                if fmt == 'svg':
                    xml = ET.parse(path).getroot()
                    detail['editable_text_elements'] = len(xml.findall('.//{http://www.w3.org/2000/svg}text'))
                    if detail['editable_text_elements'] == 0:
                        raise ValueError('SVG contains no editable text; inspect font settings')
                if fmt == 'png':
                    from PIL import Image
                    with Image.open(path) as img:
                        detail['pixels'] = list(img.size)
                    expected = np.array(fig.get_size_inches()) * dpi
                    if np.max(np.abs(expected - detail['pixels'])) > 1:
                        raise ValueError('PNG size differs from physical canvas contract')
                manifest['outputs'].append(detail)
        payload = json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)
        (temp / targets[-1].name).write_text(payload, encoding='utf-8')
        for target in targets:
            if target.exists() and (not overwrite):
                raise FileExistsError(f'output appeared during export: {target}')
        for target in targets:
            if overwrite:
                (temp / target.name).replace(target)
            else:
                (temp / target.name).rename(target)
    return manifest
