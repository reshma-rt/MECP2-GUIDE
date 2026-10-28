"""Frozen utilities derived from the tested ranking adapter; no ranking routine."""
import contextlib, csv, ctypes, datetime as dt, hashlib, io, json, math, os, re, sys
from pathlib import Path
from collections import Counter
MUTATIONS = {'143546', '143559', '143597', '2505626', '3237211', '372650', '4279604', '548706'}

def digest(data):
    return hashlib.sha256(data).hexdigest()

def require(condition, message):
    if not condition:
        raise ValueError(message)

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='microseconds')

def strict_json(data):

    def unique(pairs):
        result = {}
        for k, v in pairs:
            require(k not in result, 'Duplicate JSON key: ' + k)
            result[k] = v
        return result
    return json.loads(data, object_pairs_hook=unique, parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON: ' + x)))

def table(data, required=None):
    reader = csv.DictReader(io.StringIO(data.decode('utf-8-sig'), newline=''))
    fields = reader.fieldnames or []
    require(len(fields) == len(set(fields)), 'Duplicate/ambiguous CSV columns')
    require(all(fields), 'Empty header')
    if required is not None:
        require(fields == required, 'Exact schema mismatch; no aliases or fallback allowed')
    rows = list(reader)
    require(all((None not in r and all((v is not None for v in r.values())) for r in rows)), 'Malformed CSV row')
    return (rows, fields)

def number(value, unit=False):
    x = float(value)
    require(math.isfinite(x), 'Nonfinite input')
    require(not unit or 0 <= x <= 1, 'Value outside [0,1]; no clamping')
    return x

def integer(value):
    require(isinstance(value, str) and bool(re.fullmatch('\\d+', value)), 'Expected nonnegative integer string')
    return int(value)

def indexed(rows):
    require(all((r.get('Guide_ID') for r in rows)), 'Missing Guide_ID')
    require(len(rows) == len({r['Guide_ID'] for r in rows}), 'Duplicate Guide_ID')
    return {r['Guide_ID']: r for r in rows}

def cohort(rows):
    ids = indexed(rows)
    require(set(ids) == {'T1G%04d' % i for i in range(1, 377)}, 'Incomplete/extra frozen guide cohort')
    counts = Counter((r['VariationID'] for r in rows))
    require(set(counts) == MUTATIONS and set(counts.values()) == {47}, 'Incomplete/mismatched mutation groups')
    return ids

@contextlib.contextmanager
def locked_input(path):
    require(path.is_file() and (not path.is_symlink()), 'Input unavailable or symlink: ' + str(path))
    for ancestor in [path, *path.parents]:
        st = ancestor.lstat()
        require(not getattr(st, 'st_file_attributes', 0) & 1024, 'Reparse path rejected: ' + str(ancestor))
    if os.name == 'nt':
        import msvcrt
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.CreateFileW.argtypes = [ctypes.c_wchar_p, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_void_p, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_void_p]
        kernel.CreateFileW.restype = ctypes.c_void_p
        handle = kernel.CreateFileW(str(path), 2147483648, 1, None, 3, 128, None)
        require(handle != ctypes.c_void_p(-1).value, 'Cannot lock read-only input: ' + str(path))
        fd = msvcrt.open_osfhandle(handle, os.O_RDONLY | os.O_BINARY)
        stream = os.fdopen(fd, 'rb')
    else:
        stream = path.open('rb')
    with stream:
        yield stream

def write_json(path, obj):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')

def write_csv(path, rows, fields):
    with path.open('x', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='raise')
        w.writeheader()
        w.writerows(rows)
