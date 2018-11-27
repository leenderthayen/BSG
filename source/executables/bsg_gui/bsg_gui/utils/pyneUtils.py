time_conv_dict = {'as': 1e-18,
                  'attosec': 1e-18,
                  'attosecond': 1e-18,
                  'attoseconds': 1e-18,
                  'fs': 1e-15,
                  'femtosec': 1e-15,
                  'femtosecond': 1e-15,
                  'femtoseconds': 1e-15,
                  'ps': 1e-12,
                  'picosec': 1e-12,
                  'picosecond': 1e-12,
                  'picoseconds': 1e-12,
                  'ns': 1e-9,
                  'nanosec': 1e-9,
                  'nanosecond': 1e-9,
                  'nanoseconds': 1e-9,
                  'us': 1e-6,
                  'microsec': 1e-6,
                  'microsecond': 1e-6,
                  'microseconds': 1e-6,
                  'ms': 1e-3,
                  'millisec': 1e-3,
                  'millisecond': 1e-3,
                  'milliseconds': 1e-3,
                  's': 1.0,
                  'sec': 1.0,
                  'second': 1.0,
                  'seconds': 1.0,
                  'm': 60.0,
                  'min': 60.0,
                  'minute': 60.0,
                  'minutes': 60.0,
                  'h': 3600.0,
                  'hour': 3600.0,
                  'hours': 3600.0,
                  'd': 86400.0,
                  'day': 86400.0,
                  'days': 86400.0,
                  'w': 86400.0*7.0,
                  'week': 86400.0*7.0,
                  'weeks': 86400.0*7.0,
                  'y': 86400.0*365.25,
                  'year': 86400.0*365.25,
                  'years': 86400.0*365.25,
                  'c': 86400.0*365.25*100,
                  'century': 86400.0*365.25*100,
                  'centuries': 86400.0*365.25*100,
                  }

def to_id(nuc):
    nuc = nuc.strip()
    n = _nucleus.match(nuc)
    nucid = 0
    if n:
        A = int(n.group(1))
        Z = int(atoms.index(n.group(2).capitalize()))+1
        nucid = (1000 * Z + A) * 1000
    return nucid

def id_from_level(nuc, level, levellist, special=' '):
    nostate = int(nuc / 1000) * 1000
    if not levellist:
        return None
    ret_id = nuc
    minDif = 1e9
    for l in levellist:
        if int(l[0] / 1000 ) * 1000 == nostate:
             if abs(l[5] - level) < minDif and special == l[-1]:
                 minDif = abs(l[5] - level)
                 ret_id = l[0]
    return ret_id
