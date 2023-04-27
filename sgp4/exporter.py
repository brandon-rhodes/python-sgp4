"""Export orbit data to a Two-Line-Element representation.

Contributed by Egemen Imre https://github.com/egemenimre

"""
from math import pi
from sgp4.io import compute_checksum
from sgp4.conveniences import sat_epoch_datetime

# Define constants
_deg2rad = pi / 180.0  # 0.0174532925199433
_xpdotp = 1440.0 / (2.0 * pi)  # 229.1831180523293

def export_tle(satrec):
    """Generate the TLE for a given `Satrec` object; returns two strings."""

    # --------------------- Start generating line 1 ---------------------

    # Build the list by appending successive items
    pieces = ['1 ', satrec.satnum_str]
    append = pieces.append

    # Add classification code (use "U" if empty)
    classification = getattr(satrec, 'classification', 'U')
    append(classification.strip() or 'U')
    append(' ')

    # Add int'l designator and pad to 8 chars
    intldesg = getattr(satrec, 'intldesg', '')
    append('{0:8} '.format(intldesg))

    # Add epoch year and days in YYDDD.DDDDDDDD format
    epochyr = satrec.epochyr
    # Undo non-standard 4-digit year for old satrec objects
    epochyr %= 100
    append(str(epochyr).zfill(2) + "{:012.8f}".format(satrec.epochdays) + " ")

    # Add First Time Derivative of the Mean Motion (don't use "+")
    append("{0: 8.8f}".format(satrec.ndot * (_xpdotp * 1440.0)).replace("0", "", 1) + " ")

    # Add Second Time Derivative of Mean Motion
    append(_abbreviate_rate(satrec.nddot * _xpdotp * 20736000.0, '-0'))

    # Add BSTAR
    append(_abbreviate_rate(satrec.bstar * 10.0, '+0'))

    # Add Ephemeris Type and Element Number
    ephtype = getattr(satrec, 'ephtype', 0)
    elnum = getattr(satrec, 'elnum', 0)
    append('{0} {1:4}'.format(ephtype, elnum))

    # Join all the parts and add the Checksum
    line1 = ''.join(pieces)
    line1 += str(compute_checksum(line1))

    # --------------------- Start generating line 2 ---------------------

    # Reset the str array
    pieces = ['2 ', satrec.satnum_str]
    append = pieces.append

    # Add the inclination (deg)
    if not 0 <= satrec.inclo <= pi:
        raise ValueError("Inclination must be between 0 and pi, got %r", satrec.inclo)
    append(' {0:8.4f} '.format(satrec.inclo / _deg2rad))

    # Add the RAAN (deg)
    if not 0 <= satrec.nodeo <= 2 * pi:
        raise ValueError("RAAN must be between 0 and 2 pi, got %r", satrec.nodeo)
    append("{0:8.4f}".format(satrec.nodeo / _deg2rad).rjust(8, " ") + " ")

    # Add the eccentricity (delete the leading zero an decimal point)
    append("{0:8.7f}".format(satrec.ecco).replace("0.", "") + " ")

    # Add the Argument of Perigee (deg)
    if not 0 <= satrec.argpo <= 2 * pi:
        raise ValueError("Argument of Perigee must be between 0 and 2 pi, got %r", satrec.argpo)
    append("{0:8.4f}".format(satrec.argpo / _deg2rad).rjust(8, " ") + " ")

    # Add the Mean Anomaly (deg)
    if not 0 <= satrec.mo <= 2 * pi:
        raise ValueError("Mean Anomaly must be between 0 and 2 pi, got %r", satrec.mo)
    append("{0:8.4f}".format(satrec.mo / _deg2rad).rjust(8, " ") + " ")

    # Add the Mean Motion (revs/day)
    append("{0:11.8f}".format(satrec.no_kozai * _xpdotp).rjust(8, " "))

    # Add the rev number at epoch
    append(str(satrec.revnum).rjust(5))

    # Join all the parts and add the Checksum
    line2 = ''.join(pieces)
    line2 += str(compute_checksum(line2))

    return line1, line2

def _abbreviate_rate(value, zero_exponent_string):
    return (
        '{0: 4.4e} '.format(value)
        .replace('.', '')
        .replace('e+00', zero_exponent_string)
        .replace('e-0', '-')
        .replace('e+0', '+')
    )

def export_omm(satrec, object_name):
    launch_year = int(satrec.intldesg[:2])
    launch_year += 1900 + (launch_year < 57) * 100
    object_id = '{0}-{1}'.format(launch_year, satrec.intldesg[2:])

    return {
        "OBJECT_NAME": object_name,
        "OBJECT_ID": object_id,
        "CENTER_NAME": "EARTH",
        "REF_FRAME": "TEME",
        "TIME_SYSTEM": "UTC",
        "MEAN_ELEMENT_THEORY": "SGP4",
        "EPOCH": sat_epoch_datetime(satrec).strftime('%Y-%m-%dT%H:%M:%S.%f'),
        "MEAN_MOTION": satrec.no_kozai * _xpdotp,
        "ECCENTRICITY": satrec.ecco,
        "INCLINATION": satrec.inclo / _deg2rad,
        "RA_OF_ASC_NODE": satrec.nodeo / _deg2rad,
        "ARG_OF_PERICENTER": satrec.argpo / _deg2rad,
        "MEAN_ANOMALY": satrec.mo / _deg2rad,
        "EPHEMERIS_TYPE": satrec.ephtype,
        "CLASSIFICATION_TYPE": satrec.classification,
        "NORAD_CAT_ID": satrec.satnum,
        "ELEMENT_SET_NO": satrec.elnum,
        "REV_AT_EPOCH": satrec.revnum,
        "BSTAR": satrec.bstar,
        "MEAN_MOTION_DOT": satrec.ndot * (_xpdotp * 1440.0),
        "MEAN_MOTION_DDOT": satrec.nddot * (_xpdotp * 1440.0 * 1440),
    }
