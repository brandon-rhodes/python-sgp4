"""Export orbit data to a Two-Line-Element representation.

Contributed by Egemen Imre https://github.com/egemenimre

"""
from math import pi
from sgp4.io import compute_checksum

# Define constants
_deg2rad = pi / 180.0  # 0.0174532925199433
_xpdotp = 1440.0 / (2.0 * pi)  # 229.1831180523293

def export_tle(satrec):
    """Generate the TLE for a given `Satrec` object; returns two strings."""

    # --------------------- Start generating line 1 ---------------------

    # Build the list by appending successive items
    pieces = ["1 "]
    append = pieces.append

    # Pad the `satnum` entry with zeros
    if len(str(satrec.satnum)) <= 5:
        append(str(satrec.satnum).zfill(5))

    # Add classification code (use "U" if empty)
    append((satrec.classification.strip() or "U") + " ")

    # Add int'l designator and pad to 8 chars
    append(satrec.intldesg.ljust(8, " ") + " ")

    # Add epoch year and days in YYDDD.DDDDDDDD format
    append(str(satrec.epochyr).zfill(2) + "{:012.8f}".format(satrec.epochdays) + " ")

    # Add First Time Derivative of the Mean Motion (don't use "+")
    append("{0: 8.8f}".format(satrec.ndot * (_xpdotp * 1440.0)).replace("0", "", 1) + " ")

    # Add Second Time Derivative of Mean Motion (don't use "+")
    # Multiplication with 10 is a hack to get the exponent right
    append("{0: 4.4e}".format((satrec.nddot * (_xpdotp * 1440.0 * 1440)) * 10).replace(".", "")
                        .replace("e+00", "-0").replace("e-0", "-") + " ")

    # Add BSTAR
    # Multiplication with 10 is a hack to get the exponent right
    append("{0: 4.4e}".format(satrec.bstar * 10).replace(".", "").replace("e+00", "+0").replace("e-0", "-") + " ")

    # Add Ephemeris Type and Element Number
    append("{} ".format(satrec.ephtype) + str(satrec.elnum).rjust(4, " "))

    # Join all the parts and add the Checksum
    line1 = ''.join(pieces)
    line1 += str(compute_checksum(line1))

    # --------------------- Start generating line 2 ---------------------

    # Reset the str array
    pieces = ["2 "]
    append = pieces.append

    # Pad the `satnum` entry with zeros
    if len(str(satrec.satnum)) <= 5:
        append(str(satrec.satnum).zfill(5) + " ")

    # Add the inclination (deg)
    append("{0:8.4f}".format(satrec.inclo / _deg2rad).rjust(8, " ") + " ")

    # Add the RAAN (deg)
    append("{0:8.4f}".format(satrec.nodeo / _deg2rad).rjust(8, " ") + " ")

    # Add the eccentricity (delete the leading zero an decimal point)
    append("{0:8.7f}".format(satrec.ecco).replace("0.", "") + " ")

    # Add the Argument of Perigee (deg)
    append("{0:8.4f}".format(satrec.argpo / _deg2rad).rjust(8, " ") + " ")

    # Add the Mean Anomaly (deg)
    append("{0:8.4f}".format(satrec.mo / _deg2rad).rjust(8, " ") + " ")

    # Add the Mean Motion (revs/day)
    append("{0:11.8f}".format(satrec.no_kozai * _xpdotp).rjust(8, " "))

    # Add the rev number at epoch
    append(str(satrec.revnum).rjust(5))

    # Join all the parts and add the Checksum
    line2 = ''.join(pieces)
    line2 += str(compute_checksum(line2))

    return line1, line2
