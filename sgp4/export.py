"""
Export orbit data to a Two-Line-Element representation.
"""
from sgp4.io import compute_checksum

from sgp4.model import Satrec


def export_tle(satrec):
    """
     Generate the TLE for the data in the `Satrec` object.

    Forms the two string lines of the TLE and returns them as a Tuple.

    Parameters
    ----------
    satrec : Satrec
        Object that holds the orbit information

    Returns
    -------
    (line1, line2) : (str, str)
    """

    # Define constants
    from math import pi
    deg2rad = pi / 180.0  # 0.0174532925199433
    xpdotp = 1440.0 / (2.0 * pi)  # 229.1831180523293

    # --------------------- Start generating line 1 ---------------------

    # Build the list by appending successive items
    line_as_list = ["1 "]

    # Pad the `satnum` entry with zeros
    if len(str(satrec.satnum)) <= 5:
        line_as_list.append(str(satrec.satnum).zfill(5))

    # Add classification code (use "U" if empty)
    line_as_list.append((satrec.classification.strip() or "U") + " ")

    # Add int'l designator and pad to 8 chars
    line_as_list.append(satrec.intldesg.ljust(8, " ") + " ")

    # Add epoch year and days in YYDDD.DDDDDDDD format
    line_as_list.append(str(satrec.epochyr).zfill(2) + "{:012.8f}".format(satrec.epochdays) + " ")

    # Add First Time Derivative of the Mean Motion (don't use "+")
    line_as_list.append("{0: 8.8f}".format(satrec.ndot * (xpdotp * 1440.0)).replace("0", "", 1) + " ")

    # Add Second Time Derivative of Mean Motion (don't use "+")
    # Multiplication with 10 is a hack to get the exponent right
    line_as_list.append("{0: 4.4e}".format((satrec.nddot * (xpdotp * 1440.0 * 1440)) * 10).replace(".", "")
                        .replace("e+00", "-0").replace("e-0", "-") + " ")

    # Add BSTAR
    # Multiplication with 10 is a hack to get the exponent right
    line_as_list.append("{0: 4.4e}".format(satrec.bstar * 10).replace(".", "").replace("e+00", "+0").replace("e-0", "-") + " ")

    # Add Ephemeris Type and Element Number
    line_as_list.append("{} ".format(satrec.ephtype) + str(satrec.elnum).rjust(4, " "))

    # Join all the parts and add the Checksum
    line1 = ''.join(line_as_list)
    line1 += str(compute_checksum(line1))

    # --------------------- Start generating line 2 ---------------------

    # Reset the str array
    line_as_list = ["2 "]

    # Pad the `satnum` entry with zeros
    if len(str(satrec.satnum)) <= 5:
        line_as_list.append(str(satrec.satnum).zfill(5) + " ")

    # Add the inclination (deg)
    line_as_list.append("{0:8.4f}".format(satrec.inclo / deg2rad).rjust(8, " ") + " ")

    # Add the RAAN (deg)
    line_as_list.append("{0:8.4f}".format(satrec.nodeo / deg2rad).rjust(8, " ") + " ")

    # Add the eccentricity (delete the leading zero an decimal point)
    line_as_list.append("{0:8.7f}".format(satrec.ecco).replace("0.", "") + " ")

    # Add the Argument of Perigee (deg)
    line_as_list.append("{0:8.4f}".format(satrec.argpo / deg2rad).rjust(8, " ") + " ")

    # Add the Mean Anomaly (deg)
    line_as_list.append("{0:8.4f}".format(satrec.mo / deg2rad).rjust(8, " ") + " ")

    # Add the Mean Motion (revs/day)
    line_as_list.append("{0:11.8f}".format(satrec.no_kozai * xpdotp).rjust(8, " "))

    # Add the rev number at epoch
    line_as_list.append(str(satrec.revnum).zfill(4))

    # Join all the parts and add the Checksum
    line2 = ''.join(line_as_list)
    line2 += str(compute_checksum(line2))

    return line1, line2
