from . import vallado_cpp

class Satrec(vallado_cpp.Satrec):
    """High-speed computation of satellite positions and velocities."""

    __slots__ = ()

    def sgp4_array(self, jd, fr):
        """Compute positions and velocities for the times in a NumPy array.

        Given NumPy arrays ``jd`` and ``fr`` of the same length that
        supply the whole part and the fractional part of one or more
        Julian dates, return a tuple ``(e, r, v)`` of three vectors:

        * ``e``: nonzero for any dates that produced errors, 0 otherwise.
        * ``r``: position vectors in kilometers.
        * ``v``: velocity vectors in kilometers per second.

        """
        jd = jd.astype('float64', copy=False)
        fr = fr.astype('float64', copy=False)

        eshape = jd.shape
        fshape = eshape[0], 3

        array = type(jd)
        e = array(eshape, 'uint8')
        r = array(fshape, 'float64')
        v = array(fshape, 'float64')
        self._sgp4(jd, fr, e, r, v)
        return e, r, v

class SatrecArray(vallado_cpp.SatrecArray):
    """High-speed satellite array for computing positions and velocities."""

    __slots__ = ()

    def sgp4(self, jd, fr):
        """Compute positions and velocities for the satellites in this array.

        Given NumPy scalars or arrays ``jd`` and ``fr`` supplying the
        whole part and the fractional part of one or more Julian dates,
        return a tuple ``(e, r, v)`` of three vectors:

        * ``e``: nonzero for any dates that produced errors, 0 otherwise.
        * ``r``: position vectors in kilometers.
        * ``v``: velocity vectors in kilometers per second.

        The first dimension of each output vector has the same length as
        this satellite array, the second dimension the same length as
        the input date arrays, and the third dimension has length 3.

        """
        jd = jd.astype('float64', copy=False)
        fr = fr.astype('float64', copy=False)

        ilength = len(self)
        jlength, = jd.shape
        eshape = ilength, jlength
        fshape = ilength, jlength, 3

        array = type(jd)
        e = array(eshape, 'uint8')
        r = array(fshape, 'float64')
        v = array(fshape, 'float64')
        self._sgp4(jd, fr, e, r, v)
        return e, r, v
