#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "SGP4.h"
#include "structmember.h"

/* Whether scanf() prefers commas as the decimal point: true means that
   we cannot count on the current locale to interpret numbers correctly
   when the Vallado C++ code calls `sscanf()`. */

static bool switch_locale;

/* Satrec object that wraps a single raw SGP4 struct, and a Satrec array
   that can broadcast into NumPy arrays. */

typedef struct {
    PyObject_HEAD
    elsetrec satrec;
} SatrecObject;

typedef struct {
    PyObject_VAR_HEAD
    elsetrec satrec[0];
} SatrecArrayObject;

/* Support routine that is used to support NumPy array broadcasting for
   both individual satellite objects and also arrays. */

static PyObject *
_vectorized_sgp4(PyObject *args, elsetrec *raw_satrec_array, int imax)
{
    PyObject *jd_arg, *fr_arg, *e_arg, *r_arg, *v_arg;
    Py_buffer jd_buf, fr_buf, e_buf, r_buf, v_buf;
    PyObject *rv = NULL;

    // To prepare for "cleanup:" below.
    jd_buf.buf = fr_buf.buf = e_buf.buf = r_buf.buf = v_buf.buf = NULL;

    if (!PyArg_ParseTuple(args, "OOOOO:sgp4",
                          &jd_arg, &fr_arg, &e_arg, &r_arg, &v_arg))
        return NULL;

    if (PyObject_GetBuffer(jd_arg, &jd_buf, PyBUF_SIMPLE)) goto cleanup;
    if (PyObject_GetBuffer(fr_arg, &fr_buf, PyBUF_SIMPLE)) goto cleanup;
    if (PyObject_GetBuffer(e_arg, &e_buf, PyBUF_WRITABLE)) goto cleanup;
    if (PyObject_GetBuffer(r_arg, &r_buf, PyBUF_WRITABLE)) goto cleanup;
    if (PyObject_GetBuffer(v_arg, &v_buf, PyBUF_WRITABLE)) goto cleanup;

    if (jd_buf.len != fr_buf.len) {
        PyErr_SetString(PyExc_ValueError, "jd and fr must have the same shape");
        goto cleanup;
    }

    // This extra block allows the "goto" statements above to jump
    // across these further variable declarations.
    {
        Py_ssize_t jmax = jd_buf.len / sizeof(double);

        if ((r_buf.len != (Py_ssize_t) sizeof(double) * imax * jmax * 3) ||
            (v_buf.len != (Py_ssize_t) sizeof(double) * imax * jmax * 3) ||
            (e_buf.len != (Py_ssize_t) sizeof(uint8_t) * imax * jmax)) {
            PyErr_SetString(PyExc_ValueError, "bad output array dimension");
            goto cleanup;
        }

        double *jd = (double*) jd_buf.buf;
        double *fr = (double*) fr_buf.buf;
        double *r = (double*) r_buf.buf;
        double *v = (double*) v_buf.buf;
        uint8_t *e = (uint8_t*) e_buf.buf;

#pragma omp parallel for
        for (Py_ssize_t i=0; i < imax; i++) {
            elsetrec &satrec = raw_satrec_array[i];
            for (Py_ssize_t j=0; j < jmax; j++) {
                double t = (jd[j] - satrec.jdsatepoch) * 1440.0
                         + (fr[j] - satrec.jdsatepochF) * 1440.0;
                Py_ssize_t k1 = i * jmax + j;
                Py_ssize_t k3 = 3 * k1;
                SGP4Funcs::sgp4(satrec, t, r + k3, v + k3);
                e[k1] = (uint8_t) satrec.error;
                if (satrec.error && satrec.error < 6) {
                    r[k3] = r[k3+1] = r[k3+2] = v[k3] = v[k3+1] = v[k3+2] = NAN;
                }
            }
        }
    }

    Py_INCREF(Py_None);
    rv = Py_None;
 cleanup:
    if (jd_buf.buf) PyBuffer_Release(&jd_buf);
    if (fr_buf.buf) PyBuffer_Release(&fr_buf);
    if (r_buf.buf) PyBuffer_Release(&r_buf);
    if (v_buf.buf) PyBuffer_Release(&v_buf);
    if (e_buf.buf) PyBuffer_Release(&e_buf);
    return rv;
}

/* Details of the "Satrec" satellite object. */

static PyObject *
Satrec_twoline2rv(PyTypeObject *cls, PyObject *args)
{
    char *string1, *string2, line1[130], line2[130];
    gravconsttype whichconst = wgs72;
    double dummy;

    if (!PyArg_ParseTuple(args,
            "ss|i:twoline2rv",
            &string1, &string2, &whichconst))
        return NULL;

    // Copy both lines, since twoline2rv() might write to both buffers.
    strncpy(line1, string1, 130);
    strncpy(line2, string2, 130);

    // Truncate the line before any checksum characters, if present;
    // otherwise the C++ scanf() interprets each checksum character as
    // part of the preceding value.
    line1[68] = '\0';
    line2[68] = '\0';

    /* Allocate the new object. */
    SatrecObject *self = (SatrecObject*) cls->tp_alloc(cls, 0);
    if (!self)
        return NULL;

    /* Correct for locales that use a comma as the decimal point, since
       users report that the scanf() function on macOS is sensitive to
       locale when parsing floats.  This operation is not thread-safe,
       but we have not released the GIL. */
    float f;
    sscanf("1,5", "%f", &f);
    switch_locale = (f == 1.5);

    char *old_locale = NULL;
    if (switch_locale)
        old_locale = setlocale(LC_NUMERIC, "C");

    /* Leading spaces in a catalog number make scanf() in the official
       code consume the Classification letter as part of the catalog
       number.  (The first character of the International Designator
       then gets consumed as the Classification instead.)  But no
       parsing error is reported, which is bad for users, so let's avoid
       the situation by adding leading zeros ourselves. */
    for (int i=2; i<7; i++) {
        if (line1[i] == ' ') line1[i] = '0';
        if (line2[i] == ' ') line2[i] = '0';
    }

    /* Call the official routine. */
    SGP4Funcs::twoline2rv(line1, line2, ' ', ' ', 'i', whichconst,
                          dummy, dummy, dummy, self->satrec);

    /* Usability bonus: round the fractional day to exactly the eight
       digits specified in the TLE. */
    self->satrec.jdsatepochF = round(self->satrec.jdsatepochF * 1e8) / 1e8;

    /* To avoid having scanf() interpret the "intldesg" as zero or as
       several fields, the C++ library changes spaces to periods and
       underscores.  Let's convert them back to avoid surprising users
       and to match our Python implementation. */
    if (self->satrec.intldesg[0] == '.')
         self->satrec.intldesg[0] = ' ';
    for (int i=1; i<11; i++)
         if (self->satrec.intldesg[i] == '_')
              self->satrec.intldesg[i] = ' ';

    /* Restore previous locale. */
    if (switch_locale)
        setlocale(LC_NUMERIC, old_locale);

    return (PyObject*) self;
}

static PyObject *
Satrec_sgp4init(PyObject *self, PyObject *args)
{
    char satnum_str[6];
    int whichconst;  /* "int" rather than "gravconsttype" so we know size */
    int opsmode;     /* "int" rather than "char" because "C" needs an int */
    long int satnum;
    double epoch, bstar, ndot, nddot;
    double ecco, argpo, inclo, mo, no_kozai, nodeo;

    if (!PyArg_ParseTuple(args, "iCldddddddddd:sgp4init", &whichconst,
                          &opsmode, &satnum, &epoch, &bstar, &ndot, &nddot,
                          &ecco, &argpo, &inclo, &mo, &no_kozai, &nodeo))
        return NULL;

    // See https://www.space-track.org/documentation#tle-alpha5
    if (satnum < 100000) {
        snprintf(satnum_str, 6, "%ld", satnum);
    } else {
        char c = 'A' + satnum / 10000 - 10;
        if (c > 'I') c++;
        if (c > 'O') c++;
        satnum_str[0] = c;
        snprintf(satnum_str + 1, 5, "%04ld", satnum % 10000);
    }

    elsetrec &satrec = ((SatrecObject*) self)->satrec;

    SGP4Funcs::sgp4init((gravconsttype) whichconst, opsmode, satnum_str, epoch,
                        bstar, ndot, nddot, ecco, argpo, inclo, mo, no_kozai,
                        nodeo, satrec);

    /* Populate date fields that SGP4Funcs::twoline2rv would set. */
    int y, m, d, H, M;
    double S, jan0jd, jan0fr /* always comes out 0.0 */;
    SGP4Funcs::invjday_SGP4(2433281.5, epoch, y, m, d, H, M, S);
    SGP4Funcs::jday_SGP4(y, 1, 0, 0, 0, 0.0, jan0jd, jan0fr);
    satrec.epochyr = y % 1000;
    satrec.epochdays = 2433281.5 - jan0jd + epoch;
    satrec.jdsatepochF = modf(epoch, &satrec.jdsatepoch);
    satrec.jdsatepoch += 2433281.5;

    /* Return true as sgp4init does, satrec.error contains any error codes */

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
Satrec_sgp4_tsince(PyObject *self, PyObject *args)
{
    double tsince, r[3], v[3];

    if (!PyArg_ParseTuple(args, "d:sgp4_tsince", &tsince))
        return NULL;
    elsetrec &satrec = ((SatrecObject*) self)->satrec;
    SGP4Funcs::sgp4(satrec, tsince, r, v);
    if (satrec.error && satrec.error < 6)
        r[0] = r[1] = r[2] = v[0] = v[1] = v[2] = NAN;
    return Py_BuildValue("i(fff)(fff)", satrec.error,
                         r[0], r[1], r[2], v[0], v[1], v[2]);
}

static PyObject *
Satrec_sgp4(PyObject *self, PyObject *args)
{
    double jd, fr, r[3], v[3];
    if (!PyArg_ParseTuple(args, "dd:sgp4", &jd, &fr))
        return NULL;
    elsetrec &satrec = ((SatrecObject*) self)->satrec;
    double tsince = (jd - satrec.jdsatepoch) * 1440.0
                  + (fr - satrec.jdsatepochF) * 1440.0;
    SGP4Funcs::sgp4(satrec, tsince, r, v);
    if (satrec.error && satrec.error < 6)
        r[0] = r[1] = r[2] = v[0] = v[1] = v[2] = NAN;
    return Py_BuildValue("i(fff)(fff)", satrec.error,
                         r[0], r[1], r[2], v[0], v[1], v[2]);
}

static PyObject *
Satrec__sgp4(PyObject *self, PyObject *args)
{
    SatrecObject *obj = (SatrecObject*) self;
    elsetrec *raw_satrec_array = &(obj->satrec);
    Py_ssize_t imax = 1;
    return _vectorized_sgp4(args, raw_satrec_array, imax);
}

static PyMethodDef Satrec_methods[] = {
    {"twoline2rv", (PyCFunction)Satrec_twoline2rv, METH_VARARGS | METH_CLASS,
     PyDoc_STR("Initialize the record from two lines of TLE text and an optional gravity constant.")},
    {"sgp4init", (PyCFunction)Satrec_sgp4init, METH_VARARGS,
     PyDoc_STR("Initialize the record from orbital elements.")},
    {"sgp4", (PyCFunction)Satrec_sgp4, METH_VARARGS,
     PyDoc_STR("Given a modified julian date, return position and velocity.")},
    {"_sgp4", (PyCFunction)Satrec__sgp4, METH_VARARGS,
     PyDoc_STR("Given an array of modified julian dates, return position and velocity arrays.")},
    {"sgp4_tsince", (PyCFunction)Satrec_sgp4_tsince, METH_VARARGS,
     PyDoc_STR("Given minutes since epoch, return position and velocity.")},
    {NULL, NULL}
};

#define O(member) offsetof(SatrecObject, satrec.member)

static PyMemberDef Satrec_members[] = {
    /* Listed in the order they appear in a TLE record. */

    {"operationmode", T_CHAR, O(operationmode), READONLY,
     PyDoc_STR("Operation mode: 'a' legacy AFSPC, or 'i' improved.")},
    {"jdsatepoch", T_DOUBLE, O(jdsatepoch), 0,
     PyDoc_STR("Julian date of epoch, day number (see jdsatepochF).")},
    {"jdsatepochF", T_DOUBLE, O(jdsatepochF), 0,
     PyDoc_STR("Julian date of epoch, fraction of day (see jdsatepoch).")},
    {"classification", T_CHAR, O(classification), 0,
     "Usually U=Unclassified, C=Classified, or S=Secret."},
    /* intldesg: inline character array; see Satrec_getset. */
    {"epochyr", T_INT, O(epochyr), 0,
     PyDoc_STR("Year of this element set's epoch (see epochdays). Not set by sgp4init().")},
    {"epochdays", T_DOUBLE, O(epochdays), 0,
     PyDoc_STR("Day of the year of this element set's epoch (see epochyr). Not set by sgp4init().")},
    {"ndot", T_DOUBLE, O(ndot), READONLY,
     PyDoc_STR("Ballistic Coefficient in revs/day.")},
    {"nddot", T_DOUBLE, O(nddot), READONLY,
     PyDoc_STR("Second Derivative of Mean Motion in revs/day^3.")},
    {"bstar", T_DOUBLE, O(bstar), READONLY,
     PyDoc_STR("Drag Term in inverse Earth radii.")},
    {"ephtype", T_INT, O(ephtype), 0,
     PyDoc_STR("Ephemeris type (should be 0 in published TLEs).")},
    {"elnum", T_LONG, O(elnum), 0,
     PyDoc_STR("Element set number.")},
    {"inclo", T_DOUBLE, O(inclo), READONLY,
     PyDoc_STR("Inclination in radians.")},
    {"nodeo", T_DOUBLE, O(nodeo), READONLY,
     PyDoc_STR("Right ascension of ascending node in radians.")},
    {"ecco", T_DOUBLE, O(ecco), READONLY,
     PyDoc_STR("Eccentricity.")},
    {"argpo", T_DOUBLE, O(argpo), READONLY,
     PyDoc_STR("Argument of perigee in radians.")},
    {"mo", T_DOUBLE, O(mo), READONLY,
     PyDoc_STR("Mean anomaly in radians.")},
    {"no_kozai", T_DOUBLE, O(no_kozai), READONLY,
     PyDoc_STR("Mean motion in radians per minute.")},
    {"revnum", T_LONG, O(revnum), 0,
     PyDoc_STR("Integer revolution number at the epoch.")},

    /* For compatibility with the old struct members, also accept the
       plain name "no". */

    {"no", T_DOUBLE, O(no_kozai), READONLY,
     PyDoc_STR("Alias for the more carefully named ``no_kozai``.")},

    /* Derived values that do not appear explicitly in the TLE. */

    {"method", T_CHAR, O(method), READONLY,
     PyDoc_STR("Method, either 'n' near earth or 'd' deep space.")},
    {"error", T_INT, O(method), READONLY,
     PyDoc_STR("Error code (1-6) documented in sgp4()")},
    {"a", T_DOUBLE, O(a), READONLY,
     PyDoc_STR("semi-major axis")},
    {"altp", T_DOUBLE, O(altp), READONLY,
     PyDoc_STR("altitude of perigee")},
    {"alta", T_DOUBLE, O(alta), READONLY,
     PyDoc_STR("altitude of perigee")},

    /* Single averaged mean elements */

    {"am", T_DOUBLE, O(am), READONLY,
     PyDoc_STR("am: Average semi-major axis")},
    {"em", T_DOUBLE, O(em), READONLY,
     PyDoc_STR("em: Average eccentricity")},
    {"im", T_DOUBLE, O(im), READONLY,
     PyDoc_STR("im: Average inclination")},
    {"Om", T_DOUBLE, O(Om), READONLY,
     PyDoc_STR("Om: Average right ascension of ascending node")},
    {"om", T_DOUBLE, O(om), READONLY,
     PyDoc_STR("om: Average argument of perigee")},
    {"mm", T_DOUBLE, O(mm), READONLY,
     PyDoc_STR("mm: Average mean anomaly")},
    {"nm", T_DOUBLE, O(nm), READONLY,
     PyDoc_STR("nm: Average mean motion")},

    /* Gravity-constant dependent values (initialized by sgp4init() */

    {"tumin", T_DOUBLE, O(tumin), READONLY,
     PyDoc_STR("minutes in one time unit")},
    {"mu", T_DOUBLE, O(mus), READONLY,
     PyDoc_STR("Earth gravitational parameter")},
    {"radiusearthkm", T_DOUBLE, O(radiusearthkm), READONLY,
     PyDoc_STR("radius of the earth in km")},
    {"xke", T_DOUBLE, O(xke), READONLY,
     PyDoc_STR("reciprocal of tumin")},
    {"j2", T_DOUBLE, O(j2), READONLY,
     PyDoc_STR("un-normalized zonal harmonic j2 value")},
    {"j3", T_DOUBLE, O(j3), READONLY,
     PyDoc_STR("un-normalized zonal harmonic j3 value")},
    {"j4", T_DOUBLE, O(j4), READONLY,
     PyDoc_STR("un-normalized zonal harmonic j4 value")},
    {"j3oj2", T_DOUBLE, O(j3oj2), READONLY,
     PyDoc_STR("j3 divided by j2")},

    /* Other convenience variables (some required by propagation.py) */

    {"t", T_DOUBLE, O(t), READONLY,
     PyDoc_STR("Last tsince input to sgp4()")},
    {"mdot", T_DOUBLE, O(mdot), READONLY,
     PyDoc_STR("mean anomaly dot (rate)")},
    {"argpdot", T_DOUBLE, O(argpdot), READONLY,
     PyDoc_STR("argument of perigee dot (rate)")},
    {"nodedot", T_DOUBLE, O(nodedot), READONLY,
     PyDoc_STR("right ascension of ascending node dot (rate)")},
    {"gsto", T_DOUBLE, O(gsto), READONLY,
     PyDoc_STR("gsto: greenwich sidereal time")},

    {NULL}
};

#undef O

static PyObject *
get_intldesg(SatrecObject *self, void *closure)
{
    int length = 0;
    char *s = self->satrec.intldesg;
    while (length <= 7 && s[length])
        length++;
    while (length && s[length-1] == ' ')
        length--;
    return PyUnicode_FromStringAndSize(s, length);
}

static int
set_intldesg(SatrecObject *self, PyObject *value, void *closure)
{
    if (!PyUnicode_Check(value))
        return -1;
    const char *s = PyUnicode_AsUTF8(value);
    if (!s)
        return -1;
    strncpy(self->satrec.intldesg, s, 11);
    return 0;
}

static PyObject *
get_satnum(SatrecObject *self, void *closure)
{
    char *s = self->satrec.satnum;
    long n;
    if (strlen(s) < 5 || s[0] <= '9')
        n = atol(s);
    else if (s[0] <= 'I')
        n = (s[0] - 'A' + 10) * 10000 + atol(s + 1);
    else if (s[0] <= 'O')
        n = (s[0] - 'A' + 9) * 10000 + atol(s + 1);
    else
        n = (s[0] - 'A' + 8) * 10000 + atol(s + 1);
    return PyLong_FromLong(n);
}

static PyGetSetDef Satrec_getset[] = {
    {"intldesg", (getter)get_intldesg, (setter)set_intldesg,
     PyDoc_STR("International Designator: a string of up to 7 characters"
               " from the first line of the TLE that typically provides"
               " two digits for the launch year, a 3-digit launch number,"
               " and one or two letters for which piece of the launch.")},
    {"satnum", (getter)get_satnum, NULL,
     PyDoc_STR("Satellite number, from characters 3-7 of each TLE line.")},
    {NULL},
};

static PyTypeObject SatrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    /* See the module initialization function at the bottom of this file. */
};

/* Details of the SatrecArray. */

static Py_ssize_t
Satrec_len(PyObject *self) {
    return ((SatrecArrayObject*)self)->ob_base.ob_size;
}

static PySequenceMethods SatrecArray_as_sequence = {
    Satrec_len
};

static PyObject *
SatrecArray_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *sequence;
    if (!PyArg_ParseTuple(args, "O:SatrecArray", &sequence))
        return NULL;

    Py_ssize_t length = PySequence_Length(sequence);
    if (length == -1)
        return NULL;

    return type->tp_alloc(type, length);
}


static int
SatrecArray_init(SatrecArrayObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *sequence;
    if (!PyArg_ParseTuple(args, "O:SatrecArray", &sequence))
        return -1;

    Py_ssize_t length = PySequence_Length(sequence);
    if (length == -1)
        return -1;

    for (Py_ssize_t i=0; i<length; i++) {
        PyObject *item = PySequence_GetItem(sequence, i);
        if (!item)
            return -1;
        if (!PyObject_IsInstance(item, (PyObject*) &SatrecType)) {
            PyErr_Format(PyExc_ValueError, "every item must be a Satrec,"
                         " but element %d is: %R", i, item);
            Py_DECREF(item);
            return -1;
        }
        self->satrec[i] = ((SatrecObject *)item)->satrec;
        Py_DECREF(item);
    }
    return 0;
}

static PyObject *
SatrecArray_sgp4(PyObject *self, PyObject *args)
{
    SatrecArrayObject *satrec_array = (SatrecArrayObject*) self;
    elsetrec *raw_satrec_array = &(satrec_array->satrec[0]);
    Py_ssize_t imax = satrec_array->ob_base.ob_size;
    return _vectorized_sgp4(args, raw_satrec_array, imax);
}

static PyMethodDef SatrecArray_methods[] = {
    {"_sgp4", (PyCFunction)SatrecArray_sgp4, METH_VARARGS,
     PyDoc_STR("Given array of minutes since epoch, write"
               " positions, velocities, and errors to arrays.")},
    {NULL, NULL}
};

static PyTypeObject SatrecArrayType = {
    PyVarObject_HEAD_INIT(NULL, sizeof(elsetrec))
    /* See the module initialization function at the bottom of this file. */
};

/* The module that ties it all together. */
static PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "sgp4.vallado_cpp",
    "Official C++ SGP4 implementation.",
    -1
};

#include <stdio.h>

PyMODINIT_FUNC
PyInit_vallado_cpp(void)
{
    SatrecType.tp_name = "sgp4.vallado_cpp.Satrec";
    SatrecType.tp_basicsize = sizeof(SatrecObject);
    SatrecArrayType.tp_itemsize = 0;
    SatrecType.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    SatrecType.tp_doc = "SGP4 satellite record.";
    SatrecType.tp_methods = Satrec_methods;
    SatrecType.tp_members = Satrec_members;
    SatrecType.tp_getset = Satrec_getset;
    SatrecType.tp_new = PyType_GenericNew;

    if (PyType_Ready(&SatrecType) < 0)
        return NULL;

    SatrecArrayType.tp_name = "sgp4.vallado_cpp.SatrecArray";
    SatrecArrayType.tp_basicsize = sizeof(SatrecArrayObject);
    SatrecArrayType.tp_itemsize = sizeof(elsetrec);
    SatrecArrayType.tp_as_sequence = &SatrecArray_as_sequence;
    SatrecArrayType.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    SatrecArrayType.tp_doc = "SGP4 array of satellites.";
    SatrecArrayType.tp_methods = SatrecArray_methods;
    SatrecArrayType.tp_init = (initproc) SatrecArray_init;
    SatrecArrayType.tp_new = SatrecArray_new;

    if (PyType_Ready(&SatrecArrayType) < 0)
        return NULL;

    PyObject *m = PyModule_Create(&module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&SatrecType);
    if (PyModule_AddObject(m, "Satrec", (PyObject *) &SatrecType) < 0) {
        Py_DECREF(&SatrecType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&SatrecArrayType);
    if (PyModule_AddObject(m, "SatrecArray", (PyObject *) &SatrecArrayType) < 0) {
        Py_DECREF(&SatrecArrayType);
        Py_DECREF(&SatrecType);
        Py_DECREF(m);
        return NULL;
    }

    if (PyModule_AddIntConstant(m, "WGS72", (int)wgs72) ||
        PyModule_AddIntConstant(m, "WGS72OLD", (int)wgs72old) ||
        PyModule_AddIntConstant(m, "WGS84", (int)wgs84))
        return NULL;

    return m;
}
