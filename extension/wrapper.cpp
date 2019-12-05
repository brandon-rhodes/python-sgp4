#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "SGP4.h"
#include "structmember.h"

/* Satrec object that wraps a single raw SGP4 struct. */

typedef struct {
    PyObject_HEAD
    elsetrec satrec;
} SatrecObject;

static PyObject *
Satrec_twoline2rv(PyObject *cls, PyObject *args)
{
    char *string1, *string2, line1[130], line2[130];
    double dummy;

    if (!PyArg_ParseTuple(args, "ss:twoline2rv", &string1, &string2))
        return NULL;

    // Copy both lines, since twoline2rv() might write to both buffers.
    strncpy(line1, string1, 130);
    strncpy(line2, string2, 130);
    line1[129] = '\0';
    line2[129] = '\0';

    PyTypeObject *type = (PyTypeObject*) cls;
    SatrecObject *self = (SatrecObject*) PyObject_New(SatrecObject, type);
    if (!self)
        return NULL;

    SGP4Funcs::twoline2rv(line1, line2, ' ', ' ', 'i', wgs72,
                          dummy, dummy, dummy, self->satrec);

    return (PyObject*) self;
}

static PyObject *
Satrec_sgp4(PyObject *self, PyObject *args)
{
    double tsince, r[3], v[3];
    if (!PyArg_ParseTuple(args, "d:sgp4", &tsince))
        return NULL;
    elsetrec &satrec = ((SatrecObject *)self)->satrec;
    SGP4Funcs::sgp4(satrec, tsince, r, v);
    if (satrec.error && satrec.error < 6)
        r[0] = r[1] = r[2] = v[0] = v[1] = v[2] = NAN;
    return Py_BuildValue("(fff)(fff)i", r[0], r[1], r[2], v[0], v[1], v[2],
                         satrec.error);
}

static PyMethodDef Satrec_methods[] = {
    {"twoline2rv", (PyCFunction)Satrec_twoline2rv, METH_VARARGS | METH_CLASS,
     PyDoc_STR("Initialize the record from two lines of TLE text.")},
    {"sgp4", (PyCFunction)Satrec_sgp4, METH_VARARGS,
     PyDoc_STR("Given minutes since epoch, return a position and velocity.")},
    {NULL, NULL}
};

#define O(member) offsetof(SatrecObject, satrec.member)

static PyMemberDef Satrec_members[] = {
    /* Listed in the order they appear in a TLE record. */

    {"satnum", T_INT, O(satnum), READONLY,
     PyDoc_STR("Satellite number, from characters 3-7 of each TLE line.")},
    {"classification", T_CHAR, O(classification), READONLY,
     "Usually U=Unclassified, C=Classified, or S=Secret."},
    /* intldesg: inline character array; see Satrec_getset. */
    {"epochyr", T_INT, O(epochyr), READONLY,
     PyDoc_STR("Year of this element set's epoch (see epochdays).")},
    {"epochdays", T_DOUBLE, O(epochdays), READONLY,
     PyDoc_STR("Day of the year of this element set's epoch (see epochyr).")},
    {"ndot", T_DOUBLE, O(ndot), READONLY,
     PyDoc_STR("Ballistic Coefficient in revs/day.")},
    {"nddot", T_DOUBLE, O(nddot), READONLY,
     PyDoc_STR("Second Derivative of Mean Motion in revs/day^3.")},
    {"bstar", T_DOUBLE, O(bstar), READONLY,
     PyDoc_STR("Drag Term in inverse Earth radii.")},
    {"ephtype", T_INT, O(ephtype), READONLY,
     PyDoc_STR("Ephemeris type (should be 0 in published TLEs).")},
    {"elnum", T_INT, O(elnum), READONLY,
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
    {"revnum", T_LONG, O(revnum), READONLY,
     PyDoc_STR("Integer revolution number at the epoch.")},

    /* Derived values that do not appear explicitly in the TLE. */

    {"method", T_CHAR, O(method), READONLY,
     PyDoc_STR("Method, either 'n' near earth or 'd' deep space.")},
    {NULL}
};

#undef O

static PyObject *
get_intldesg(SatrecObject *self, void *closure)
{
  return PyUnicode_FromStringAndSize(self->satrec.intldesg, 6);
}

static PyGetSetDef Satrec_getset[] = {
    {"intldesg", (getter)get_intldesg, NULL,
     PyDoc_STR("A string copied from the first TLE line that typically"
               " provides the last two digits of the launch year, a"
               " 3-digit launch number, and a letter A-Z specifying"
               " which piece of the launch.")},
    {NULL},
};

static PyTypeObject SatrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    /* See the module initialization function at the bottom of this file. */
};

/* Satrec array that can broadcast into NumPy arrays. */

typedef struct {
    PyObject_VAR_HEAD
    elsetrec satrec[0];
} SatrecArrayObject;

static Py_ssize_t
Satrec_len(PyObject *self) {
    return ((SatrecArrayObject*)self)->ob_base.ob_size;
}

static PySequenceMethods SatrecArray_as_sequence = {
    sq_length : Satrec_len,
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

/* This, and the definition for which it is a forward reference,
   originally said "static PyTypeObject SatrecArrayType;" but that gave
   "error: redefinition" during compilation, thus the switch to API. */
PyAPI_DATA(PyTypeObject) SatrecArrayType;

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
        if (item->ob_type != &SatrecType) {
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
    PyObject *jd_arg, *fr_arg, *r_arg, *v_arg, *e_arg;
    Py_buffer jd_buf, fr_buf, r_buf, v_buf, e_buf;
    PyObject *rv = NULL;

    // To prepare for "cleanup:" below.
    jd_buf.buf = fr_buf.buf = r_buf.buf = v_buf.buf = e_buf.buf = NULL;

    if (!PyArg_ParseTuple(args, "OOOOO:sgp4",
                          &jd_arg, &fr_arg, &r_arg, &v_arg, &e_arg))
        return NULL;

    if (PyObject_GetBuffer(jd_arg, &jd_buf, PyBUF_SIMPLE)) goto cleanup;
    if (PyObject_GetBuffer(fr_arg, &fr_buf, PyBUF_SIMPLE)) goto cleanup;
    if (PyObject_GetBuffer(r_arg, &r_buf, PyBUF_WRITABLE)) goto cleanup;
    if (PyObject_GetBuffer(v_arg, &v_buf, PyBUF_WRITABLE)) goto cleanup;
    if (PyObject_GetBuffer(e_arg, &e_buf, PyBUF_WRITABLE)) goto cleanup;

    if (jd_buf.len != fr_buf.len) {
        PyErr_SetString(PyExc_ValueError, "jd and fr must have the same shape");
        goto cleanup;
    }

    // This extra block allows the "goto" statements above to jump
    // across these further variable declarations.
    {
        SatrecArrayObject *satrec_array = (SatrecArrayObject*) self;
        Py_ssize_t imax = ((SatrecArrayObject*) self)->ob_base.ob_size;
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
            elsetrec *satrec = satrec_array->satrec + i;
            for (Py_ssize_t j=0; j < jmax; j++) {
                double t = (jd[j] - satrec->jdsatepoch) * 1440.0
                         + (fr[j] - satrec->jdsatepochF) * 1440.0;
                Py_ssize_t k = i * jmax + j;
                SGP4Funcs::sgp4(*satrec, t, r + k*3, v + k*3);
                if (satrec->error && satrec->error < 6) {
                    r[k] = r[k+1] = r[k+2] = v[k] = v[k+1] = v[k+2] = NAN;
                }
                e[k] = (uint8_t) satrec->error;
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

static PyMethodDef SatrecArray_methods[] = {
    {"sgp4", (PyCFunction)SatrecArray_sgp4, METH_VARARGS,
     PyDoc_STR("Given array of minutes since epoch, write"
               " positions, velocities, and errors to arrays.")},
    {NULL, NULL}
};

PyAPI_DATA(PyTypeObject) SatrecArrayType = {
    PyVarObject_HEAD_INIT(NULL, sizeof(elsetrec))
    /* See the module initialization function at the bottom of this file. */
};

/* The module that ties it all together. */

static PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    m_name : "sgp4.vallado_cpp",
    m_doc : "Official C++ SGP4 implementation.",
    m_size : -1,
};

PyMODINIT_FUNC
PyInit_vallado_cpp(void)
{
    SatrecType.tp_name = "sgp4.vallado_cpp.Satrec";
    SatrecType.tp_basicsize = sizeof(SatrecObject);
    SatrecType.tp_flags = Py_TPFLAGS_DEFAULT;
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
    SatrecArrayType.tp_flags = Py_TPFLAGS_DEFAULT;
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

    return m;
}
