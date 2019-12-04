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

    SGP4Funcs::twoline2rv(line1, line2, ' ', ' ', 'i', wgs84,
                          dummy, dummy, dummy, self->satrec);

    return (PyObject*) self;
}

static PyObject *
Satrec_sgp4(PyObject *self, PyObject *args)
{
    double tsince, r[3], v[3];
    if (!PyArg_ParseTuple(args, "d:sgp4", &tsince))
        return NULL;
    SGP4Funcs::sgp4(((SatrecObject *)self)->satrec, tsince, r, v);
    return Py_BuildValue("(fff)(fff)i", r[0], r[1], r[2], v[0], v[1], v[2],
                         ((SatrecObject*)self)->satrec.error);
}

static PyMethodDef Satrec_methods[] = {
    {"twoline2rv", (PyCFunction)Satrec_twoline2rv, METH_VARARGS | METH_CLASS,
     PyDoc_STR("Initialize the record from two lines of TLE text.")},
    {"sgp4", (PyCFunction)Satrec_sgp4, METH_VARARGS,
     PyDoc_STR("Given minutes since epoch, return a position and velocity.")},
    {NULL, NULL}
};

static PyMemberDef Satrec_members[] = {
    {"satnum", T_INT, offsetof(SatrecObject, satrec.satnum), READONLY,
     PyDoc_STR("Satellite number (characters 3-7 of each TLE line).")},
    {"method", T_CHAR, offsetof(SatrecObject, satrec.method), READONLY,
     PyDoc_STR("Method, either 'n' near earth or 'd' deep space.")},
    /* TODO: expose other elements that are loaded from the TLE */
    {NULL}
};

static PyTypeObject SatrecType = {
    PyObject_HEAD_INIT(NULL)
    tp_name : "sgp4.vallado_cpp.Satrec",
    tp_basicsize : sizeof(SatrecObject),
    tp_flags : Py_TPFLAGS_DEFAULT,
    tp_doc : "SGP4 satellite record.",
    tp_methods : Satrec_methods,
    tp_members : Satrec_members,
    tp_new : PyType_GenericNew,
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
    PyObject *tsince_arg, *r_arg, *v_arg, *e_arg;
    Py_buffer tsince, r, v, e;
    PyObject *rv = NULL;

    tsince.buf = r.buf = v.buf = e.buf = NULL;  // see "cleanup:" below

    if (!PyArg_ParseTuple(args, "OOOO:sgp4",
                          &tsince_arg, &r_arg, &v_arg, &e_arg))
        return NULL;

    if (PyObject_GetBuffer(tsince_arg, &tsince, PyBUF_SIMPLE)) goto cleanup;
    if (PyObject_GetBuffer(r_arg, &r, PyBUF_WRITABLE)) goto cleanup;
    if (PyObject_GetBuffer(v_arg, &v, PyBUF_WRITABLE)) goto cleanup;
    if (PyObject_GetBuffer(e_arg, &e, PyBUF_WRITABLE)) goto cleanup;

    // This extra block allows "goto" to cross declarations.
    {
        SatrecArrayObject *satrec_array = (SatrecArrayObject*) self;
        Py_ssize_t imax = ((SatrecArrayObject*) self)->ob_base.ob_size;
        Py_ssize_t jmax = tsince.len / sizeof(double);

        double *tp = (double*) tsince.buf;
        double *rp = (double*) r.buf;
        double *vp = (double*) v.buf;
        uint8_t *ep = (uint8_t*) e.buf;

        if ((r.len != (Py_ssize_t) sizeof(double) * imax * jmax * 3) ||
            (v.len != (Py_ssize_t) sizeof(double) * imax * jmax * 3) ||
            (e.len != (Py_ssize_t) sizeof(uint8_t) * imax * jmax)) {
            PyErr_SetString(PyExc_ValueError, "bad output array dimension");
            goto cleanup;
        }

        for (Py_ssize_t i=0; i < imax; i++) {
            for (Py_ssize_t j=0; j < jmax; j++) {
                SGP4Funcs::sgp4(satrec_array->satrec[i], tp[j], rp, vp);
                rp += 3;
                vp += 3;
                *(ep++) = (uint8_t) satrec_array->satrec[i].error;
            }
        }
    }

    Py_INCREF(Py_None);
    rv = Py_None;
 cleanup:
    if (tsince.buf) PyBuffer_Release(&tsince);
    if (r.buf) PyBuffer_Release(&r);
    if (v.buf) PyBuffer_Release(&v);
    if (e.buf) PyBuffer_Release(&e);
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
    tp_name : "sgp4.vallado_cpp.SatrecArray",
    tp_basicsize : sizeof(SatrecArrayObject),
    tp_itemsize : sizeof(elsetrec),
    tp_as_sequence : &SatrecArray_as_sequence,
    tp_flags : Py_TPFLAGS_DEFAULT,
    tp_doc : "SGP4 array of satellites.",
    tp_methods : SatrecArray_methods,
    tp_init : (initproc) SatrecArray_init,
    tp_new : SatrecArray_new,
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
    PyObject *m;
    if (PyType_Ready(&SatrecType) < 0)
        return NULL;

    if (PyType_Ready(&SatrecArrayType) < 0)
        return NULL;

    m = PyModule_Create(&module);
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
