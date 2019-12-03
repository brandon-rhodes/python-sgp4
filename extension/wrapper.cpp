#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "SGP4.h"
#include "structmember.h"

typedef struct {
    PyObject_HEAD
    elsetrec satrec;
} SatrecObject;

static PyObject *
Satrec_twoline2rv(PyObject *self, PyObject *args)
{
    char *string1, *string2, line1[130], line2[130];
    double dummy;

    if (!PyArg_ParseTuple(args, "ss:twoline2rv", &string1, &string2))
        return NULL;

    // Make copies, since twoline2rv() writes to both strings.
    strncpy(line1, string1, 130);
    strncpy(line2, string2, 130);
    line1[129] = '\0';
    line2[129] = '\0';

    SGP4Funcs::twoline2rv(line1, line2, ' ', ' ', 'i', wgs84,
                          dummy, dummy, dummy,
                          ((SatrecObject*)self)->satrec);
    Py_RETURN_NONE;
}

static PyObject *
Satrec_sgp4(PyObject *self, PyObject *args)
{
    double tsince, r[3], v[3];
    if (!PyArg_ParseTuple(args, "d:sgp4", &tsince))
        return NULL;
    SGP4Funcs::sgp4(((SatrecObject*)self)->satrec, tsince, r, v);
    return Py_BuildValue("(fff)(fff)", r[0], r[1], r[2], v[0], v[1], v[2]);
}

static PyMethodDef Satrec_methods[] = {
    {"twoline2rv", (PyCFunction)Satrec_twoline2rv, METH_VARARGS,
     PyDoc_STR("Initialize the record from two lines of TLE text.")},
    {"sgp4", (PyCFunction)Satrec_sgp4, METH_VARARGS,
     PyDoc_STR("Given minutes since epoch, return a position and velocity.")},
    {NULL, NULL}
};

static PyMemberDef Satrec_members[] = {
    {"satnum", T_INT, offsetof(SatrecObject, satrec.satnum), READONLY,
     PyDoc_STR("Satellite number (characters 3-7 of each TLE line).")},
    {"method", T_CHAR, offsetof(SatrecObject, satrec.method), READONLY,
     PyDoc_STR("Method 'n' near earth or 'd' deep space.")},
    {NULL}
};

static PyTypeObject SatrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    tp_name : "sgp4.vallado_cpp.Satrec",
    tp_basicsize : sizeof(SatrecObject),
    tp_itemsize : 0,
    tp_flags : Py_TPFLAGS_DEFAULT,
    tp_doc : "SGP4 satellite record",
    tp_methods : Satrec_methods,
    tp_members : Satrec_members,
    tp_new : PyType_GenericNew,
};

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

    m = PyModule_Create(&module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&SatrecType);
    if (PyModule_AddObject(m, "Satrec", (PyObject *) &SatrecType) < 0) {
        Py_DECREF(&SatrecType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
