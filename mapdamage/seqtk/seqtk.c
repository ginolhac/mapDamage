/* The MIT License

   Copyright (c) 20082-2012 by Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int
mapdamage_set(PyObject *dict, const char *key, PyObject *value)
{
    if (value == NULL) {
        return 1;
    }

    int result = PyDict_SetItemString(dict, key, value);
    Py_DECREF(value);

    return result;
}

int
mapdamage_set_long(PyObject *dict, const char *key, long value)
{
    return mapdamage_set(dict, key, PyLong_FromLong(value));
}

static PyObject *
mapdamage_stk_comp(PyObject *self, PyObject *args)
{
    (void)self;

    char *filename = NULL;
    if (PyArg_ParseTuple(args, "s", &filename) == 0) {
        return NULL;
    }

    gzFile fp = gzopen(filename, "r");
    if (fp == NULL) {
        PyErr_SetFromErrnoWithFilename(PyExc_OSError, filename);
        return NULL;
    }

    // Cannot fail; kseq_init will segfault if calloc fails
    kseq_t *seq = kseq_init(fp);

    long int l = 0;
    int sig = 1;
    PyObject *result = PyList_New(0);
    if (!result) {
        goto Cleanup;
    }

    long counts[256];
    while (!(sig = PyErr_CheckSignals()) && (l = kseq_read(seq)) >= 0) {
        memset(counts, 0, 256 * sizeof(long));
        for (long i = 0; i < l; ++i) {
            counts[(unsigned char)seq->seq.s[i]]++;
        }

        // If we break it will have been due to an Py* failure
        sig = 1;

        PyObject *stats = PyDict_New();
        if (stats == NULL || PyList_Append(result, stats) != 0) {
            Py_XDECREF(stats);
            break;
        }

        PyObject *name = PyUnicode_FromString(seq->name.s);
        if (mapdamage_set(stats, "name", name) != 0 ||
            mapdamage_set_long(stats, "len", l) != 0 ||
            mapdamage_set_long(stats, "A", counts['A'] + counts['a']) != 0 ||
            mapdamage_set_long(stats, "C", counts['C'] + counts['c']) != 0 ||
            mapdamage_set_long(stats, "G", counts['G'] + counts['g']) != 0 ||
            mapdamage_set_long(stats, "T", counts['T'] + counts['t']) != 0) {
            Py_DECREF(stats);
            break;
        }
    }

Cleanup:
    kseq_destroy(seq);
    const int gzerr = gzclose(fp);

    if (gzerr != Z_OK || sig || l != -1) {
        Py_XDECREF(result);

        if (gzerr != Z_OK) {
            switch (gzerr) {
                case Z_STREAM_ERROR:
                    PyErr_SetString(PyExc_RuntimeError,
                                    "Stream error while reading gzip file");
                    break;
                case Z_ERRNO:
                    PyErr_SetFromErrnoWithFilename(PyExc_OSError, filename);
                    break;
                case Z_BUF_ERROR:
                    PyErr_SetString(PyExc_RuntimeError,
                                    "Buffer error while reading gzip file");
                    break;
                default:
                    PyErr_SetString(PyExc_RuntimeError,
                                    "Unexpected zlib error");
            }
        }
        else if (l == -3) {
            PyErr_SetFromErrnoWithFilename(PyExc_OSError, filename);
        }
        else if (l == -2) {
            PyErr_SetString(PyExc_RuntimeError, "Malformed file");
        }

        return NULL;
    }

    return result;
}

static PyMethodDef SeqtkMethods[] = {
    {"comp", mapdamage_stk_comp, METH_VARARGS,
     "Wrapper around seqtk stk_comp function"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef seqtkmodule = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "seqtk",
    .m_doc = "Python interface to seqtk functions",
    .m_size = 0,
    .m_methods = SeqtkMethods,
    .m_slots = NULL,
    .m_traverse = NULL,
    .m_clear = NULL,
    .m_free = NULL,
};

PyMODINIT_FUNC
PyInit_seqtk(void)
{
    return PyModule_Create(&seqtkmodule);
}
