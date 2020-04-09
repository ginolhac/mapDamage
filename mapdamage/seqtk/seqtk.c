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

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* constant table */

unsigned char seq_nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static PyObject *stk_comp(PyObject *self, PyObject *args)
{
    gzFile fp;
    kseq_t *seq;
    int l, sig;
    char *filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename))
    {
        return NULL;
    }

    fp = gzopen(filename, "r");
    if (!fp)
    {
        return PyErr_SetFromErrnoWithFilename(PyExc_OSError, filename);
    }

    PyObject *results = PyList_New(0);

    seq = kseq_init(fp);
    while (!(sig = PyErr_CheckSignals()) && (l = kseq_read(seq)) >= 0)
    {
        long counts[256];
        memset(counts, 0, 256 * sizeof(long));

        for (long i = 0; i < l; ++i)
        {
            const unsigned char na = seq->seq.s[i];

            ++counts[na];
        }

        PyObject *result = PyDict_New();
        PyDict_SetItemString(result, "name", PyUnicode_FromString(seq->name.s));
        PyDict_SetItemString(result, "len", PyLong_FromLong(l));
        PyDict_SetItemString(result, "A", PyLong_FromLong(counts['a'] + counts['A']));
        PyDict_SetItemString(result, "C", PyLong_FromLong(counts['c'] + counts['C']));
        PyDict_SetItemString(result, "G", PyLong_FromLong(counts['g'] + counts['G']));
        PyDict_SetItemString(result, "T", PyLong_FromLong(counts['t'] + counts['T']));
        PyList_Append(results, result);
    }

    kseq_destroy(seq);
    gzclose(fp);

    return sig ? NULL : results;
}

static PyMethodDef SeqtkMethods[] = {
    {"comp", stk_comp, METH_VARARGS, "Wrapper around seqtk stk_comp function"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef seqtkmodule = {
    PyModuleDef_HEAD_INIT,
    "seqtk",
    "Python interface to seqtk functions",
    -1,
    SeqtkMethods,
};

PyMODINIT_FUNC PyInit_seqtk(void)
{
    return PyModule_Create(&seqtkmodule);
}
