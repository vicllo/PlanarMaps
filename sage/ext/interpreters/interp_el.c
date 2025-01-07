/* Automatically generated by /home/victor/sage/src/sage_setup/autogen/interpreters/generator.py.  Do not edit! */
#include <Python.h>


#include "wrapper_el.h"

#define CHECK(x) do_check(&(x), domain)

static inline int do_check(PyObject **x, PyObject *domain) {
  if (*x == NULL) return 0;
  PyObject *new_x = el_check_element(*x, domain);
  Py_DECREF(*x);
  *x = new_x;
  if (*x == NULL) return 0;
  return 1;
}

PyObject* interp_el(PyObject** args,
        PyObject** constants,
        PyObject** stack,
        PyObject* domain,
        int* code) {
  while (1) {
    switch (*code++) {
    case 0: /* load_arg */
      {
        int ai0 = *code++;
        PyObject* i0 = args[ai0];
        PyObject* o0;
        o0 = i0; Py_INCREF(o0);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 1: /* load_const */
      {
        int ai0 = *code++;
        PyObject* i0 = constants[ai0];
        PyObject* o0;
        o0 = i0; Py_INCREF(o0);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 2: /* return */
      {
        PyObject* i0 = *--stack;
        *stack = NULL;
        return i0;
      }
      break;
    case 3: /* py_call */
      {
        int ai0 = *code++;
        PyObject* i0 = constants[ai0];
        int n_i1 = *code++;
        stack -= n_i1;
        PyObject** i1 = stack;
        PyObject* o0;

        PyObject *py_args = PyTuple_New(n_i1);
        if (py_args == NULL) goto error;
        int i;
        for (i = 0; i < n_i1; i++) {
          PyObject *arg = i1[i];
          PyTuple_SET_ITEM(py_args, i, arg);
          i1[i] = NULL;
        }
        o0 = PyObject_CallObject(i0, py_args);
        Py_DECREF(py_args);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 4: /* add */
      {
        PyObject* i1 = *--stack;
        *stack = NULL;
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Add(i0, i1);
        Py_DECREF(i0);
        Py_DECREF(i1);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 5: /* sub */
      {
        PyObject* i1 = *--stack;
        *stack = NULL;
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Subtract(i0, i1);
        Py_DECREF(i0);
        Py_DECREF(i1);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 6: /* mul */
      {
        PyObject* i1 = *--stack;
        *stack = NULL;
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Multiply(i0, i1);
        Py_DECREF(i0);
        Py_DECREF(i1);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 7: /* div */
      {
        PyObject* i1 = *--stack;
        *stack = NULL;
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_TrueDivide(i0, i1);
        Py_DECREF(i0);
        Py_DECREF(i1);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 8: /* floordiv */
      {
        PyObject* i1 = *--stack;
        *stack = NULL;
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_FloorDivide(i0, i1);
        Py_DECREF(i0);
        Py_DECREF(i1);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 9: /* pow */
      {
        PyObject* i1 = *--stack;
        *stack = NULL;
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Power(i0, i1, Py_None);
        Py_DECREF(i0);
        Py_DECREF(i1);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 10: /* ipow */
      {
        int ai1 = *code++;
        PyObject* i1 = constants[ai1];
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Power(i0, i1, Py_None);
        Py_DECREF(i0);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 11: /* neg */
      {
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Negative(i0);
        Py_DECREF(i0);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 12: /* invert */
      {
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Invert(i0);
        Py_DECREF(i0);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    case 13: /* abs */
      {
        PyObject* i0 = *--stack;
        *stack = NULL;
        PyObject* o0;
        o0 = PyNumber_Absolute(i0);
        Py_DECREF(i0);
        if (!CHECK(o0)) {
          Py_XDECREF(o0);
          goto error;
        }
        *stack++ = o0;
      }
      break;
    }
  }
error:
  return NULL;
}

