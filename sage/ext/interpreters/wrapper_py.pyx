# Automatically generated by /home/victor/sage/src/sage_setup/autogen/interpreters/generator.py.  Do not edit!
from cpython.number cimport PyNumber_TrueDivide


from cpython.ref cimport PyObject
cdef extern from "Python.h":
    void Py_DECREF(PyObject *o)
    void Py_INCREF(PyObject *o)
    void Py_CLEAR(PyObject *o)

    object PyList_New(Py_ssize_t len)
    ctypedef struct PyListObject:
        PyObject **ob_item

    ctypedef struct PyTupleObject:
        PyObject **ob_item

from cysignals.memory cimport check_allocarray, sig_free

from sage.ext.fast_callable cimport Wrapper

cdef extern from "interp_py.c":
    object interp_py(PyObject** args,
        PyObject** constants,
        PyObject** stack,
        int* code)
cdef class Wrapper_py(Wrapper):
    # attributes are declared in corresponding .pxd file

    def __init__(self, args):
        Wrapper.__init__(self, args, metadata)
        cdef int i
        cdef int count
        count = args['args']
        self._n_args = count
        val = args['constants']
        self._n_constants = len(val)
        self._list_constants = PyList_New(self._n_constants)
        self._constants = (<PyListObject *>self._list_constants).ob_item
        for i in range(len(val)):
            self._constants[i] = <PyObject *>val[i]; Py_INCREF(self._constants[i])
        count = args['stack']
        self._n_stack = count
        self._list_stack = PyList_New(self._n_stack)
        self._stack = (<PyListObject *>self._list_stack).ob_item
        val = args['code']
        self._n_code = len(val)
        self._code = <int*>check_allocarray(self._n_code, sizeof(int))
        for i in range(len(val)):
            self._code[i] = val[i]

    def __dealloc__(self):
        cdef int i
        if self._code:
            sig_free(self._code)

    def __call__(self, *args):
        if self._n_args != len(args): raise ValueError
        try:
            return interp_py((<PyTupleObject*>args).ob_item
                , self._constants
                , self._stack
                , self._code
                )
        except BaseException:
            for i in range(self._n_stack):
                Py_CLEAR(self._stack[i])
            raise


from sage.ext.fast_callable import CompilerInstrSpec, InterpreterMetadata
metadata = InterpreterMetadata(by_opname={
  'load_arg':
  (CompilerInstrSpec(0, 1, ['args']), 0),
  'load_const':
  (CompilerInstrSpec(0, 1, ['constants']), 1),
  'return':
  (CompilerInstrSpec(1, 0, []), 2),
  'py_call':
  (CompilerInstrSpec(0, 1, ['constants', 'n_inputs']), 3),
  'add':
  (CompilerInstrSpec(2, 1, []), 4),
  'sub':
  (CompilerInstrSpec(2, 1, []), 5),
  'mul':
  (CompilerInstrSpec(2, 1, []), 6),
  'div':
  (CompilerInstrSpec(2, 1, []), 7),
  'floordiv':
  (CompilerInstrSpec(2, 1, []), 8),
  'pow':
  (CompilerInstrSpec(2, 1, []), 9),
  'ipow':
  (CompilerInstrSpec(1, 1, ['constants']), 10),
  'neg':
  (CompilerInstrSpec(1, 1, []), 11),
  'invert':
  (CompilerInstrSpec(1, 1, []), 12),
  'abs':
  (CompilerInstrSpec(1, 1, []), 13),
 },
 by_opcode=[
  ('load_arg',
   CompilerInstrSpec(0, 1, ['args'])),
  ('load_const',
   CompilerInstrSpec(0, 1, ['constants'])),
  ('return',
   CompilerInstrSpec(1, 0, [])),
  ('py_call',
   CompilerInstrSpec(0, 1, ['constants', 'n_inputs'])),
  ('add',
   CompilerInstrSpec(2, 1, [])),
  ('sub',
   CompilerInstrSpec(2, 1, [])),
  ('mul',
   CompilerInstrSpec(2, 1, [])),
  ('div',
   CompilerInstrSpec(2, 1, [])),
  ('floordiv',
   CompilerInstrSpec(2, 1, [])),
  ('pow',
   CompilerInstrSpec(2, 1, [])),
  ('ipow',
   CompilerInstrSpec(1, 1, ['constants'])),
  ('neg',
   CompilerInstrSpec(1, 1, [])),
  ('invert',
   CompilerInstrSpec(1, 1, [])),
  ('abs',
   CompilerInstrSpec(1, 1, [])),
 ],
 ipow_range=True)
