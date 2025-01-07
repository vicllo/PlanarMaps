# Automatically generated by /home/victor/sage/src/sage_setup/autogen/interpreters/generator.py.  Do not edit!

from cpython.ref cimport PyObject

from sage.ext.fast_callable cimport Wrapper

# This is to work around a header incompatibility with PARI using
# "I" as variable conflicting with the complex "I".
# If we cimport pari earlier, we avoid this problem.
cimport cypari2.types

# We need the type double_complex to work around
#   http://trac.cython.org/ticket/869
# so this is a bit hackish.
cdef extern from "complex.h":
    ctypedef double double_complex "double complex"

cdef class Wrapper_cdf(Wrapper):
    cdef int _n_args
    cdef double_complex* _args
    cdef int _n_constants
    cdef double_complex* _constants
    cdef object _list_py_constants
    cdef int _n_py_constants
    cdef PyObject** _py_constants
    cdef int _n_stack
    cdef double_complex* _stack
    cdef int _n_code
    cdef int* _code
    cdef bint call_c(self,
                     double_complex* args,
                     double_complex* result) except 0
