/* Automatically generated by /home/victor/sage/src/sage_setup/autogen/interpreters/generator.py.  Do not edit! */
#include <Python.h>


#include <mpc.h>
#include "wrapper_cc.h"

int interp_cc(mpc_t* args,
        mpc_t retval,
        mpc_t* constants,
        PyObject** py_constants,
        mpc_t* stack,
        int* code,
        PyObject* domain) {
  while (1) {
    switch (*code++) {
    case 0: /* load_arg */
      {
        int ai0 = *code++;
        mpc_ptr i0 = args[ai0];
        mpc_ptr o0 = *stack++;
        mpc_set(o0, i0, MPC_RNDNN);
      }
      break;
    case 1: /* load_const */
      {
        int ai0 = *code++;
        mpc_ptr i0 = constants[ai0];
        mpc_ptr o0 = *stack++;
        mpc_set(o0, i0, MPC_RNDNN);
      }
      break;
    case 2: /* return */
      {
        mpc_ptr i0 = *--stack;
        mpc_set(retval, i0, MPC_RNDNN);
        return 1;
      }
      break;
    case 3: /* py_call */
      {
        int ai0 = *code++;
        PyObject* i0 = py_constants[ai0];
        int n_i1 = *code++;
        stack -= n_i1;
        mpc_t* i1 = stack;
        mpc_ptr o0 = *stack++;

          if (!cc_py_call_helper(domain, i0, n_i1, i1, o0)) {
          goto error;
        }
      }
      break;
    case 4: /* add */
      {
        mpc_ptr i1 = *--stack;
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_add(o0, i0, i1, MPC_RNDNN);
      }
      break;
    case 5: /* sub */
      {
        mpc_ptr i1 = *--stack;
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_sub(o0, i0, i1, MPC_RNDNN);
      }
      break;
    case 6: /* mul */
      {
        mpc_ptr i1 = *--stack;
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_mul(o0, i0, i1, MPC_RNDNN);
      }
      break;
    case 7: /* div */
      {
        mpc_ptr i1 = *--stack;
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_div(o0, i0, i1, MPC_RNDNN);
      }
      break;
    case 8: /* pow */
      {
        mpc_ptr i1 = *--stack;
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_pow(o0, i0, i1, MPC_RNDNN);
      }
      break;
    case 9: /* ipow */
      {
        int i1 = *code++;
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_pow_si(o0, i0, i1, MPC_RNDNN);
      }
      break;
    case 10: /* neg */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_neg(o0, i0, MPC_RNDNN);
      }
      break;
    case 11: /* log */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_log(o0, i0, MPC_RNDNN);
      }
      break;
    case 12: /* log10 */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_log10(o0, i0, MPC_RNDNN);
      }
      break;
    case 13: /* exp */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_exp(o0, i0, MPC_RNDNN);
      }
      break;
    case 14: /* cos */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_cos(o0, i0, MPC_RNDNN);
      }
      break;
    case 15: /* sin */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_sin(o0, i0, MPC_RNDNN);
      }
      break;
    case 16: /* tan */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_tan(o0, i0, MPC_RNDNN);
      }
      break;
    case 17: /* acos */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_acos(o0, i0, MPC_RNDNN);
      }
      break;
    case 18: /* asin */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_asin(o0, i0, MPC_RNDNN);
      }
      break;
    case 19: /* atan */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_atan(o0, i0, MPC_RNDNN);
      }
      break;
    case 20: /* cosh */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_cosh(o0, i0, MPC_RNDNN);
      }
      break;
    case 21: /* sinh */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_sinh(o0, i0, MPC_RNDNN);
      }
      break;
    case 22: /* tanh */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_tanh(o0, i0, MPC_RNDNN);
      }
      break;
    case 23: /* acosh */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_acosh(o0, i0, MPC_RNDNN);
      }
      break;
    case 24: /* asinh */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_asinh(o0, i0, MPC_RNDNN);
      }
      break;
    case 25: /* atanh */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_atanh(o0, i0, MPC_RNDNN);
      }
      break;
    case 26: /* invert */
      {
        mpc_ptr i0 = *--stack;
        mpc_ptr o0 = *stack++;
        mpc_ui_div(o0, 1, i0, MPC_RNDNN);
      }
      break;
    }
  }
error:
  return 0;
}

