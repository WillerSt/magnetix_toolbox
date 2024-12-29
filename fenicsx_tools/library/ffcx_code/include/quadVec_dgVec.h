
// This code conforms with the UFC specification version 2018.2.0.dev0
// and was automatically generated by FFCx version 0.9.0.
//
// This code was generated with the following options:
//
//  {'epsilon': 1e-14,
//   'output_directory': '.',
//   'profile': False,
//   'scalar_type': 'float64',
//   'sum_factorization': False,
//   'table_atol': 1e-09,
//   'table_rtol': 1e-06,
//   'ufl_file': ['A_form_NR.py',
//                'curl_evaluation.py',
//                'curl_projection.py',
//                'curlcg2_dgVec.py',
//                'magForce_virtualWork.py',
//                'magForce_virtualWork_bak.py',
//                'quadSca_dgSca.py',
//                'quadVec_dgVec.py',
//                'quadVec_dgVec_rotated.py'],
//   'verbosity': 30,
//   'visualise': False}

#pragma once
#include <ufcx.h>

#ifdef __cplusplus
extern "C" {
#endif

extern ufcx_integral integral_d11e9f13eb903294512ed3d00b424a0714fce208;

extern ufcx_integral integral_bdbfd4159291cd8be0a94f63631b7b6758205638;

extern ufcx_form form_4c1705aee59c677a88bb5dd5e4a2ccd4a48de863;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_quadVec_dgVec_a;


extern ufcx_form form_fa45f608b560b8e93b1a7e7a1525f41585c26dd3;

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* form_quadVec_dgVec_L;


#ifdef __cplusplus
}
#endif
