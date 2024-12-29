
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

#include <math.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>
#include <ufcx.h>



// Code for integral integral_27b8ebf634fc8f2fd34f37c92a7bf7b615a168b3

void tabulate_tensor_integral_27b8ebf634fc8f2fd34f37c92a7bf7b615a168b3(double* restrict A,
                                    const double* restrict w,
                                    const double* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation)
{
// Quadrature rules
static const double weights_48e[3] = {0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
// Precomputed values of basis functions and precomputations
// FE* dimensions: [permutation][entities][points][dofs]
static const double FE1_C0_D10_Q48e[1][1][1][3] = {{{{-1.0, 1.0, 0.0}}}};
static const double FE1_C1_D01_Q48e[1][1][1][3] = {{{{-1.0, 0.0, 1.0}}}};
static const double FE2_C0_Q48e[1][1][3][3] = {{{{0.6666666666666666, 0.1666666666666666, 0.1666666666666667},
  {0.1666666666666667, 0.1666666666666666, 0.6666666666666665},
  {0.1666666666666667, 0.6666666666666665, 0.1666666666666667}}}};
// ------------------------ 
// Section: Jacobian
// Inputs: coordinate_dofs, FE1_C1_D01_Q48e, FE1_C0_D10_Q48e
// Outputs: J_c3, J_c0, J_c2, J_c1
double J_c0 = 0.0;
double J_c3 = 0.0;
double J_c1 = 0.0;
double J_c2 = 0.0;
{
  for (int ic = 0; ic < 3; ++ic)
  {
    J_c0 += coordinate_dofs[(ic) * 3] * FE1_C0_D10_Q48e[0][0][0][ic];
    J_c3 += coordinate_dofs[(ic) * 3 + 1] * FE1_C1_D01_Q48e[0][0][0][ic];
    J_c1 += coordinate_dofs[(ic) * 3] * FE1_C1_D01_Q48e[0][0][0][ic];
    J_c2 += coordinate_dofs[(ic) * 3 + 1] * FE1_C0_D10_Q48e[0][0][0][ic];
  }
}
// ------------------------ 
double sp_48e_0 = J_c0 * J_c3;
double sp_48e_1 = J_c1 * J_c2;
double sp_48e_2 = -sp_48e_1;
double sp_48e_3 = sp_48e_0 + sp_48e_2;
double sp_48e_4 = fabs(sp_48e_3);
for (int iq = 0; iq < 3; ++iq)
{
  // ------------------------ 
  // Section: Intermediates
  // Inputs: 
  // Outputs: fw0
  double fw0 = 0;
  {
    fw0 = sp_48e_4 * weights_48e[iq];
  }
  // ------------------------ 
  // ------------------------ 
  // Section: Tensor Computation
  // Inputs: fw0, FE2_C0_Q48e
  // Outputs: A
  {
    double temp_0[3] = {0};
    for (int j = 0; j < 3; ++j)
    {
      temp_0[j] = fw0 * FE2_C0_Q48e[0][0][iq][j];
    }
    for (int j = 0; j < 3; ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        A[3 * (i) + (j)] += FE2_C0_Q48e[0][0][iq][i] * temp_0[j];
      }
    }
  }
  // ------------------------ 
}

}



ufcx_integral integral_27b8ebf634fc8f2fd34f37c92a7bf7b615a168b3 =
{
  .enabled_coefficients = NULL,
  .tabulate_tensor_float32 = NULL,
  .tabulate_tensor_float64 = tabulate_tensor_integral_27b8ebf634fc8f2fd34f37c92a7bf7b615a168b3,
  .tabulate_tensor_complex64 = NULL,
  .tabulate_tensor_complex128 = NULL,
  .needs_facet_permutations = false,
  .coordinate_element_hash = UINT64_C(0),
};

// End of code for integral integral_27b8ebf634fc8f2fd34f37c92a7bf7b615a168b3

// Code for integral integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc

void tabulate_tensor_integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc(double* restrict A,
                                    const double* restrict w,
                                    const double* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation)
{
// Quadrature rules
static const double weights_48e[3] = {0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
// Precomputed values of basis functions and precomputations
// FE* dimensions: [permutation][entities][points][dofs]
static const double FE0_C0_Q48e[1][1][3][3] = {{{{1.0, 0.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0}}}};
static const double FE2_C0_D10_Q48e[1][1][1][3] = {{{{-1.0, 1.0, 0.0}}}};
static const double FE2_C1_D01_Q48e[1][1][1][3] = {{{{-1.0, 0.0, 1.0}}}};
static const double FE3_C0_Q48e[1][1][3][3] = {{{{0.6666666666666666, 0.1666666666666666, 0.1666666666666667},
  {0.1666666666666667, 0.1666666666666666, 0.6666666666666665},
  {0.1666666666666667, 0.6666666666666665, 0.1666666666666667}}}};
// ------------------------ 
// Section: Jacobian
// Inputs: coordinate_dofs, FE2_C0_D10_Q48e, FE2_C1_D01_Q48e
// Outputs: J_c3, J_c0, J_c2, J_c1
double J_c0 = 0.0;
double J_c3 = 0.0;
double J_c1 = 0.0;
double J_c2 = 0.0;
{
  for (int ic = 0; ic < 3; ++ic)
  {
    J_c0 += coordinate_dofs[(ic) * 3] * FE2_C0_D10_Q48e[0][0][0][ic];
    J_c3 += coordinate_dofs[(ic) * 3 + 1] * FE2_C1_D01_Q48e[0][0][0][ic];
    J_c1 += coordinate_dofs[(ic) * 3] * FE2_C1_D01_Q48e[0][0][0][ic];
    J_c2 += coordinate_dofs[(ic) * 3 + 1] * FE2_C0_D10_Q48e[0][0][0][ic];
  }
}
// ------------------------ 
double sp_48e_0 = J_c0 * J_c3;
double sp_48e_1 = J_c1 * J_c2;
double sp_48e_2 = -sp_48e_1;
double sp_48e_3 = sp_48e_0 + sp_48e_2;
double sp_48e_4 = fabs(sp_48e_3);
for (int iq = 0; iq < 3; ++iq)
{
  // ------------------------ 
  // Section: Coefficient
  // Inputs: w, FE0_C0_Q48e
  // Outputs: w0
  double w0 = 0.0;
  {
    for (int ic = 0; ic < 3; ++ic)
    {
      w0 += w[ic] * FE0_C0_Q48e[0][0][iq][ic];
    }
  }
  // ------------------------ 
  // ------------------------ 
  // Section: Intermediates
  // Inputs: w0
  // Outputs: fw0
  double fw0 = 0;
  {
    double sv_48e_0 = sp_48e_4 * w0;
    fw0 = sv_48e_0 * weights_48e[iq];
  }
  // ------------------------ 
  // ------------------------ 
  // Section: Tensor Computation
  // Inputs: fw0, FE3_C0_Q48e
  // Outputs: A
  {
    for (int i = 0; i < 3; ++i)
    {
      A[(i)] += fw0 * FE3_C0_Q48e[0][0][iq][i];
    }
  }
  // ------------------------ 
}

}

bool enabled_coefficients_integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc[1] = {1};

ufcx_integral integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc =
{
  .enabled_coefficients = enabled_coefficients_integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc,
  .tabulate_tensor_float32 = NULL,
  .tabulate_tensor_float64 = tabulate_tensor_integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc,
  .tabulate_tensor_complex64 = NULL,
  .tabulate_tensor_complex128 = NULL,
  .needs_facet_permutations = false,
  .coordinate_element_hash = UINT64_C(0),
};

// End of code for integral integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc

// Code for form form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a


uint64_t finite_element_hashes_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a[2] = {UINT64_C(16933917890882691060), UINT64_C(16933917890882691060)};
int form_integral_offsets_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a[4] = {0, 1, 1, 1};
static ufcx_integral* form_integrals_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a[1] = {&integral_27b8ebf634fc8f2fd34f37c92a7bf7b615a168b3};
int form_integral_ids_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a[1] = {-1};




ufcx_form form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a =
{

  .signature = "7eb460d3a1efc7ffb81c6563e2406cf057469e3dcb9934a000cf50dbea1bd707f5a0ea6d4f85b93df15a05c25cf3c7e134dc9e6b6ef7a4f3f59762c7faa8853e",
  .rank = 2,
  .num_coefficients = 0,
  .num_constants = 0,
  .original_coefficient_positions = NULL,

  .coefficient_name_map = NULL,
  .constant_name_map = NULL,

  .finite_element_hashes = finite_element_hashes_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a,

  .form_integrals = form_integrals_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a,
  .form_integral_ids = form_integral_ids_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a,
  .form_integral_offsets = form_integral_offsets_form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a
};

// Alias name
ufcx_form* form_quadSca_dgSca_a = &form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a;

// End of code for form form_a8e9c11af2d97cc34a6a12ea8d9a1d7d58353f7a

// Code for form form_b600b95f942933f56dc3ab71c00249e163f88a6e

int original_coefficient_position_form_b600b95f942933f56dc3ab71c00249e163f88a6e[1] = {0};
uint64_t finite_element_hashes_form_b600b95f942933f56dc3ab71c00249e163f88a6e[2] = {UINT64_C(16933917890882691060), UINT64_C(0)};
int form_integral_offsets_form_b600b95f942933f56dc3ab71c00249e163f88a6e[4] = {0, 1, 1, 1};
static ufcx_integral* form_integrals_form_b600b95f942933f56dc3ab71c00249e163f88a6e[1] = {&integral_80ff0b6cdb052a809490a8d4aba914d05151fbfc};
int form_integral_ids_form_b600b95f942933f56dc3ab71c00249e163f88a6e[1] = {-1};

static const char* coefficient_names_form_b600b95f942933f56dc3ab71c00249e163f88a6e[1] = {"sol"};


ufcx_form form_b600b95f942933f56dc3ab71c00249e163f88a6e =
{

  .signature = "7005b7059cb321b993b0d0a6af242f0f4e59c263e0f516620b107640e2d9c015d47dac24eeac2291b609f06388db564a16731c32adcbe717b5d44950499bbfb4",
  .rank = 1,
  .num_coefficients = 1,
  .num_constants = 0,
  .original_coefficient_positions = original_coefficient_position_form_b600b95f942933f56dc3ab71c00249e163f88a6e,

  .coefficient_name_map = coefficient_names_form_b600b95f942933f56dc3ab71c00249e163f88a6e,
  .constant_name_map = NULL,

  .finite_element_hashes = finite_element_hashes_form_b600b95f942933f56dc3ab71c00249e163f88a6e,

  .form_integrals = form_integrals_form_b600b95f942933f56dc3ab71c00249e163f88a6e,
  .form_integral_ids = form_integral_ids_form_b600b95f942933f56dc3ab71c00249e163f88a6e,
  .form_integral_offsets = form_integral_offsets_form_b600b95f942933f56dc3ab71c00249e163f88a6e
};

// Alias name
ufcx_form* form_quadSca_dgSca_L = &form_b600b95f942933f56dc3ab71c00249e163f88a6e;

// End of code for form form_b600b95f942933f56dc3ab71c00249e163f88a6e
