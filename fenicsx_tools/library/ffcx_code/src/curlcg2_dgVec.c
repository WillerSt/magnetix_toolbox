
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



// Code for integral integral_340d03dd207ec82ce66bd783f0575057141494e3

void tabulate_tensor_integral_340d03dd207ec82ce66bd783f0575057141494e3(double* restrict A,
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
static const double FE1_C0_Q48e[1][1][3][3] = {{{{0.6666666666666666, 0.1666666666666666, 0.1666666666666667},
  {0.1666666666666667, 0.1666666666666666, 0.6666666666666665},
  {0.1666666666666667, 0.6666666666666665, 0.1666666666666667}}}};
static const double FE3_C0_D10_Q48e[1][1][1][3] = {{{{-1.0, 1.0, 0.0}}}};
static const double FE3_C1_D01_Q48e[1][1][1][3] = {{{{-1.0, 0.0, 1.0}}}};
// ------------------------ 
// Section: Jacobian
// Inputs: coordinate_dofs, FE3_C1_D01_Q48e, FE3_C0_D10_Q48e
// Outputs: J_c3, J_c0, J_c2, J_c1
double J_c0 = 0.0;
double J_c3 = 0.0;
double J_c1 = 0.0;
double J_c2 = 0.0;
{
  for (int ic = 0; ic < 3; ++ic)
  {
    J_c0 += coordinate_dofs[(ic) * 3] * FE3_C0_D10_Q48e[0][0][0][ic];
    J_c3 += coordinate_dofs[(ic) * 3 + 1] * FE3_C1_D01_Q48e[0][0][0][ic];
    J_c1 += coordinate_dofs[(ic) * 3] * FE3_C1_D01_Q48e[0][0][0][ic];
    J_c2 += coordinate_dofs[(ic) * 3 + 1] * FE3_C0_D10_Q48e[0][0][0][ic];
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
  // Inputs: FE1_C0_Q48e, fw0
  // Outputs: A
  {
    double temp_0[3] = {0};
    for (int j = 0; j < 3; ++j)
    {
      temp_0[j] = fw0 * FE1_C0_Q48e[0][0][iq][j];
    }
    double temp_1[3] = {0};
    for (int j = 0; j < 3; ++j)
    {
      temp_1[j] = fw0 * FE1_C0_Q48e[0][0][iq][j];
    }
    for (int j = 0; j < 3; ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        A[6 * (2 * (i)) + 2 * (j)] += FE1_C0_Q48e[0][0][iq][i] * temp_0[j];
        A[6 * (2 * (i) + 1) + (2 * (j) + 1)] += FE1_C0_Q48e[0][0][iq][i] * temp_1[j];
      }
    }
  }
  // ------------------------ 
}

}



ufcx_integral integral_340d03dd207ec82ce66bd783f0575057141494e3 =
{
  .enabled_coefficients = NULL,
  .tabulate_tensor_float32 = NULL,
  .tabulate_tensor_float64 = tabulate_tensor_integral_340d03dd207ec82ce66bd783f0575057141494e3,
  .tabulate_tensor_complex64 = NULL,
  .tabulate_tensor_complex128 = NULL,
  .needs_facet_permutations = false,
  .coordinate_element_hash = UINT64_C(0),
};

// End of code for integral integral_340d03dd207ec82ce66bd783f0575057141494e3

// Code for integral integral_5614b9784d838f2d431c02c16ace102d5a2fcac5

void tabulate_tensor_integral_5614b9784d838f2d431c02c16ace102d5a2fcac5(double* restrict A,
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
static const double FE1_C0_Q48e[1][1][3][3] = {{{{0.6666666666666666, 0.1666666666666666, 0.1666666666666667},
  {0.1666666666666667, 0.1666666666666666, 0.6666666666666665},
  {0.1666666666666667, 0.6666666666666665, 0.1666666666666667}}}};
static const double FE2_C0_D01_Q48e[1][1][1][3] = {{{{-1.0, 0.0, 1.0}}}};
static const double FE3_C0_D10_Q48e[1][1][1][3] = {{{{-1.0, 1.0, 0.0}}}};
// ------------------------ 
// Section: Coefficient
// Inputs: FE2_C0_D01_Q48e, w, FE3_C0_D10_Q48e
// Outputs: w0_d01, w0_d10
double w0_d01 = 0.0;
double w0_d10 = 0.0;
{
  for (int ic = 0; ic < 3; ++ic)
  {
    w0_d01 += w[ic] * FE2_C0_D01_Q48e[0][0][0][ic];
    w0_d10 += w[ic] * FE3_C0_D10_Q48e[0][0][0][ic];
  }
}
// ------------------------ 
// ------------------------ 
// Section: Jacobian
// Inputs: coordinate_dofs, FE2_C0_D01_Q48e, FE3_C0_D10_Q48e
// Outputs: J_c3, J_c0, J_c2, J_c1
double J_c0 = 0.0;
double J_c3 = 0.0;
double J_c1 = 0.0;
double J_c2 = 0.0;
{
  for (int ic = 0; ic < 3; ++ic)
  {
    J_c0 += coordinate_dofs[(ic) * 3] * FE3_C0_D10_Q48e[0][0][0][ic];
    J_c3 += coordinate_dofs[(ic) * 3 + 1] * FE2_C0_D01_Q48e[0][0][0][ic];
    J_c1 += coordinate_dofs[(ic) * 3] * FE2_C0_D01_Q48e[0][0][0][ic];
    J_c2 += coordinate_dofs[(ic) * 3 + 1] * FE3_C0_D10_Q48e[0][0][0][ic];
  }
}
// ------------------------ 
double sp_48e_0 = J_c0 * J_c3;
double sp_48e_1 = J_c1 * J_c2;
double sp_48e_2 = -sp_48e_1;
double sp_48e_3 = sp_48e_0 + sp_48e_2;
double sp_48e_4 = J_c0 / sp_48e_3;
double sp_48e_5 = w0_d01 * sp_48e_4;
double sp_48e_6 = -J_c1;
double sp_48e_7 = sp_48e_6 / sp_48e_3;
double sp_48e_8 = w0_d10 * sp_48e_7;
double sp_48e_9 = sp_48e_5 + sp_48e_8;
double sp_48e_10 = J_c3 / sp_48e_3;
double sp_48e_11 = w0_d10 * sp_48e_10;
double sp_48e_12 = -J_c2;
double sp_48e_13 = sp_48e_12 / sp_48e_3;
double sp_48e_14 = w0_d01 * sp_48e_13;
double sp_48e_15 = sp_48e_11 + sp_48e_14;
double sp_48e_16 = -sp_48e_15;
double sp_48e_17 = fabs(sp_48e_3);
double sp_48e_18 = sp_48e_9 * sp_48e_17;
double sp_48e_19 = sp_48e_16 * sp_48e_17;
for (int iq = 0; iq < 3; ++iq)
{
  // ------------------------ 
  // Section: Intermediates
  // Inputs: 
  // Outputs: fw0, fw1
  double fw0 = 0;
  double fw1 = 0;
  {
    fw0 = sp_48e_18 * weights_48e[iq];
    fw1 = sp_48e_19 * weights_48e[iq];
  }
  // ------------------------ 
  // ------------------------ 
  // Section: Tensor Computation
  // Inputs: FE1_C0_Q48e, fw0, fw1
  // Outputs: A
  {
    for (int i = 0; i < 3; ++i)
    {
      A[2 * (i)] += fw0 * FE1_C0_Q48e[0][0][iq][i];
      A[(2 * (i) + 1)] += fw1 * FE1_C0_Q48e[0][0][iq][i];
    }
  }
  // ------------------------ 
}

}

bool enabled_coefficients_integral_5614b9784d838f2d431c02c16ace102d5a2fcac5[1] = {1};

ufcx_integral integral_5614b9784d838f2d431c02c16ace102d5a2fcac5 =
{
  .enabled_coefficients = enabled_coefficients_integral_5614b9784d838f2d431c02c16ace102d5a2fcac5,
  .tabulate_tensor_float32 = NULL,
  .tabulate_tensor_float64 = tabulate_tensor_integral_5614b9784d838f2d431c02c16ace102d5a2fcac5,
  .tabulate_tensor_complex64 = NULL,
  .tabulate_tensor_complex128 = NULL,
  .needs_facet_permutations = false,
  .coordinate_element_hash = UINT64_C(0),
};

// End of code for integral integral_5614b9784d838f2d431c02c16ace102d5a2fcac5

// Code for form form_64bf52c22ac0750f54afbd25b068865275d58e4a


uint64_t finite_element_hashes_form_64bf52c22ac0750f54afbd25b068865275d58e4a[2] = {UINT64_C(0), UINT64_C(0)};
int form_integral_offsets_form_64bf52c22ac0750f54afbd25b068865275d58e4a[4] = {0, 1, 1, 1};
static ufcx_integral* form_integrals_form_64bf52c22ac0750f54afbd25b068865275d58e4a[1] = {&integral_340d03dd207ec82ce66bd783f0575057141494e3};
int form_integral_ids_form_64bf52c22ac0750f54afbd25b068865275d58e4a[1] = {-1};




ufcx_form form_64bf52c22ac0750f54afbd25b068865275d58e4a =
{

  .signature = "9e95023e5dfb66f09e282134d7e976f2129319018d684c642ee8f65430a449e0a2c749c8927e9cb689405fddd8ce029b837ddf6d0fecaaa7f0c23364cb95cda4",
  .rank = 2,
  .num_coefficients = 0,
  .num_constants = 0,
  .original_coefficient_positions = NULL,

  .coefficient_name_map = NULL,
  .constant_name_map = NULL,

  .finite_element_hashes = finite_element_hashes_form_64bf52c22ac0750f54afbd25b068865275d58e4a,

  .form_integrals = form_integrals_form_64bf52c22ac0750f54afbd25b068865275d58e4a,
  .form_integral_ids = form_integral_ids_form_64bf52c22ac0750f54afbd25b068865275d58e4a,
  .form_integral_offsets = form_integral_offsets_form_64bf52c22ac0750f54afbd25b068865275d58e4a
};

// Alias name
ufcx_form* form_curlcg2_dgVec_a = &form_64bf52c22ac0750f54afbd25b068865275d58e4a;

// End of code for form form_64bf52c22ac0750f54afbd25b068865275d58e4a

// Code for form form_324460ce69424d3a150f6e5d9ce1ba67acbec627

int original_coefficient_position_form_324460ce69424d3a150f6e5d9ce1ba67acbec627[1] = {0};
uint64_t finite_element_hashes_form_324460ce69424d3a150f6e5d9ce1ba67acbec627[2] = {UINT64_C(0), UINT64_C(16933917890882727400)};
int form_integral_offsets_form_324460ce69424d3a150f6e5d9ce1ba67acbec627[4] = {0, 1, 1, 1};
static ufcx_integral* form_integrals_form_324460ce69424d3a150f6e5d9ce1ba67acbec627[1] = {&integral_5614b9784d838f2d431c02c16ace102d5a2fcac5};
int form_integral_ids_form_324460ce69424d3a150f6e5d9ce1ba67acbec627[1] = {-1};

static const char* coefficient_names_form_324460ce69424d3a150f6e5d9ce1ba67acbec627[1] = {"sol"};


ufcx_form form_324460ce69424d3a150f6e5d9ce1ba67acbec627 =
{

  .signature = "25ea82f55911a89c846ce58ca1fa176587972f543170322cc1109544c63279453158043b97f199360b9912847d59725e2c7aca7d1569ddea8dc50c785cb5574a",
  .rank = 1,
  .num_coefficients = 1,
  .num_constants = 0,
  .original_coefficient_positions = original_coefficient_position_form_324460ce69424d3a150f6e5d9ce1ba67acbec627,

  .coefficient_name_map = coefficient_names_form_324460ce69424d3a150f6e5d9ce1ba67acbec627,
  .constant_name_map = NULL,

  .finite_element_hashes = finite_element_hashes_form_324460ce69424d3a150f6e5d9ce1ba67acbec627,

  .form_integrals = form_integrals_form_324460ce69424d3a150f6e5d9ce1ba67acbec627,
  .form_integral_ids = form_integral_ids_form_324460ce69424d3a150f6e5d9ce1ba67acbec627,
  .form_integral_offsets = form_integral_offsets_form_324460ce69424d3a150f6e5d9ce1ba67acbec627
};

// Alias name
ufcx_form* form_curlcg2_dgVec_L = &form_324460ce69424d3a150f6e5d9ce1ba67acbec627;

// End of code for form form_324460ce69424d3a150f6e5d9ce1ba67acbec627
