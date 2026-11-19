/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file NL1_ww3.cpp
 * @brief Implementation of the single-point Discrete Interaction Approximation
 * (NL1) nonlinear interaction source term routine for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the NL1 nonlinear
 * wave-wave interaction source term (W3SNL1) for one spatial point, with a
 * C-compatible interface called from the WAVEWATCH III (WW3) Fortran module
 * w3snl1md.F90. The routine, the input and output spectrum (wave action
 * density N(k,theta)) and the journal references are documented in
 * NL1_ww3.h.
 *
 * @copyright © 2026 National Weather Service, National Oceanic and Atmospheric
 * Administration. WAVEWATCH IV (TM) and WW4 (TM) are trademarks of the National
 * Weather Service.
 * The National Weather Service (NWS) often uses Generative SI (GenSI, formerly
 * identified as GenAI) for code development and refactoring. Whenever GenSI is
 * used, NWS requires a full human review of code before it is added to its
 * repositories.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted and adapted from WAVEWATCH III (WW3)
 * w3snl1md.F90 (W3SNL1) by Hendrik L. Tolman.
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */

// Local
#include "NL1_ww3.h"

// System
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace ww4_nl1 {

// --- Constants --------------------------------------------------------------
/// Lower limit on the exponent of the shallow water factor (avoids underflow
/// in the exponential, value as in the WW3 Fortran code).
constexpr float kMinDepthExponent = -1.0e15F;

// --- compute_w3snl1 ---------------------------------------------------------
// Documented in NL1_ww3.h.
void compute_w3snl1(
    const std::span<const float> action_density,
    const std::span<const float> group_velocity,
    const float mean_relative_depth, const std::span<float> source_term,
    const std::span<float> source_diagonal, const std::size_t num_directions,
    const std::size_t num_extended_frequencies,
    const std::size_t num_auxiliary_bins,
    const std::span<const float> radian_frequency, const float tail_factor,
    const float dia_coefficient, const float depth_coefficient_1,
    const float depth_coefficient_2, const float depth_coefficient_3,
    const float depth_conversion_factor, const float min_relative_depth,
    const float inv_two_pi, const std::span<const int> plus_address_q1_c1,
    const std::span<const int> plus_address_q1_c2,
    const std::span<const int> plus_address_q1_c3,
    const std::span<const int> plus_address_q1_c4,
    const std::span<const int> minus_address_q1_c1,
    const std::span<const int> minus_address_q1_c2,
    const std::span<const int> minus_address_q1_c3,
    const std::span<const int> minus_address_q1_c4,
    const std::span<const int> plus_address_q2_c1,
    const std::span<const int> plus_address_q2_c2,
    const std::span<const int> plus_address_q2_c3,
    const std::span<const int> plus_address_q2_c4,
    const std::span<const int> minus_address_q2_c1,
    const std::span<const int> minus_address_q2_c2,
    const std::span<const int> minus_address_q2_c3,
    const std::span<const int> minus_address_q2_c4,
    const std::span<const int> fold_address_c1_q1,
    const std::span<const int> fold_address_c1_q2,
    const std::span<const int> fold_address_c2_q1,
    const std::span<const int> fold_address_c2_q2,
    const std::span<const int> fold_address_c3_q1,
    const std::span<const int> fold_address_c3_q2,
    const std::span<const int> fold_address_c4_q1,
    const std::span<const int> fold_address_c4_q2,
    const std::span<const int> fold_address_c5_q1,
    const std::span<const int> fold_address_c5_q2,
    const std::span<const int> fold_address_c6_q1,
    const std::span<const int> fold_address_c6_q2,
    const std::span<const int> fold_address_c7_q1,
    const std::span<const int> fold_address_c7_q2,
    const std::span<const int> fold_address_c8_q1,
    const std::span<const int> fold_address_c8_q2, const float lambda_weight_1,
    const float lambda_weight_2, const float lambda_weight_3,
    const std::span<const float> frequency_scaling, const float interp_weight_1,
    const float interp_weight_2, const float interp_weight_3,
    const float interp_weight_4, const float interp_weight_5,
    const float interp_weight_6, const float interp_weight_7,
    const float interp_weight_8, const float squared_weight_1,
    const float squared_weight_2, const float squared_weight_3,
    const float squared_weight_4, const float squared_weight_5,
    const float squared_weight_6, const float squared_weight_7,
    const float squared_weight_8) noexcept {
  const std::size_t num_frequencies = group_velocity.size();
  const std::size_t num_spectrum = action_density.size();
  const std::size_t num_interaction_bins = frequency_scaling.size();

  // The auxiliary arrays use the 1-based Fortran addresses 1 - nth to nspecy
  // (or nspecx). Address `a` is stored at element a + address_shift.
  const int address_shift = static_cast<int>(num_directions) - 1;
  const auto element = [address_shift](const int address) {
    const int shifted = address + address_shift;
    return static_cast<std::size_t>(shifted);
  };

  // === Section: Proportionality constant ==================================
  // Includes the shallow water factor
  // D = 1 + c_1 / (k d) [1 - c_2 k d] exp(-c_3 k d).
  const float relative_depth = std::max(
      depth_conversion_factor * mean_relative_depth, min_relative_depth);
  const float depth_exponent =
      std::max(kMinDepthExponent, depth_coefficient_3 * relative_depth);
  const float constant =
      dia_coefficient *
      (1.0F + (((depth_coefficient_1 / relative_depth) *
                (1.0F - (depth_coefficient_2 * relative_depth))) *
               std::exp(depth_exponent)));

  // === Section: Auxiliary spectrum and arrays =============================
  // Auxiliary spectrum F(f_r,theta), with empty bins at low frequencies.
  std::vector<float> aux_spectrum(num_directions + num_auxiliary_bins, 0.0F);
  // Conversion factor c_g / (2 pi sigma) between N(k,theta) and F(f_r,theta).
  std::vector<float> conversion_factor(num_spectrum, 0.0F);

  for (std::size_t ifr = 0; ifr < num_frequencies; ++ifr) {
    const float factor =
        (inv_two_pi / radian_frequency[ifr]) * group_velocity[ifr];
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      const std::size_t is = ith + (ifr * num_directions);
      aux_spectrum[element(static_cast<int>(is) + 1)] =
          action_density[is] / factor;
      conversion_factor[is] = factor;
    }
  }

  // Tail beyond the discrete spectrum.
  for (std::size_t ifr = num_frequencies; ifr < num_extended_frequencies;
       ++ifr) {
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      const int address = static_cast<int>(ith + (ifr * num_directions)) + 1;
      aux_spectrum[element(address)] =
          aux_spectrum[element(address - static_cast<int>(num_directions))] *
          tail_factor;
    }
  }

  // Interaction arrays, with empty bins at low frequencies.
  const std::size_t interaction_size = num_directions + num_interaction_bins;
  std::vector<float> source_quad1(interaction_size, 0.0F);
  std::vector<float> source_quad2(interaction_size, 0.0F);
  std::vector<float> diag_quad1_center(interaction_size, 0.0F);
  std::vector<float> diag_quad1_plus(interaction_size, 0.0F);
  std::vector<float> diag_quad1_minus(interaction_size, 0.0F);
  std::vector<float> diag_quad2_center(interaction_size, 0.0F);
  std::vector<float> diag_quad2_plus(interaction_size, 0.0F);
  std::vector<float> diag_quad2_minus(interaction_size, 0.0F);

  // Zero the empty bins at addresses 1 - nth to 0.
  for (int address = 1 - static_cast<int>(num_directions); address <= 0;
       ++address) {
    aux_spectrum[element(address)] = 0.0F;
    source_quad1[element(address)] = 0.0F;
    source_quad2[element(address)] = 0.0F;
    diag_quad1_center[element(address)] = 0.0F;
    diag_quad1_plus[element(address)] = 0.0F;
    diag_quad1_minus[element(address)] = 0.0F;
    diag_quad2_center[element(address)] = 0.0F;
    diag_quad2_plus[element(address)] = 0.0F;
    diag_quad2_minus[element(address)] = 0.0F;
  }

  // === Section: Interactions for the extended spectrum ====================
  for (std::size_t is = 0; is < num_interaction_bins; ++is) {
    const std::size_t address = element(static_cast<int>(is) + 1);

    // --- Energy at interacting bins
    const float e_center = aux_spectrum[address];
    const float e_plus1 =
        ((((interp_weight_1 * aux_spectrum[element(plus_address_q1_c1[is])]) +
           (interp_weight_2 * aux_spectrum[element(plus_address_q1_c2[is])])) +
          (interp_weight_3 * aux_spectrum[element(plus_address_q1_c3[is])])) +
         (interp_weight_4 * aux_spectrum[element(plus_address_q1_c4[is])]));
    const float e_minus1 =
        ((((interp_weight_5 * aux_spectrum[element(minus_address_q1_c1[is])]) +
           (interp_weight_6 * aux_spectrum[element(minus_address_q1_c2[is])])) +
          (interp_weight_7 * aux_spectrum[element(minus_address_q1_c3[is])])) +
         (interp_weight_8 * aux_spectrum[element(minus_address_q1_c4[is])]));
    const float e_plus2 =
        ((((interp_weight_1 * aux_spectrum[element(plus_address_q2_c1[is])]) +
           (interp_weight_2 * aux_spectrum[element(plus_address_q2_c2[is])])) +
          (interp_weight_3 * aux_spectrum[element(plus_address_q2_c3[is])])) +
         (interp_weight_4 * aux_spectrum[element(plus_address_q2_c4[is])]));
    const float e_minus2 =
        ((((interp_weight_5 * aux_spectrum[element(minus_address_q2_c1[is])]) +
           (interp_weight_6 * aux_spectrum[element(minus_address_q2_c2[is])])) +
          (interp_weight_7 * aux_spectrum[element(minus_address_q2_c3[is])])) +
         (interp_weight_8 * aux_spectrum[element(minus_address_q2_c4[is])]));

    // --- Contribution to interactions
    const float scaled_constant = constant * frequency_scaling[is];
    const float factor = scaled_constant * e_center;

    const float source1_a =
        e_center * ((e_plus1 * lambda_weight_1) + (e_minus1 * lambda_weight_2));
    const float source1_b =
        source1_a - ((e_plus1 * e_minus1) * lambda_weight_3);
    const float source2_a =
        e_center * ((e_plus2 * lambda_weight_1) + (e_minus2 * lambda_weight_2));
    const float source2_b =
        source2_a - ((e_plus2 * e_minus2) * lambda_weight_3);

    source_quad1[address] = factor * source1_b;
    source_quad2[address] = factor * source2_b;

    diag_quad1_center[address] = scaled_constant * (source1_a + source1_b);
    diag_quad1_plus[address] =
        factor * ((lambda_weight_1 * e_center) - (lambda_weight_3 * e_minus1));
    diag_quad1_minus[address] =
        factor * ((lambda_weight_2 * e_center) - (lambda_weight_3 * e_plus1));

    diag_quad2_center[address] = scaled_constant * (source2_a + source2_b);
    diag_quad2_plus[address] =
        factor * ((lambda_weight_1 * e_center) - (lambda_weight_3 * e_minus2));
    diag_quad2_minus[address] =
        factor * ((lambda_weight_2 * e_center) - (lambda_weight_3 * e_plus2));
  }

  // === Section: Source and diagonal terms =================================
  for (std::size_t is = 0; is < num_spectrum; ++is) {
    const std::size_t address = element(static_cast<int>(is) + 1);

    source_term[is] =
        conversion_factor[is] *
        (((((((((-2.0F * (source_quad1[address] + source_quad2[address])) +
                (interp_weight_1 *
                 (source_quad1[element(fold_address_c1_q1[is])] +
                  source_quad2[element(fold_address_c1_q2[is])]))) +
               (interp_weight_2 *
                (source_quad1[element(fold_address_c2_q1[is])] +
                 source_quad2[element(fold_address_c2_q2[is])]))) +
              (interp_weight_3 *
               (source_quad1[element(fold_address_c3_q1[is])] +
                source_quad2[element(fold_address_c3_q2[is])]))) +
             (interp_weight_4 *
              (source_quad1[element(fold_address_c4_q1[is])] +
               source_quad2[element(fold_address_c4_q2[is])]))) +
            (interp_weight_5 *
             (source_quad1[element(fold_address_c5_q1[is])] +
              source_quad2[element(fold_address_c5_q2[is])]))) +
           (interp_weight_6 *
            (source_quad1[element(fold_address_c6_q1[is])] +
             source_quad2[element(fold_address_c6_q2[is])]))) +
          (interp_weight_7 * (source_quad1[element(fold_address_c7_q1[is])] +
                              source_quad2[element(fold_address_c7_q2[is])]))) +
         (interp_weight_8 * (source_quad1[element(fold_address_c8_q1[is])] +
                             source_quad2[element(fold_address_c8_q2[is])])));

    source_diagonal[is] =
        (((((((((-2.0F *
                 (diag_quad1_center[address] + diag_quad2_center[address])) +
                (squared_weight_1 *
                 (diag_quad1_plus[element(fold_address_c1_q1[is])] +
                  diag_quad2_plus[element(fold_address_c1_q2[is])]))) +
               (squared_weight_2 *
                (diag_quad1_plus[element(fold_address_c2_q1[is])] +
                 diag_quad2_plus[element(fold_address_c2_q2[is])]))) +
              (squared_weight_3 *
               (diag_quad1_plus[element(fold_address_c3_q1[is])] +
                diag_quad2_plus[element(fold_address_c3_q2[is])]))) +
             (squared_weight_4 *
              (diag_quad1_plus[element(fold_address_c4_q1[is])] +
               diag_quad2_plus[element(fold_address_c4_q2[is])]))) +
            (squared_weight_5 *
             (diag_quad1_minus[element(fold_address_c5_q1[is])] +
              diag_quad2_minus[element(fold_address_c5_q2[is])]))) +
           (squared_weight_6 *
            (diag_quad1_minus[element(fold_address_c6_q1[is])] +
             diag_quad2_minus[element(fold_address_c6_q2[is])]))) +
          (squared_weight_7 *
           (diag_quad1_minus[element(fold_address_c7_q1[is])] +
            diag_quad2_minus[element(fold_address_c7_q2[is])]))) +
         (squared_weight_8 *
          (diag_quad1_minus[element(fold_address_c8_q1[is])] +
           diag_quad2_minus[element(fold_address_c8_q2[is])])));
  }
}

} // namespace ww4_nl1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3snl1_cpp -------------------------------------------------------------
// Documented in NL1_ww3.h.
void w3snl1_cpp(const float *a, const float *cg, float kdmean, float *s,
                float *d, int nk, int nth, int nspec, int nfrhgh, int nspecx,
                int nspecy, const float *sig, float fachfe, float snlc1,
                float snls1, float snls2, float snls3, float kdcon, float kdmn,
                float tpiinv, const int *ip11, const int *ip12, const int *ip13,
                const int *ip14, const int *im11, const int *im12,
                const int *im13, const int *im14, const int *ip21,
                const int *ip22, const int *ip23, const int *ip24,
                const int *im21, const int *im22, const int *im23,
                const int *im24, const int *ic11, const int *ic12,
                const int *ic21, const int *ic22, const int *ic31,
                const int *ic32, const int *ic41, const int *ic42,
                const int *ic51, const int *ic52, const int *ic61,
                const int *ic62, const int *ic71, const int *ic72,
                const int *ic81, const int *ic82, float dal1, float dal2,
                float dal3, const float *af11, float awg1, float awg2,
                float awg3, float awg4, float awg5, float awg6, float awg7,
                float awg8, float swg1, float swg2, float swg3, float swg4,
                float swg5, float swg6, float swg7, float swg8) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  const auto nspecz = static_cast<std::size_t>(nspec);
  const auto nfrhghz = static_cast<std::size_t>(nfrhgh);
  const auto nspecxz = static_cast<std::size_t>(nspecx);
  const auto nspecyz = static_cast<std::size_t>(nspecy);
  ww4_nl1::compute_w3snl1(
      std::span<const float>(a, nspecz), std::span<const float>(cg, nkz),
      kdmean, std::span<float>(s, nspecz), std::span<float>(d, nspecz), nthz,
      nfrhghz, nspecyz, std::span<const float>(sig, nkz), fachfe, snlc1, snls1,
      snls2, snls3, kdcon, kdmn, tpiinv, std::span<const int>(ip11, nspecxz),
      std::span<const int>(ip12, nspecxz), std::span<const int>(ip13, nspecxz),
      std::span<const int>(ip14, nspecxz), std::span<const int>(im11, nspecxz),
      std::span<const int>(im12, nspecxz), std::span<const int>(im13, nspecxz),
      std::span<const int>(im14, nspecxz), std::span<const int>(ip21, nspecxz),
      std::span<const int>(ip22, nspecxz), std::span<const int>(ip23, nspecxz),
      std::span<const int>(ip24, nspecxz), std::span<const int>(im21, nspecxz),
      std::span<const int>(im22, nspecxz), std::span<const int>(im23, nspecxz),
      std::span<const int>(im24, nspecxz), std::span<const int>(ic11, nspecz),
      std::span<const int>(ic12, nspecz), std::span<const int>(ic21, nspecz),
      std::span<const int>(ic22, nspecz), std::span<const int>(ic31, nspecz),
      std::span<const int>(ic32, nspecz), std::span<const int>(ic41, nspecz),
      std::span<const int>(ic42, nspecz), std::span<const int>(ic51, nspecz),
      std::span<const int>(ic52, nspecz), std::span<const int>(ic61, nspecz),
      std::span<const int>(ic62, nspecz), std::span<const int>(ic71, nspecz),
      std::span<const int>(ic72, nspecz), std::span<const int>(ic81, nspecz),
      std::span<const int>(ic82, nspecz), dal1, dal2, dal3,
      std::span<const float>(af11, nspecxz), awg1, awg2, awg3, awg4, awg5, awg6,
      awg7, awg8, swg1, swg2, swg3, swg4, swg5, swg6, swg7, swg8);
}

} // extern "C"
