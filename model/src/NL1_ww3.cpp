/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file NL1_ww3.cpp
 * @brief Implementation of the single-point Discrete Interaction Approximation
 * (NL1) nonlinear interaction routines for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the NL1 nonlinear wave-wave
 * interaction source term (W3SNL1) and its initialization (INSNL1) for one
 * spatial point, with C-compatible interfaces called from the WAVEWATCH III
 * (WW3) Fortran module w3snl1md.F90. The routines, the input and output
 * spectrum (wave action density N(k,theta)) and the journal references are
 * documented in NL1_ww3.h.
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
 * w3snl1md.F90 (W3SNL1, INSNL1) by Hendrik L. Tolman.
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

namespace {

/// Square of x, as the expansion of x**2 by gfortran.
constexpr float square(const float value) noexcept { return value * value; }

/// Fourth power of x, as the expansion of x**4 by gfortran (x**2 squared).
constexpr float fourth_power(const float value) noexcept {
  return square(square(value));
}

/// Eleventh power of x, with the multiplication chain gfortran uses for the
/// constant exponent in x**11: x**11 = x**6 * x**5, x**6 = x**3 * x**3,
/// x**5 = x**3 * x**2, x**3 = x**2 * x. The chain defines the rounding, so it
/// must not be changed.
constexpr float eleventh_power(const float value) noexcept {
  const float power2 = value * value;
  const float power3 = power2 * value;
  const float power5 = power3 * power2;
  const float power6 = power3 * power3;
  return power6 * power5;
}

/// Integer power of x for a run-time exponent, with the square-and-multiply
/// algorithm of libgcc (__powisf2) that gfortran uses for x**n with variable
/// n. The algorithm defines the rounding, so it must not be changed.
float integer_power(float value, const int exponent) noexcept {
  unsigned int remaining = exponent < 0
                               ? 0U - static_cast<unsigned int>(exponent)
                               : static_cast<unsigned int>(exponent);
  float result = ((remaining % 2U) != 0U) ? value : 1.0F;
  while ((remaining >>= 1U) != 0U) {
    value = value * value;
    if ((remaining % 2U) != 0U) {
      result = result * value;
    }
  }
  return exponent < 0 ? 1.0F / result : result;
}

/// Frequency indices and weights of the quadruplet components at
/// (1 + lambda) f_r and (1 - lambda) f_r (section 4 of INSNL1).
struct FrequencyGeometry {
  int plus;          ///< IFRP, index of the lower neighbor at (1 + lambda) f_r.
  int plus1;         ///< IFRP1 = IFRP + 1.
  int minus;         ///< IFRM, index of the upper neighbor at (1 - lambda) f_r.
  int minus1;        ///< IFRM1 = IFRM - 1.
  float plus_weight; ///< WFRP.
  float plus1_weight;  ///< WFRP1.
  float minus_weight;  ///< WFRM.
  float minus1_weight; ///< WFRM1.
};

/// Computes the frequency geometry for the relative frequency ratio and
/// lambda (section 4 of INSNL1).
FrequencyGeometry
compute_frequency_geometry(const float frequency_ratio,
                           const float quadruplet_lambda) noexcept {
  FrequencyGeometry geometry{};
  const float log_ratio = std::log(frequency_ratio);

  geometry.plus =
      static_cast<int>(std::log(1.0F + quadruplet_lambda) / log_ratio);
  geometry.plus1 = geometry.plus + 1;
  geometry.plus_weight = ((1.0F + quadruplet_lambda) -
                          integer_power(frequency_ratio, geometry.plus)) /
                         (integer_power(frequency_ratio, geometry.plus1) -
                          integer_power(frequency_ratio, geometry.plus));
  geometry.plus1_weight = 1.0F - geometry.plus_weight;

  geometry.minus =
      static_cast<int>(std::log(1.0F - quadruplet_lambda) / log_ratio);
  geometry.minus1 = geometry.minus - 1;
  geometry.minus_weight = (integer_power(frequency_ratio, geometry.minus) -
                           (1.0F - quadruplet_lambda)) /
                          (integer_power(frequency_ratio, geometry.minus) -
                           integer_power(frequency_ratio, geometry.minus1));
  geometry.minus1_weight = 1.0F - geometry.minus_weight;
  return geometry;
}

} // namespace

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

// --- compute_insnl1_ranges --------------------------------------------------
// Documented in NL1_ww3.h.
void compute_insnl1_ranges(const std::size_t num_frequencies,
                           const std::size_t num_directions,
                           const float frequency_ratio,
                           const float quadruplet_lambda,
                           std::size_t &num_extended_frequencies,
                           std::size_t &num_interaction_frequencies,
                           std::size_t &num_auxiliary_bins,
                           std::size_t &num_interaction_bins) noexcept {
  // === Section: Frequency indices =========================================
  const FrequencyGeometry geometry =
      compute_frequency_geometry(frequency_ratio, quadruplet_lambda);

  // === Section: Range of calculations =====================================
  const auto num_frequencies_int = static_cast<int>(num_frequencies);
  num_extended_frequencies = static_cast<std::size_t>(
      (num_frequencies_int + geometry.plus1) - geometry.minus1);
  num_interaction_frequencies =
      static_cast<std::size_t>(num_frequencies_int - geometry.minus1);
  num_auxiliary_bins = num_extended_frequencies * num_directions;
  num_interaction_bins = num_interaction_frequencies * num_directions;
}

// --- compute_insnl1 ---------------------------------------------------------
// Documented in NL1_ww3.h.
void compute_insnl1(
    const std::span<const float> radian_frequency,
    const std::size_t num_directions, const float direction_increment,
    const float frequency_ratio, const float quadruplet_lambda,
    const float inv_two_pi, float &lambda_weight_1, float &lambda_weight_2,
    float &lambda_weight_3, const std::span<int> plus_address_q1_c1,
    const std::span<int> plus_address_q1_c2,
    const std::span<int> plus_address_q1_c3,
    const std::span<int> plus_address_q1_c4,
    const std::span<int> minus_address_q1_c1,
    const std::span<int> minus_address_q1_c2,
    const std::span<int> minus_address_q1_c3,
    const std::span<int> minus_address_q1_c4,
    const std::span<int> plus_address_q2_c1,
    const std::span<int> plus_address_q2_c2,
    const std::span<int> plus_address_q2_c3,
    const std::span<int> plus_address_q2_c4,
    const std::span<int> minus_address_q2_c1,
    const std::span<int> minus_address_q2_c2,
    const std::span<int> minus_address_q2_c3,
    const std::span<int> minus_address_q2_c4,
    const std::span<int> fold_address_c1_q1,
    const std::span<int> fold_address_c1_q2,
    const std::span<int> fold_address_c2_q1,
    const std::span<int> fold_address_c2_q2,
    const std::span<int> fold_address_c3_q1,
    const std::span<int> fold_address_c3_q2,
    const std::span<int> fold_address_c4_q1,
    const std::span<int> fold_address_c4_q2,
    const std::span<int> fold_address_c5_q1,
    const std::span<int> fold_address_c5_q2,
    const std::span<int> fold_address_c6_q1,
    const std::span<int> fold_address_c6_q2,
    const std::span<int> fold_address_c7_q1,
    const std::span<int> fold_address_c7_q2,
    const std::span<int> fold_address_c8_q1,
    const std::span<int> fold_address_c8_q2,
    const std::span<float> frequency_scaling, float &interp_weight_1,
    float &interp_weight_2, float &interp_weight_3, float &interp_weight_4,
    float &interp_weight_5, float &interp_weight_6, float &interp_weight_7,
    float &interp_weight_8, float &squared_weight_1, float &squared_weight_2,
    float &squared_weight_3, float &squared_weight_4, float &squared_weight_5,
    float &squared_weight_6, float &squared_weight_7,
    float &squared_weight_8) noexcept {
  const std::size_t num_frequencies = radian_frequency.size();
  const std::size_t num_spectrum = num_directions * num_frequencies;
  const std::size_t num_interaction_bins = frequency_scaling.size();
  const auto nth = static_cast<int>(num_directions);
  const auto nfr = static_cast<int>(num_frequencies);

  // === Section: Internal angles of the quadruplet =========================
  const float lamm2 = square(1.0F - quadruplet_lambda);
  const float lamp2 = square(1.0F + quadruplet_lambda);
  const float delta_theta3 =
      std::acos((((square(lamm2) + 4.0F) - square(lamp2)) / (4.0F * lamm2)));
  const float delta_theta4 =
      std::asin(-((std::sin(delta_theta3) * lamm2) / lamp2));

  // === Section: Lambda dependent weight factors ===========================
  lambda_weight_1 = 1.0F / fourth_power(1.0F + quadruplet_lambda);
  lambda_weight_2 = 1.0F / fourth_power(1.0F - quadruplet_lambda);
  lambda_weight_3 = (2.0F * lambda_weight_1) * lambda_weight_2;

  // === Section: Directional indices =======================================
  const float theta_plus = std::abs(delta_theta4 / direction_increment);
  const int ithp = static_cast<int>(theta_plus);
  const int ithp1 = ithp + 1;
  const float wthp = theta_plus - static_cast<float>(ithp);
  const float wthp1 = 1.0F - wthp;

  const float theta_minus = std::abs(delta_theta3 / direction_increment);
  const int ithm = static_cast<int>(theta_minus);
  const int ithm1 = ithm + 1;
  const float wthm = theta_minus - static_cast<float>(ithm);
  const float wthm1 = 1.0F - wthm;

  // === Section: Frequency indices =========================================
  const FrequencyGeometry geometry =
      compute_frequency_geometry(frequency_ratio, quadruplet_lambda);
  const int ifrp = geometry.plus;
  const int ifrp1 = geometry.plus1;
  const int ifrm = geometry.minus;
  const int ifrm1 = geometry.minus1;

  // === Section: Range of calculations =====================================
  const int nfrchg = nfr - ifrm1;

  // === Section: Spectral addresses ========================================
  // Relative frequency addresses IF1 ... IF8 (size nfrchg) and directional
  // addresses IT1 ... IT8 (size nth), both 1-based as in the Fortran.
  std::vector<std::vector<int>> frequency_addresses(8);
  std::vector<std::vector<int>> direction_addresses(8);
  for (auto &addresses : frequency_addresses) {
    addresses.resize(static_cast<std::size_t>(nfrchg));
  }
  for (auto &addresses : direction_addresses) {
    addresses.resize(num_directions);
  }
  // IFn(ifr) and ITn(ith) with 1-based ifr, ith and n.
  const auto frequency_address = [&frequency_addresses](const int number,
                                                        const int ifr) -> int {
    return frequency_addresses[static_cast<std::size_t>(number - 1)]
                              [static_cast<std::size_t>(ifr - 1)];
  };
  const auto direction_address = [&direction_addresses](const int number,
                                                        const int ith) -> int {
    return direction_addresses[static_cast<std::size_t>(number - 1)]
                              [static_cast<std::size_t>(ith - 1)];
  };

  for (int ifr = 1; ifr <= nfrchg; ++ifr) {
    const auto at = static_cast<std::size_t>(ifr - 1);
    frequency_addresses[0][at] = ifr + ifrp;
    frequency_addresses[1][at] = ifr + ifrp1;
    frequency_addresses[2][at] = std::max(0, ifr + ifrm);
    frequency_addresses[3][at] = std::max(0, ifr + ifrm1);
    frequency_addresses[4][at] = std::max(0, ifr - ifrp);
    frequency_addresses[5][at] = std::max(0, ifr - ifrp1);
    frequency_addresses[6][at] = ifr - ifrm;
    frequency_addresses[7][at] = ifr - ifrm1;
  }

  for (int ith = 1; ith <= nth; ++ith) {
    const auto at = static_cast<std::size_t>(ith - 1);
    direction_addresses[0][at] = ith + ithp;
    direction_addresses[1][at] = ith + ithp1;
    direction_addresses[2][at] = ith + ithm;
    direction_addresses[3][at] = ith + ithm1;
    direction_addresses[4][at] = ith - ithp;
    direction_addresses[5][at] = ith - ithp1;
    direction_addresses[6][at] = ith - ithm;
    direction_addresses[7][at] = ith - ithm1;
    // Wrap around the closed directional space
    for (std::size_t number = 0; number < 4; ++number) {
      if (direction_addresses[number][at] > nth) {
        direction_addresses[number][at] -= nth;
      }
    }
    for (std::size_t number = 4; number < 8; ++number) {
      if (direction_addresses[number][at] < 1) {
        direction_addresses[number][at] += nth;
      }
    }
  }

  // Addresses of the interacting components (extended spectrum)
  for (std::size_t element = 0; element < num_interaction_bins; ++element) {
    const int isp = static_cast<int>(element) + 1;
    const int ifr = 1 + ((isp - 1) / nth);
    const int ith = 1 + ((isp - 1) % nth);
    plus_address_q1_c1[element] =
        direction_address(2, ith) + ((frequency_address(2, ifr) - 1) * nth);
    plus_address_q1_c2[element] =
        direction_address(1, ith) + ((frequency_address(2, ifr) - 1) * nth);
    plus_address_q1_c3[element] =
        direction_address(2, ith) + ((frequency_address(1, ifr) - 1) * nth);
    plus_address_q1_c4[element] =
        direction_address(1, ith) + ((frequency_address(1, ifr) - 1) * nth);
    minus_address_q1_c1[element] =
        direction_address(8, ith) + ((frequency_address(4, ifr) - 1) * nth);
    minus_address_q1_c2[element] =
        direction_address(7, ith) + ((frequency_address(4, ifr) - 1) * nth);
    minus_address_q1_c3[element] =
        direction_address(8, ith) + ((frequency_address(3, ifr) - 1) * nth);
    minus_address_q1_c4[element] =
        direction_address(7, ith) + ((frequency_address(3, ifr) - 1) * nth);
    plus_address_q2_c1[element] =
        direction_address(6, ith) + ((frequency_address(2, ifr) - 1) * nth);
    plus_address_q2_c2[element] =
        direction_address(5, ith) + ((frequency_address(2, ifr) - 1) * nth);
    plus_address_q2_c3[element] =
        direction_address(6, ith) + ((frequency_address(1, ifr) - 1) * nth);
    plus_address_q2_c4[element] =
        direction_address(5, ith) + ((frequency_address(1, ifr) - 1) * nth);
    minus_address_q2_c1[element] =
        direction_address(4, ith) + ((frequency_address(4, ifr) - 1) * nth);
    minus_address_q2_c2[element] =
        direction_address(3, ith) + ((frequency_address(4, ifr) - 1) * nth);
    minus_address_q2_c3[element] =
        direction_address(4, ith) + ((frequency_address(3, ifr) - 1) * nth);
    minus_address_q2_c4[element] =
        direction_address(3, ith) + ((frequency_address(3, ifr) - 1) * nth);
  }

  // Addresses for folding the interactions back to the discrete spectrum
  for (std::size_t element = 0; element < num_spectrum; ++element) {
    const int isp = static_cast<int>(element) + 1;
    const int ifr = 1 + ((isp - 1) / nth);
    const int ith = 1 + ((isp - 1) % nth);
    fold_address_c1_q1[element] =
        direction_address(6, ith) + ((frequency_address(6, ifr) - 1) * nth);
    fold_address_c1_q2[element] =
        direction_address(2, ith) + ((frequency_address(6, ifr) - 1) * nth);
    fold_address_c2_q1[element] =
        direction_address(5, ith) + ((frequency_address(6, ifr) - 1) * nth);
    fold_address_c2_q2[element] =
        direction_address(1, ith) + ((frequency_address(6, ifr) - 1) * nth);
    fold_address_c3_q1[element] =
        direction_address(6, ith) + ((frequency_address(5, ifr) - 1) * nth);
    fold_address_c3_q2[element] =
        direction_address(2, ith) + ((frequency_address(5, ifr) - 1) * nth);
    fold_address_c4_q1[element] =
        direction_address(5, ith) + ((frequency_address(5, ifr) - 1) * nth);
    fold_address_c4_q2[element] =
        direction_address(1, ith) + ((frequency_address(5, ifr) - 1) * nth);
    fold_address_c5_q1[element] =
        direction_address(4, ith) + ((frequency_address(8, ifr) - 1) * nth);
    fold_address_c5_q2[element] =
        direction_address(8, ith) + ((frequency_address(8, ifr) - 1) * nth);
    fold_address_c6_q1[element] =
        direction_address(3, ith) + ((frequency_address(8, ifr) - 1) * nth);
    fold_address_c6_q2[element] =
        direction_address(7, ith) + ((frequency_address(8, ifr) - 1) * nth);
    fold_address_c7_q1[element] =
        direction_address(4, ith) + ((frequency_address(7, ifr) - 1) * nth);
    fold_address_c7_q2[element] =
        direction_address(8, ith) + ((frequency_address(7, ifr) - 1) * nth);
    fold_address_c8_q1[element] =
        direction_address(3, ith) + ((frequency_address(7, ifr) - 1) * nth);
    fold_address_c8_q2[element] =
        direction_address(7, ith) + ((frequency_address(7, ifr) - 1) * nth);
  }

  // === Section: Scaling array (f**11) =====================================
  for (int ifr = 1; ifr <= nfr; ++ifr) {
    const auto at = static_cast<std::size_t>(ifr - 1);
    const float scaling = eleventh_power(radian_frequency[at] * inv_two_pi);
    for (int ith = 1; ith <= nth; ++ith) {
      frequency_scaling[static_cast<std::size_t>(ith + ((ifr - 1) * nth) - 1)] =
          scaling;
    }
  }

  float frequency = radian_frequency[num_frequencies - 1] * inv_two_pi;
  for (int ifr = nfr + 1; ifr <= nfrchg; ++ifr) {
    frequency = frequency * frequency_ratio;
    const float scaling = eleventh_power(frequency);
    for (int ith = 1; ith <= nth; ++ith) {
      frequency_scaling[static_cast<std::size_t>(ith + ((ifr - 1) * nth) - 1)] =
          scaling;
    }
  }

  // === Section: Interpolation weights =====================================
  interp_weight_1 = wthp * geometry.plus_weight;
  interp_weight_2 = wthp1 * geometry.plus_weight;
  interp_weight_3 = wthp * geometry.plus1_weight;
  interp_weight_4 = wthp1 * geometry.plus1_weight;
  interp_weight_5 = wthm * geometry.minus_weight;
  interp_weight_6 = wthm1 * geometry.minus_weight;
  interp_weight_7 = wthm * geometry.minus1_weight;
  interp_weight_8 = wthm1 * geometry.minus1_weight;

  squared_weight_1 = square(interp_weight_1);
  squared_weight_2 = square(interp_weight_2);
  squared_weight_3 = square(interp_weight_3);
  squared_weight_4 = square(interp_weight_4);
  squared_weight_5 = square(interp_weight_5);
  squared_weight_6 = square(interp_weight_6);
  squared_weight_7 = square(interp_weight_7);
  squared_weight_8 = square(interp_weight_8);
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

// --- insnl1_ranges_cpp ------------------------------------------------------
// Documented in NL1_ww3.h.
void insnl1_ranges_cpp(const int nk, const int nth, const float xfr,
                       const float lam, int *const nfrhgh, int *const nfrchg,
                       int *const nspecy, int *const nspecx) noexcept {
  std::size_t num_extended_frequencies = 0;
  std::size_t num_interaction_frequencies = 0;
  std::size_t num_auxiliary_bins = 0;
  std::size_t num_interaction_bins = 0;
  ww4_nl1::compute_insnl1_ranges(
      static_cast<std::size_t>(nk), static_cast<std::size_t>(nth), xfr, lam,
      num_extended_frequencies, num_interaction_frequencies, num_auxiliary_bins,
      num_interaction_bins);
  *nfrhgh = static_cast<int>(num_extended_frequencies);
  *nfrchg = static_cast<int>(num_interaction_frequencies);
  *nspecy = static_cast<int>(num_auxiliary_bins);
  *nspecx = static_cast<int>(num_interaction_bins);
}

// --- insnl1_cpp -------------------------------------------------------------
// Documented in NL1_ww3.h.
void insnl1_cpp(int nk, int nth, float dth, float xfr, float lam, float tpiinv,
                int nspecx, const float *sig, float *dal1, float *dal2,
                float *dal3, int *ip11, int *ip12, int *ip13, int *ip14,
                int *im11, int *im12, int *im13, int *im14, int *ip21,
                int *ip22, int *ip23, int *ip24, int *im21, int *im22,
                int *im23, int *im24, int *ic11, int *ic12, int *ic21,
                int *ic22, int *ic31, int *ic32, int *ic41, int *ic42,
                int *ic51, int *ic52, int *ic61, int *ic62, int *ic71,
                int *ic72, int *ic81, int *ic82, float *af11, float *awg1,
                float *awg2, float *awg3, float *awg4, float *awg5, float *awg6,
                float *awg7, float *awg8, float *swg1, float *swg2, float *swg3,
                float *swg4, float *swg5, float *swg6, float *swg7,
                float *swg8) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  const auto nspecz = nthz * nkz;
  const auto nspecxz = static_cast<std::size_t>(nspecx);
  ww4_nl1::compute_insnl1(
      std::span<const float>(sig, nkz), nthz, dth, xfr, lam, tpiinv, *dal1,
      *dal2, *dal3, std::span<int>(ip11, nspecxz),
      std::span<int>(ip12, nspecxz), std::span<int>(ip13, nspecxz),
      std::span<int>(ip14, nspecxz), std::span<int>(im11, nspecxz),
      std::span<int>(im12, nspecxz), std::span<int>(im13, nspecxz),
      std::span<int>(im14, nspecxz), std::span<int>(ip21, nspecxz),
      std::span<int>(ip22, nspecxz), std::span<int>(ip23, nspecxz),
      std::span<int>(ip24, nspecxz), std::span<int>(im21, nspecxz),
      std::span<int>(im22, nspecxz), std::span<int>(im23, nspecxz),
      std::span<int>(im24, nspecxz), std::span<int>(ic11, nspecz),
      std::span<int>(ic12, nspecz), std::span<int>(ic21, nspecz),
      std::span<int>(ic22, nspecz), std::span<int>(ic31, nspecz),
      std::span<int>(ic32, nspecz), std::span<int>(ic41, nspecz),
      std::span<int>(ic42, nspecz), std::span<int>(ic51, nspecz),
      std::span<int>(ic52, nspecz), std::span<int>(ic61, nspecz),
      std::span<int>(ic62, nspecz), std::span<int>(ic71, nspecz),
      std::span<int>(ic72, nspecz), std::span<int>(ic81, nspecz),
      std::span<int>(ic82, nspecz), std::span<float>(af11, nspecxz), *awg1,
      *awg2, *awg3, *awg4, *awg5, *awg6, *awg7, *awg8, *swg1, *swg2, *swg3,
      *swg4, *swg5, *swg6, *swg7, *swg8);
}

} // extern "C"
