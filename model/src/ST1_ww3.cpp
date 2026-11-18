/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file ST1_ww3.cpp
 * @brief Implementation of the single-point WAM-3 (ST1) source term routines
 * for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementations of the ST1 mean wave parameters (W3SPR1),
 * wind input (W3SIN1) and whitecapping dissipation (W3SDS1) for one spatial
 * point, with C-compatible interfaces called from the WAVEWATCH III (WW3)
 * Fortran module w3src1md.F90. Starting point of the WW4 conversion of the
 * ST1 source terms. The routines, the input spectrum (wave action density
 * N(k,theta)) and the journal references are documented in ST1_ww3.h.
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
 * w3src1md.F90 (W3SPR1, W3SIN1, W3SDS1) by Hendrik L. Tolman.
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-18
 */

// Local
#include "ST1_ww3.h"

// System
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

namespace ww4_st1 {

// --- Constants --------------------------------------------------------------
/// Inverse of the wind speed ratio U_5 / u_* = 28 in the WAM-3 wind input
/// (1/28, value as in the WW3 Fortran code). With it, the growth factor
/// (u_* / c) cos(theta - theta_w) - 1/28 is positive only for spectral
/// components with 28 u_* cos(theta - theta_w) > c.
constexpr float kInvWindSpeedRatio = 0.035714F;

/// Lower limit on the frequency and wavenumber integrals in the mean parameters
/// (avoids division by zero).
constexpr float kMinIntegral = 1.0e-7F;

// --- compute_w3sin1 ---------------------------------------------------------
// Documented in ST1_ww3.h.
void compute_w3sin1(const std::span<const float> action_density,
                    const std::span<const float> wavenumber,
                    const float friction_velocity, const float cos_wind_dir,
                    const float sin_wind_dir,
                    const std::span<float> source_term,
                    const std::span<float> source_diagonal,
                    const float input_coefficient,
                    const std::span<const float> radian_frequency,
                    const std::span<const float> cos_direction,
                    const std::span<const float> sin_direction) noexcept {
  // === Section: Diagonal and source term ==================================
  for (std::size_t is = 0; is < action_density.size(); ++is) {
    // (u_* / c) cos(theta - theta_w) - 1/28, with c = sigma / k
    const float wind_forcing =
        (((friction_velocity * ((cos_direction[is] * cos_wind_dir) +
                                (sin_direction[is] * sin_wind_dir))) *
          wavenumber[is]) /
         radian_frequency[is]) -
        kInvWindSpeedRatio;
    source_diagonal[is] = (input_coefficient * radian_frequency[is]) *
                          std::max(0.0F, wind_forcing);
    source_term[is] = source_diagonal[is] * action_density[is];
  }
}

// --- compute_w3sds1 ---------------------------------------------------------
// Documented in ST1_ww3.h.
void compute_w3sds1(const std::span<const float> action_density,
                    const std::span<const float> wavenumber,
                    const float total_variance, const float mean_frequency,
                    const float mean_wavenumber,
                    const std::span<float> source_term,
                    const std::span<float> source_diagonal,
                    const float dissipation_coefficient) noexcept {
  // === Section: Common factor =============================================
  const float dissipation_scale =
      ((dissipation_coefficient * mean_frequency) *
       (mean_wavenumber * mean_wavenumber * mean_wavenumber)) *
      (total_variance * total_variance);

  // === Section: Diagonal and source term ==================================
  for (std::size_t is = 0; is < action_density.size(); ++is) {
    source_diagonal[is] = dissipation_scale * wavenumber[is];
    source_term[is] = source_diagonal[is] * action_density[is];
  }
}

// --- compute_w3spr1 ---------------------------------------------------------
// Documented in ST1_ww3.h.
void compute_w3spr1(const std::span<const float> action_density,
                    const std::span<const float> group_velocity,
                    const std::span<const float> wavenumber,
                    const std::span<const float> radian_frequency,
                    const std::span<const float> band_integration_factor,
                    const std::size_t num_directions,
                    const float tail_factor_energy,
                    const float tail_factor_frequency,
                    const float tail_factor_wavenumber, const float inv_two_pi,
                    float &total_variance, float &mean_frequency,
                    float &mean_wavenumber,
                    float &max_action_density) noexcept {
  const std::size_t num_wavenumbers = group_velocity.size();

  float energy_sum = 0.0F;
  float inv_frequency_sum = 0.0F;
  float inv_sqrt_wavenumber_sum = 0.0F;
  float max_action = 0.0F;
  float band_energy = 0.0F;

  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    // === Section: Integral over directions ==============================
    band_energy = 0.0F;
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      const float action = action_density[ith + (num_directions * ik)];
      band_energy += action;
      max_action = std::max(max_action, action);
    }
    // === Section: Integral over frequencies =============================
    band_energy =
        (band_energy * band_integration_factor[ik]) / group_velocity[ik];
    energy_sum += band_energy;
    inv_frequency_sum += band_energy / radian_frequency[ik];
    inv_sqrt_wavenumber_sum += band_energy / std::sqrt(wavenumber[ik]);
  }

  // === Section: Tail beyond discrete spectrum =============================
  // band_energy holds the energy of the last discrete band
  const float tail_energy =
      band_energy / band_integration_factor[num_wavenumbers - 1];
  energy_sum += tail_energy * tail_factor_energy;
  inv_frequency_sum += tail_energy * tail_factor_frequency;
  inv_sqrt_wavenumber_sum += tail_energy * tail_factor_wavenumber;

  // === Section: Final processing ==========================================
  mean_frequency =
      (inv_two_pi * energy_sum) / std::max(kMinIntegral, inv_frequency_sum);
  const float sqrt_mean_wavenumber =
      energy_sum / std::max(kMinIntegral, inv_sqrt_wavenumber_sum);
  mean_wavenumber = sqrt_mean_wavenumber * sqrt_mean_wavenumber;
  total_variance = energy_sum;
  max_action_density = max_action;
}

} // namespace ww4_st1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sin1_cpp -------------------------------------------------------------
// Documented in ST1_ww3.h.
void w3sin1_cpp(const float *const A, const float *const K, const float ustar,
                const float cosu, const float sinu, float *const S,
                float *const D, const int nspec, const float sinc1,
                const float *const sig2, const float *const ecos,
                const float *const esin) noexcept {
  const auto n = static_cast<std::size_t>(nspec);
  ww4_st1::compute_w3sin1(
      std::span<const float>(A, n), std::span<const float>(K, n), ustar, cosu,
      sinu, std::span<float>(S, n), std::span<float>(D, n), sinc1,
      std::span<const float>(sig2, n), std::span<const float>(ecos, n),
      std::span<const float>(esin, n));
}

// --- w3sds1_cpp -------------------------------------------------------------
// Documented in ST1_ww3.h.
void w3sds1_cpp(const float *const A, const float *const K, const float emean,
                const float fmean, const float wnmean, float *const S,
                float *const D, const int nspec, const float sdsc1) noexcept {
  const auto n = static_cast<std::size_t>(nspec);
  ww4_st1::compute_w3sds1(
      std::span<const float>(A, n), std::span<const float>(K, n), emean, fmean,
      wnmean, std::span<float>(S, n), std::span<float>(D, n), sdsc1);
}

// --- w3spr1_cpp -------------------------------------------------------------
// Documented in ST1_ww3.h.
void w3spr1_cpp(const float *const A, const float *const cg,
                const float *const wn, const float *const sig,
                const float *const dden, const int nth, const int nk,
                const float fte, const float ftf, const float ftwn,
                const float tpiinv, float *const emean, float *const fmean,
                float *const wnmean, float *const amax) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  ww4_st1::compute_w3spr1(
      std::span<const float>(A, nthz * nkz), std::span<const float>(cg, nkz),
      std::span<const float>(wn, nkz), std::span<const float>(sig, nkz),
      std::span<const float>(dden, nkz), nthz, fte, ftf, ftwn, tpiinv, *emean,
      *fmean, *wnmean, *amax);
}

} // extern "C"
