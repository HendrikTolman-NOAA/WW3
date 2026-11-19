/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file LN1_ww3.cpp
 * @brief Implementation of the single-point linear wind input (LN1) source
 * term routine for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the LN1 linear wind input source term
 * (W3SLN1) for one spatial point, with a C-compatible interface called from
 * the WAVEWATCH III (WW3) Fortran module w3sln1md.F90. Starting point of the
 * WW4 conversion of the LN1 source term. The routine, the output spectrum
 * (wave action density source term S(k,theta)) and the journal references are
 * documented in LN1_ww3.h.
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
 * w3sln1md.F90 (W3SLN1) by Hendrik L. Tolman.
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */

// Local
#include "LN1_ww3.h"

// System
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace ww4_ln1 {

// --- Constants --------------------------------------------------------------
/// Wind speed ratio U_5 / u_* = 28 of the WAM-3 wind input, used in the
/// Pierson-Moskowitz frequency of the filter (alpha_PM g / (28 u_*)).
constexpr float kWindSpeedRatio = 28.0F;

/// The filter frequency is capped at this multiple of the highest discrete
/// radian frequency.
constexpr float kMaxFilterFrequencyRatio = 2.0F;

/// Lower limit on sigma / sigma_filt below which the filter G is set to zero.
constexpr float kMinFrequencyRatio = 0.5F;

// --- compute_w3sln1 ---------------------------------------------------------
// Documented in LN1_ww3.h.
void compute_w3sln1(
    const std::span<const float> wavenumber, const float cutoff_frequency,
    const float friction_velocity, const float wind_direction,
    const std::span<float> source_term, const std::size_t num_directions,
    const std::span<const float> cos_direction,
    const std::span<const float> sin_direction,
    const std::span<const float> radian_frequency,
    const float linear_input_coefficient, const float pm_frequency_factor,
    const float high_frequency_factor, const float gravity) noexcept {
  const std::size_t num_wavenumbers = wavenumber.size();

  // === Section: Set up factors ============================================
  const float cos_wind_dir = std::cos(wind_direction);
  const float sin_wind_dir = std::sin(wind_direction);

  // Directional factor max[0, cos(theta - theta_w)]^4
  std::vector<float> directional_factor(num_directions);
  for (std::size_t ith = 0; ith < num_directions; ++ith) {
    const float projection =
        std::max(0.0F, (cos_direction[ith] * cos_wind_dir) +
                           (sin_direction[ith] * sin_wind_dir));
    const float projection_sq = projection * projection;
    directional_factor[ith] = projection_sq * projection_sq;
  }

  // Amplitude and filter frequency sigma_filt
  const float friction_velocity_sq = friction_velocity * friction_velocity;
  const float amplitude =
      linear_input_coefficient * (friction_velocity_sq * friction_velocity_sq);
  const float pm_frequency =
      (pm_frequency_factor * gravity) / (kWindSpeedRatio * friction_velocity);
  const float cutoff_frequency_limited =
      high_frequency_factor *
      std::min(radian_frequency[num_wavenumbers - 1], cutoff_frequency);
  const float filter_frequency = std::min(
      std::max(pm_frequency, cutoff_frequency_limited),
      kMaxFilterFrequencyRatio * radian_frequency[num_wavenumbers - 1]);

  // === Section: Frequency dependent factor ================================
  std::vector<float> frequency_factor(num_wavenumbers);
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    const float frequency_ratio = radian_frequency[ik] / filter_frequency;
    if (frequency_ratio < kMinFrequencyRatio) {
      frequency_factor[ik] = 0.0F;
    } else {
      const float ratio_sq = frequency_ratio * frequency_ratio;
      const float ratio_4 = ratio_sq * ratio_sq;
      frequency_factor[ik] =
          (amplitude / wavenumber[ik]) * std::exp(-(1.0F / ratio_4));
    }
  }

  // === Section: Compose source term =======================================
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      source_term[ith + (num_directions * ik)] =
          frequency_factor[ik] * directional_factor[ith];
    }
  }
}

} // namespace ww4_ln1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sln1_cpp -------------------------------------------------------------
// Documented in LN1_ww3.h.
void w3sln1_cpp(const float *const K, const float fhigh, const float ustar,
                const float usdir, float *const S, const int nth, const int nk,
                const float *const ecos, const float *const esin,
                const float *const sig, const float slnc1, const float fspm,
                const float fshf, const float grav) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  ww4_ln1::compute_w3sln1(
      std::span<const float>(K, nkz), fhigh, ustar, usdir,
      std::span<float>(S, nthz * nkz), nthz, std::span<const float>(ecos, nthz),
      std::span<const float>(esin, nthz), std::span<const float>(sig, nkz),
      slnc1, fspm, fshf, grav);
}

} // extern "C"
