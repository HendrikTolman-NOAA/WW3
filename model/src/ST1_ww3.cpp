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
 * ST1 source terms. The input spectrum (wave action density
 * N(k,theta)) and the journal references are described in ST1_ww3.h. Each
 * routine has its full documentation (details, contracts and heritage) in
 * front of its definition, and a shorter version with its declaration in
 * ST1_ww3.h.
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
 * @date Last update: 2026-11-19
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
/**
 * @brief Computes the ST1 wind input source term (W3SIN1).
 * @details WAM-3 wind input of Snyder et al. (1981), as adapted by Komen et
 * al. (1984): S_in = D * N with diagonal term
 * D = C * sigma * max[0, (u_* / c) cos(theta - theta_w) - 1/28],
 * where c = sigma / k is the phase velocity and C = 28 C_in rho_a / rho_w.
 * Snyder et al. (1981) expressed the growth rate with the wind speed U_5 at
 * 5 m height, which Komen et al. (1984) replaced by U_5 = 28 u_*. Wind input
 * therefore only acts on components with c < 28 u_* cos(theta - theta_w).
 * @param[in] action_density Wave action density spectrum N(k,theta) in 1-D
 * (flattened) form (size nspec = nth * nk), see file description.
 * @param[in] wavenumber Wavenumber of each spectral bin (size nspec).
 * @param[in] friction_velocity Wind friction velocity u_*.
 * @param[in] cos_wind_dir Cosine of the wind (friction velocity) direction
 * theta_w.
 * @param[in] sin_wind_dir Sine of the wind (friction velocity) direction
 * theta_w.
 * @param[out] source_term Wind input source term in action density form (size
 * nspec).
 * @param[out] source_diagonal Diagonal term of the derivative of source_term
 * with respect to action_density (size nspec).
 * @param[in] input_coefficient Combined input constant 28 C_in rho_a / rho_w
 * (SINC1).
 * @param[in] radian_frequency Radian frequency sigma of each spectral bin
 * (SIG2, size nspec).
 * @param[in] cos_direction Cosine of the direction of each spectral bin (size
 * nspec).
 * @param[in] sin_direction Sine of the direction of each spectral bin (size
 * nspec).
 * @pre All spans have the same size nspec.
 * @post The arrays source_term and source_diagonal are populated with the
 * calculated input source terms.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SIN1 by Hendrik L. Tolman.
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
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
/**
 * @brief Computes the ST1 whitecapping dissipation source term (W3SDS1).
 * @details WAM-3 whitecapping of Komen et al. (1984), see also WAMDI Group
 * (1988): S_ds = D * N with diagonal term
 * D = C_ds sigma_hat (k / k_hat) (alpha_hat / alpha_PM)^2,
 * where alpha_hat = E k_hat^2. This reduces to
 * D = dissipation_coefficient * f_hat * k_hat^3 * E^2 * k.
 * @param[in] action_density Wave action density spectrum N(k,theta) in 1-D
 * (flattened) form (size nspec = nth * nk), see file description.
 * @param[in] wavenumber Wavenumber of each spectral bin (size nspec).
 * @param[in] total_variance Total wave variance (energy) E of the spectrum.
 * @param[in] mean_frequency Mean wave frequency f_hat (Hz).
 * @param[in] mean_wavenumber Mean wavenumber k_hat.
 * @param[out] source_term Dissipation source term in action density form (size
 * nspec).
 * @param[out] source_diagonal Diagonal term of the derivative of source_term
 * with respect to action_density (size nspec).
 * @param[in] dissipation_coefficient Combined dissipation constant
 * 2 pi C_ds / alpha_PM^2 (SDSC1).
 * @pre All spans have the same size nspec.
 * @post The arrays source_term and source_diagonal are populated with the
 * calculated dissipation source terms.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SDS1 by Hendrik L. Tolman.
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
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
/**
 * @brief Computes the mean wave parameters used by the ST1 source terms
 * (W3SPR1).
 * @details Calculates the total wave variance E, the mean frequency
 * f_hat = (1 / 2 pi) E / integral(F / sigma) and the mean wavenumber
 * k_hat = [E / integral(F k^-1/2)]^2 of the WAM-3 method (WAMDI Group, 1988),
 * plus the maximum action density, including the parametric tail beyond the
 * discrete spectrum. The summation order matches the Fortran original so
 * results are reproducible.
 * @param[in] action_density Wave action density spectrum N(k,theta) as a 2-D
 * array A(ith, ik) in Fortran (column-major) order (size nth * nk), see file
 * description.
 * @param[in] group_velocity Group velocities (size nk).
 * @param[in] wavenumber Wavenumbers (size nk).
 * @param[in] radian_frequency Radian frequencies (size nk).
 * @param[in] band_integration_factor Integration factors DTH * DSII * SIG,
 * converting action density to energy per frequency band (DDEN, size nk).
 * @param[in] num_directions Number of directions nth.
 * @param[in] tail_factor_energy Factor in tail integration of energy (FTE).
 * @param[in] tail_factor_frequency Factor in tail integration of frequency
 * (FTF).
 * @param[in] tail_factor_wavenumber Factor in tail integration of wavenumber
 * (FTWN).
 * @param[in] inv_two_pi Inverse of 2*Pi (TPIINV), passed from Fortran so that
 * results stay bit-identical to WW3.
 * @param[out] total_variance Total wave variance (energy) E.
 * @param[out] mean_frequency Mean wave frequency f_hat (Hz).
 * @param[out] mean_wavenumber Mean wavenumber k_hat.
 * @param[out] max_action_density Maximum action density in the spectrum.
 * @pre group_velocity, wavenumber, radian_frequency and band_integration_factor
 * have the same size nk > 0, and action_density has size nth * nk.
 * @post total_variance, mean_frequency, mean_wavenumber and max_action_density
 * are set.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SPR1 by Hendrik L. Tolman.
 * @date Initial: 2026-11-17
 * @date Last update: 2026-11-19
 */
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
/**
 * @brief Fortran interface to ww4_st1::compute_w3sin1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_st1::compute_w3sin1. Bound in w3src1md.F90 as
 * BIND(C, name="w3sin1_cpp").
 * @param[in] A Wave action density spectrum N(k,theta) in 1-D form (size
 * nspec).
 * @param[in] K Wavenumber of each spectral bin (size nspec).
 * @param[in] ustar Wind friction velocity u_*.
 * @param[in] cosu Cosine of the wind (friction velocity) direction.
 * @param[in] sinu Sine of the wind (friction velocity) direction.
 * @param[out] S Wind input source term (size nspec).
 * @param[out] D Diagonal term of the derivative (size nspec).
 * @param[in] nspec Number of spectral bins (nth * nk).
 * @param[in] sinc1 Combined input constant (SINC1).
 * @param[in] sig2 Radian frequency of each spectral bin (size nspec).
 * @param[in] ecos Cosine of the direction of each spectral bin (size nspec).
 * @param[in] esin Sine of the direction of each spectral bin (size nspec).
 * @pre All arrays are allocated with at least nspec elements.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
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
/**
 * @brief Fortran interface to ww4_st1::compute_w3sds1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_st1::compute_w3sds1. Bound in w3src1md.F90 as
 * BIND(C, name="w3sds1_cpp").
 * @param[in] A Wave action density spectrum N(k,theta) in 1-D form (size
 * nspec).
 * @param[in] K Wavenumber of each spectral bin (size nspec).
 * @param[in] emean Total wave variance (energy) E.
 * @param[in] fmean Mean wave frequency (Hz).
 * @param[in] wnmean Mean wavenumber.
 * @param[out] S Dissipation source term (size nspec).
 * @param[out] D Diagonal term of the derivative (size nspec).
 * @param[in] nspec Number of spectral bins (nth * nk).
 * @param[in] sdsc1 Combined dissipation constant (SDSC1).
 * @pre All arrays are allocated with at least nspec elements.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
void w3sds1_cpp(const float *const A, const float *const K, const float emean,
                const float fmean, const float wnmean, float *const S,
                float *const D, const int nspec, const float sdsc1) noexcept {
  const auto n = static_cast<std::size_t>(nspec);
  ww4_st1::compute_w3sds1(
      std::span<const float>(A, n), std::span<const float>(K, n), emean, fmean,
      wnmean, std::span<float>(S, n), std::span<float>(D, n), sdsc1);
}

// --- w3spr1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_st1::compute_w3spr1.
 * @details A is the Fortran array A(nth, nk). Output scalars are passed by
 * reference (Fortran default). Bound in w3src1md.F90 as
 * BIND(C, name="w3spr1_cpp").
 * @param[in] A Wave action density spectrum N(k,theta) as 2-D array
 * A(nth, nk) (size nth * nk).
 * @param[in] cg Group velocities (size nk).
 * @param[in] wn Wavenumbers (size nk).
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Integration factors DTH * DSII * SIG (size nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] fte Factor in tail integration of energy.
 * @param[in] ftf Factor in tail integration of frequency.
 * @param[in] ftwn Factor in tail integration of wavenumber.
 * @param[in] tpiinv Inverse of 2*Pi.
 * @param[out] emean Total wave variance (energy) E.
 * @param[out] fmean Mean wave frequency (Hz).
 * @param[out] wnmean Mean wavenumber.
 * @param[out] amax Maximum action density in the spectrum.
 * @pre nk > 0 and all arrays are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-17
 * @date Last update: 2026-11-19
 */
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
