/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file DB2_ww3.cpp
 * @brief Implementation of the single-point Thornton-Guza type depth-induced
 * breaking (DB2) source term routine for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the DB2 depth-induced breaking source
 * term for one spatial point, with a C-compatible interface for the WAVEWATCH
 * III (WW3) Fortran code. DB2 is the IWB = 2 option of the WW3 routine W3SDB1,
 * made a separate source term that is fully independent of DB1. The input and
 * output spectrum (wave action density N(k,theta)) and the sign of the source
 * term are described in DB2_ww3.h. Each routine has its full documentation
 * (details, contracts and heritage) in front of its definition, and a shorter
 * version with its declaration in DB2_ww3.h.
 *
 * References:
 * - Thornton, E. B. and R. T. Guza, 1983: Transformation of wave height
 *   distribution. J. Geophys. Res., 88, 5925-5938.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section DB1.
 *
 * @copyright © 2026 National Weather Service, National Oceanic and Atmospheric
 * Administration. WAVEWATCH IV (TM) and WW4 (TM) are trademarks of the National
 * Weather Service.
 * The National Weather Service (NWS) often uses Generative SI (GenSI, formerly
 * identified as GenAI) for code development and refactoring. Whenever GenSI is
 * used, NWS requires a full human review of code before it is added to its
 * repositories.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted and adapted from the IWB = 2 option of
 * WAVEWATCH III (WW3) w3sdb1md.F90 (W3SDB1) by J. H. G. M. Alves and Hendrik
 * L. Tolman, with contributions by A. Roland.
 * @date Initial: 2026-11-20
 * @date Last update: 2026-11-20
 */

// Local
#include "DB2_ww3.h"

// System
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>

namespace ww4_db2 {

// --- Constants --------------------------------------------------------------
/// Smallest positive normal double (TINY(1.d0) in WW3). The source term is
/// not computed for a mean wave energy below this value.
constexpr double kMinEnergy = std::numeric_limits<double>::min();

/// Small threshold THR = TINY(1.0) of the WW3 module W3PARALL (single
/// precision), converted to double for the comparison with the total energy.
constexpr auto kThreshold =
    static_cast<double>(std::numeric_limits<float>::min());

/// Ratio H_rms^2 / E = 8 for a Rayleigh distribution of wave heights.
constexpr float kRmsHeightFactor = 8.0F;

/// Single precision pi, as the WW3 constant PI (3.141592653589793 rounded to
/// single precision).
constexpr float kPi = std::numbers::pi_v<float>;
static_assert(kPi == 3.141592653589793F); // NOLINT(modernize-use-std-numbers)

/// Denominator 16 in the constant alpha sqrt(pi) / 16.
constexpr float kDenominator = 16.0F;

// --- compute_w3sdb2 ---------------------------------------------------------
/**
 * @brief Computes the DB2 Thornton-Guza type depth-induced breaking source
 * term (W3SDB1 option IWB = 2).
 * @details Source term S_db = D_db N with
 * D_db = alpha (sqrt(pi) / 16) f_m H_rms^3 / (d E), see file description in
 * DB2_ww3.h, including the sign of D_db. The total energy E is integrated in
 * double precision, first over directions and then over frequencies.
 * H_rms = sqrt(8 E_m) and the factor alpha sqrt(pi) / 16 f_m are evaluated in
 * single precision, the rest of D_db in double precision, and D_db is rounded
 * to single precision before S_db = D_db N is computed. The order of
 * operations matches the Fortran original so results are reproducible.
 *
 * If the mean wave energy is below TINY(1.d0) (in practice zero), W3SDB1
 * returns without setting S, D and LBREAK, and the arrays of the caller keep
 * their values. This routine does the same: it returns an empty optional and
 * does not write to source_term and diagonal_term.
 *
 * Differences from the IWB = 2 code in W3SDB1: the maximum wave height, the
 * breaking fraction (DB1 only) and the factor FAK, which was computed but not
 * used, are not computed.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nspec = nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] depth Mean water depth d.
 * @param[in] mean_energy Mean wave energy E_m (EMEAN).
 * @param[in] mean_frequency Mean frequency f_m in Hz (FMEAN).
 * @param[out] source_term Depth-induced breaking source term S_db as 1-D array
 * (size nspec).
 * @param[out] diagonal_term Diagonal term D_db of the source term derivative
 * as 1-D array (size nspec).
 * @param[in] num_directions Number of directions nth.
 * @param[in] integration_factor Factor DDEN = d_sigma d_theta sigma of each
 * frequency band for the energy integral (size nk).
 * @param[in] proportionality_constant Proportionality constant alpha (SDBC1).
 * @return Empty if the mean wave energy is zero (outputs not changed),
 * otherwise true if waves are breaking (LBREAK).
 * @pre group_velocity and integration_factor have the same size nk,
 * action_density, source_term and diagonal_term have size nth * nk with
 * nth = num_directions.
 * @post If a value is returned, source_term and diagonal_term are fully
 * populated.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from the IWB = 2 option of WW3 W3SDB1 by
 * J. H. G. M. Alves and Hendrik L. Tolman, with contributions by A. Roland.
 * @date Initial: 2026-11-20
 * @date Last update: 2026-11-20
 */
std::optional<bool>
compute_w3sdb2(const std::span<const float> action_density,
               const std::span<const float> group_velocity, const float depth,
               const float mean_energy, const float mean_frequency,
               const std::span<float> source_term,
               const std::span<float> diagonal_term,
               const std::size_t num_directions,
               const std::span<const float> integration_factor,
               const float proportionality_constant) noexcept {
  // === Section: No energy, nothing to compute =============================
  if (static_cast<double>(mean_energy) < kMinEnergy) {
    return std::nullopt;
  }

  // === Section: Total energy ==============================================
  const std::size_t num_wavenumbers = group_velocity.size();
  double total_energy = 0.0;
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    double band_energy = 0.0;
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      band_energy +=
          static_cast<double>(action_density[ith + (num_directions * ik)]);
    }
    band_energy = (band_energy * static_cast<double>(integration_factor[ik])) /
                  static_cast<double>(group_velocity[ik]);
    total_energy += band_energy;
  }

  // === Section: Breaking coefficient ======================================
  double coefficient = 0.0;
  if (total_energy > kThreshold) {
    const auto rms_height =
        static_cast<double>(std::sqrt(kRmsHeightFactor * mean_energy));
    const float constant =
        ((proportionality_constant * std::sqrt(kPi)) / kDenominator) *
        mean_frequency;
    const double rms_height_cubed = (rms_height * rms_height) * rms_height;
    coefficient = -(((static_cast<double>(constant) * rms_height_cubed) /
                     static_cast<double>(depth)) /
                    total_energy);
  }

  // === Section: Source term ===============================================
  const auto diagonal_value = static_cast<float>(-coefficient);
  std::ranges::fill(diagonal_term, diagonal_value);
  std::ranges::transform(
      action_density, source_term.begin(),
      [diagonal_value](const float value) { return diagonal_value * value; });

  return coefficient > 0.0;
}

} // namespace ww4_db2

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sdb2_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_db2::compute_w3sdb2.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_db2::compute_w3sdb2. For binding as BIND(C, name="w3sdb2_cpp"). The
 * breaking flag is returned as an integer, as the Fortran LBREAK is a default
 * LOGICAL, and should be set in Fortran only if the return value is not
 * negative.
 * @param[in] A Action density spectrum as 1-D array (size nth * nk).
 * @param[in] CG Group velocities (size nk).
 * @param[in] depth Mean water depth.
 * @param[in] emean Mean wave energy.
 * @param[in] fmean Mean frequency (Hz).
 * @param[out] S Depth-induced breaking source term as 1-D array (size
 * nth * nk).
 * @param[out] D Diagonal term of the derivative as 1-D array (size nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] dden Integration factor DDEN (size nk).
 * @param[in] sdbc1 Proportionality constant (SDBC1).
 * @return -1 if nothing was computed (S, D and LBREAK not changed), 1 if waves
 * are breaking, 0 otherwise.
 * @pre nth > 0, nk > 0 and all arrays are allocated with the sizes given
 * above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-20
 * @date Last update: 2026-11-20
 */
int w3sdb2_cpp(const float *const A, const float *const CG, const float depth,
               const float emean, const float fmean, float *const S,
               float *const D, const int nth, const int nk,
               const float *const dden, const float sdbc1) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  const std::optional<bool> is_breaking = ww4_db2::compute_w3sdb2(
      std::span<const float>(A, nthz * nkz), std::span<const float>(CG, nkz),
      depth, emean, fmean, std::span<float>(S, nthz * nkz),
      std::span<float>(D, nthz * nkz), nthz, std::span<const float>(dden, nkz),
      sdbc1);
  if (!is_breaking.has_value()) {
    return -1;
  }
  return *is_breaking ? 1 : 0;
}

} // extern "C"
