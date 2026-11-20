/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file DB1_ww3.cpp
 * @brief Implementation of the single-point Battjes-Janssen depth-induced
 * breaking (DB1) source term routine for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the DB1 depth-induced breaking source
 * term (W3SDB1) for one spatial point, with a C-compatible interface called
 * from the WAVEWATCH III (WW3) Fortran module w3sdb1md.F90. Starting point of
 * the WW4 conversion of the DB1 source term. The input and output spectrum
 * (wave action density N(k,theta)) is described in DB1_ww3.h. Each routine has
 * its full documentation (details, contracts and heritage) in front of its
 * definition, and a shorter version with its declaration in DB1_ww3.h.
 *
 * References:
 * - Battjes, J. A. and J. P. F. M. Janssen, 1978: Energy loss and set-up due
 *   to breaking of random waves. Proc. 16th Int. Conf. Coastal Eng., ASCE,
 *   569-587.
 * - Miche, A., 1944: Mouvements ondulatoire de la mer en profondeur croissante
 *   ou decroissante. Forme limite de la houle lors de son deferlement.
 *   Application aux digues maritimes. Annales des Ponts et Chaussees, 114,
 *   25-78, 131-164, 270-292, 369-406.
 * - Eldeberky, Y. and J. A. Battjes, 1996: Spectral modelling of wave
 *   breaking: Application to Boussinesq equations. J. Geophys. Res., 101,
 *   1253-1264.
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
 * @author Code Heritage: Converted and adapted from WAVEWATCH III (WW3)
 * w3sdb1md.F90 (W3SDB1) by J. H. G. M. Alves and Hendrik L. Tolman, with
 * contributions by A. Roland.
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */

// Local
#include "DB1_ww3.h"

// System
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>

namespace ww4_db1 {

// --- Constants --------------------------------------------------------------
/// Smallest positive normal double (TINY(1.d0) in WW3). The source term is
/// not computed for a mean wave energy below this value.
constexpr double kMinEnergy = std::numeric_limits<double>::min();

/// Small threshold THR = TINY(1.0) of the WW3 module W3PARALL (single
/// precision).
constexpr float kThreshold = std::numeric_limits<float>::min();

/// kThreshold converted to double, as in the mixed-precision comparisons of
/// W3SDB1.
constexpr double kThresholdDouble = static_cast<double>(kThreshold);

/// Ratio H_rms^2 / E = 8 for a Rayleigh distribution of wave heights.
constexpr double kRmsHeightFactor = 8.0;

/// Below this ratio b = H_rms / H_max, the first guess of Q_b is zero.
constexpr double kFirstGuessLimit = 0.5;

/// Below this ratio b = H_rms / H_max, no waves are breaking (Q_b = 0).
constexpr double kNoBreakingLimit = 0.2;

/// Number of fixed-point iterations for Q_b.
constexpr int kNumIterations = 3;

/// Breaking fraction in the degenerate regime b >= 1. Computed in single
/// precision as in W3SDB1 (1.0 - THR), which gives exactly 1.
constexpr double kFullBreakingFraction = static_cast<double>(1.0F - kThreshold);

// --- compute_w3sdb1 ---------------------------------------------------------
/**
 * @brief Computes the DB1 Battjes-Janssen depth-induced breaking source term
 * (W3SDB1).
 * @details Source term S_db = D_db N with a diagonal term D_db that is the same
 * for all spectral components, see file description. The maximum wave height
 * H_max, the ratio b^2 = H_rms^2 / H_max^2, the fraction of breaking waves Q_b
 * and D_db are computed in double precision. Q_b is obtained from a first guess
 * Q_0 = (2 b - 1)^2 (zero for b <= 0.5), one Newton-type step and three
 * fixed-point iterations Q_b = exp((Q_b - 1) / b^2). Q_b = 0 for b <= 0.2, and
 * Q_b = 1 for b >= 1. D_db is zero if b^2 or |b^2 - Q_b| does not exceed the
 * single precision threshold THR = TINY(1.0). D_db is rounded to single
 * precision before S_db = D_db N is computed. The order of operations matches
 * the Fortran original so results are reproducible.
 *
 * If the mean wave energy is below TINY(1.d0) (in practice zero), W3SDB1
 * returns without setting S, D and LBREAK, and the arrays of the caller keep
 * their values. This routine does the same: it returns an empty optional and
 * does not write to source_term and diagonal_term.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nspec = nth * nk), see file description.
 * @param[in] depth Mean water depth d.
 * @param[in] mean_energy Total wave energy E (EMEAN).
 * @param[in] mean_frequency Mean frequency f_m in Hz (FMEAN).
 * @param[in] mean_wavenumber Mean wavenumber k_m (WNMEAN).
 * @param[out] source_term Depth-induced breaking source term S_db as 1-D array
 * (size nspec).
 * @param[out] diagonal_term Diagonal term D_db of the source term derivative
 * as 1-D array (size nspec).
 * @param[in] proportionality_constant Proportionality constant alpha (SDBC1).
 * @param[in] breaker_index Breaker index gamma, ratio H_max / d (SDBC2).
 * @param[in] depth_only Use the McCowan-type limit H_max = gamma d (FDONLY).
 * @return Empty if the mean wave energy is zero (outputs not changed),
 * otherwise true if waves are breaking (LBREAK).
 * @pre action_density, source_term and diagonal_term have the same size nspec.
 * mean_wavenumber > 0 if depth_only is false.
 * @post If a value is returned, source_term and diagonal_term are fully
 * populated.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SDB1 by J. H. G. M. Alves and
 * Hendrik L. Tolman, with contributions by A. Roland.
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
std::optional<bool>
compute_w3sdb1(const std::span<const float> action_density, const float depth,
               const float mean_energy, const float mean_frequency,
               const float mean_wavenumber, const std::span<float> source_term,
               const std::span<float> diagonal_term,
               const float proportionality_constant, const float breaker_index,
               const bool depth_only) noexcept {
  // === Section: No energy, nothing to compute =============================
  if (static_cast<double>(mean_energy) < kMinEnergy) {
    return std::nullopt;
  }

  // === Section: Maximum wave height =======================================
  const double max_height =
      depth_only
          ? static_cast<double>(breaker_index) * static_cast<double>(depth)
          : (static_cast<double>(breaker_index) /
             static_cast<double>(mean_wavenumber)) *
                std::tanh(static_cast<double>(mean_wavenumber) *
                          static_cast<double>(std::max(depth, 0.0F)));

  // === Section: H_rms and ratio to H_max ==================================
  const double rms_height =
      std::sqrt(kRmsHeightFactor * static_cast<double>(mean_energy));
  double ratio_sq = 0.0;
  double ratio = 0.0;
  if (max_height > kThresholdDouble) {
    ratio_sq = (rms_height * rms_height) / (max_height * max_height);
    ratio = std::sqrt(ratio_sq);
  }

  // === Section: Fraction of breaking waves ================================
  double breaking_fraction = 0.0;
  if (ratio <= kNoBreakingLimit) {
    breaking_fraction = 0.0;
  } else if (ratio < 1.0) {
    const double first_guess =
        ratio <= kFirstGuessLimit
            ? 0.0
            : ((2.0 * ratio) - 1.0) * ((2.0 * ratio) - 1.0);
    const double arg = std::exp((first_guess - 1.0) / ratio_sq);
    breaking_fraction =
        first_guess - ((ratio_sq * (first_guess - arg)) / (ratio_sq - arg));
    for (int iteration = 0; iteration < kNumIterations; ++iteration) {
      breaking_fraction = std::exp((breaking_fraction - 1.0) / ratio_sq);
    }
  } else {
    breaking_fraction = kFullBreakingFraction;
  }

  // === Section: Breaking coefficient ======================================
  double coefficient = 0.0;
  if ((ratio_sq > kThresholdDouble) &&
      (std::abs(ratio_sq - breaking_fraction) > kThresholdDouble)) {
    const auto alpha = static_cast<double>(proportionality_constant);
    const auto frequency = static_cast<double>(mean_frequency);
    if (ratio_sq < 1.0) {
      coefficient =
          (((2.0 * alpha) * breaking_fraction) * frequency) / ratio_sq;
    } else {
      // Degenerate regime: all waves are broken, excess energy is removed
      // smoothly by the factor b^2.
      coefficient = ((2.0 * alpha) * frequency) * ratio_sq;
    }
  }

  // === Section: Source term ===============================================
  const auto diagonal_value = static_cast<float>(-coefficient);
  std::ranges::fill(diagonal_term, diagonal_value);
  std::ranges::transform(
      action_density, source_term.begin(),
      [diagonal_value](const float value) { return diagonal_value * value; });

  return coefficient > 0.0;
}

} // namespace ww4_db1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sdb1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_db1::compute_w3sdb1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_db1::compute_w3sdb1. Bound in w3sdb1md.F90 as
 * BIND(C, name="w3sdb1_cpp"). The breaking flag is returned as an integer, as
 * the Fortran LBREAK is a default LOGICAL, and is set in Fortran only if the
 * return value is not negative.
 * @param[in] A Action density spectrum as 1-D array (size nspec).
 * @param[in] depth Mean water depth.
 * @param[in] emean Total wave energy.
 * @param[in] fmean Mean frequency (Hz).
 * @param[in] wnmean Mean wavenumber.
 * @param[out] S Depth-induced breaking source term as 1-D array (size nspec).
 * @param[out] D Diagonal term of the derivative as 1-D array (size nspec).
 * @param[in] nspec Number of spectral components nth * nk.
 * @param[in] sdbc1 Proportionality constant (SDBC1).
 * @param[in] sdbc2 Breaker index, ratio H_max / d (SDBC2).
 * @param[in] fdonly 1 for the depth-only (McCowan-type) limit, 0 for the
 * Miche-type limit (FDONLY).
 * @return -1 if nothing was computed (S, D and LBREAK not changed), 1 if waves
 * are breaking, 0 otherwise.
 * @pre nspec > 0 and all arrays are allocated with size nspec.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
int w3sdb1_cpp(const float *const A, const float depth, const float emean,
               const float fmean, const float wnmean, float *const S,
               float *const D, const int nspec, const float sdbc1,
               const float sdbc2, const int fdonly) noexcept {
  const auto nspecz = static_cast<std::size_t>(nspec);
  const std::optional<bool> is_breaking = ww4_db1::compute_w3sdb1(
      std::span<const float>(A, nspecz), depth, emean, fmean, wnmean,
      std::span<float>(S, nspecz), std::span<float>(D, nspecz), sdbc1, sdbc2,
      fdonly != 0);
  if (!is_breaking.has_value()) {
    return -1;
  }
  return *is_breaking ? 1 : 0;
}

} // extern "C"
