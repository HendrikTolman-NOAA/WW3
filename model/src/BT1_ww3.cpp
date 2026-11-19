/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file BT1_ww3.cpp
 * @brief Implementation of the single-point JONSWAP bottom friction (BT1)
 * source term routine for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the BT1 JONSWAP bottom friction source
 * term (W3SBT1) for one spatial point, with a C-compatible interface called
 * from the WAVEWATCH III (WW3) Fortran module w3sbt1md.F90. Starting point of
 * the WW4 conversion of the BT1 source term. The output spectrum (wave action
 * density source term S(k,theta)) is described in BT1_ww3.h. Each routine has
 * its full documentation (details, contracts and heritage) in front of its
 * definition, and a shorter version with its declaration in BT1_ww3.h.
 *
 * References:
 * - Hasselmann, K., et al., 1973: Measurements of wind-wave growth and swell
 *   decay during the Joint North Sea Wave Project (JONSWAP). Dtsch. Hydrogr.
 *   Z., Suppl. A, 8(12), 95 pp.
 * - Bouws, E. and G. J. Komen, 1983: On the balance between growth and
 *   dissipation in an extreme depth-limited wind-sea in the southern North
 *   Sea. J. Phys. Oceanogr., 13, 1653-1658.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section BT1.
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
 * w3sbt1md.F90 (W3SBT1) by Hendrik L. Tolman.
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */

// Local
#include "BT1_ww3.h"

// System
#include <algorithm>
#include <cstddef>
#include <span>

namespace ww4_bt1 {

// --- Constants --------------------------------------------------------------
/// Limit on k d above which the water is considered deep and the bottom
/// friction is set to zero.
constexpr float kDeepWaterLimit = 6.0F;

/// Offset in the group velocity to phase velocity ratio factor
/// max[0, c_g k / sigma - 0.5].
constexpr float kRatioOffset = 0.5F;

// --- compute_w3sbt1 ---------------------------------------------------------
/**
 * @brief Computes the BT1 JONSWAP bottom friction source term (W3SBT1).
 * @details Bottom friction S_bt = D_bt N with the empirical JONSWAP diagonal
 * term D_bt(k) = (C_bt / d) max[0, c_g k / sigma - 0.5], see file description.
 * If k d > 6 for the lowest wavenumber, the water is deep for the entire
 * spectrum and source_term and diagonal_term are set to zero. Otherwise
 * D_bt is evaluated for the wavenumber bands up to, but not including, the
 * first band with k d > 6 and set to zero for all higher bands. The order of
 * operations matches the Fortran original so results are reproducible.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] wavenumber Wavenumber k of each frequency band (size nk).
 * @param[in] depth Mean water depth d.
 * @param[out] source_term Bottom friction source term S_bt in action density
 * form as 1-D array (size nth * nk).
 * @param[out] diagonal_term Diagonal term D_bt of the source term derivative
 * as 1-D array (size nth * nk).
 * @param[in] num_directions Number of directions nth.
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] friction_coefficient Combined friction constant 2 Gamma / g
 * (SBTC1).
 * @pre group_velocity, wavenumber and radian_frequency have the same size
 * nk > 0, depth > 0, and action_density, source_term and diagonal_term have
 * size nth * nk with nth = num_directions.
 * @post The arrays source_term and diagonal_term are fully populated.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SBT1 by Hendrik L. Tolman.
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
void compute_w3sbt1(const std::span<const float> action_density,
                    const std::span<const float> group_velocity,
                    const std::span<const float> wavenumber, const float depth,
                    const std::span<float> source_term,
                    const std::span<float> diagonal_term,
                    const std::size_t num_directions,
                    const std::span<const float> radian_frequency,
                    const float friction_coefficient) noexcept {
  const std::size_t num_wavenumbers = wavenumber.size();

  // === Section: Deep water ================================================
  if (depth * wavenumber[0] > kDeepWaterLimit) {
    std::ranges::fill(diagonal_term, 0.0F);
    std::ranges::fill(source_term, 0.0F);
    return;
  }

  // === Section: Shallow water =============================================
  const float factor = friction_coefficient / depth;

  // Wavenumber dependent part, and fill of the diagonal term. Bands from the
  // first band with k d > 6 upward get a zero diagonal term.
  bool is_deep_band = false;
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    is_deep_band = is_deep_band || (wavenumber[ik] * depth > kDeepWaterLimit);
    const std::size_t first = num_directions * ik;
    const float diagonal_value =
        is_deep_band
            ? 0.0F
            : factor * std::max(0.0F, ((group_velocity[ik] * wavenumber[ik]) /
                                       radian_frequency[ik]) -
                                          kRatioOffset);
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      diagonal_term[first + ith] = diagonal_value;
      source_term[first + ith] = diagonal_value * action_density[first + ith];
    }
  }
}

} // namespace ww4_bt1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sbt1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_bt1::compute_w3sbt1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_bt1::compute_w3sbt1. Bound in w3sbt1md.F90 as
 * BIND(C, name="w3sbt1_cpp").
 * @param[in] A Action density spectrum as 1-D array (size nth * nk).
 * @param[in] CG Group velocities (size nk).
 * @param[in] WN Wavenumbers (size nk).
 * @param[in] depth Mean water depth.
 * @param[out] S Bottom friction source term as 1-D array (size nth * nk).
 * @param[out] D Diagonal term of the derivative as 1-D array (size nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] sbtc1 Combined friction constant 2 Gamma / g (SBTC1).
 * @pre nk > 0 and all arrays are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
void w3sbt1_cpp(const float *const A, const float *const CG,
                const float *const WN, const float depth, float *const S,
                float *const D, const int nth, const int nk,
                const float *const sig, const float sbtc1) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  ww4_bt1::compute_w3sbt1(
      std::span<const float>(A, nthz * nkz), std::span<const float>(CG, nkz),
      std::span<const float>(WN, nkz), depth, std::span<float>(S, nthz * nkz),
      std::span<float>(D, nthz * nkz), nthz, std::span<const float>(sig, nkz),
      sbtc1);
}

} // extern "C"
