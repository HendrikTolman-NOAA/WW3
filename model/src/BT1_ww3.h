/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file BT1_ww3.h
 * @brief Declarations of the single-point JONSWAP bottom friction (BT1) source
 * term routine for WAVEWATCH IV (WW4).
 * @details Declares the C++20 BT1 routine for the JONSWAP bottom friction
 * source term (W3SBT1) for one spatial point, and its C-compatible interface
 * called from the WAVEWATCH III (WW3) Fortran module w3sbt1md.F90.
 * Implementations are in BT1_ww3.cpp. The ww4_bt1 routine uses descriptive
 * variable names. The extern "C" wrapper keeps the WW3 Fortran names (A, CG,
 * WN, DEPTH, S, D, ...) for easy cross-referencing with w3sbt1md.F90.
 *
 * Bottom friction according to the empirical JONSWAP formulation (Hasselmann
 * et al., 1973; Bouws and Komen, 1983). For the wave action density spectrum
 * N(k,theta) the source term is linear in N,
 * S_bt(k,theta) = D_bt(k) N(k,theta), with diagonal term
 * D_bt(k) = (C_bt / d) max[0, c_g k / sigma - 0.5], where
 * C_bt = 2 Gamma / g (SBTC1), d is the mean water depth, c_g the group
 * velocity, k the wavenumber and sigma the radian frequency. The empirical
 * constant Gamma is -0.038 m^2/s^3 (JONSWAP) or -0.067 m^2/s^3 (Bouws and
 * Komen, 1983). The bottom friction is set to zero for wavenumbers with
 * k d > 6 (deep water), so that D_bt(k) is non-zero only for the lowest
 * wavenumber bands, up to but not including the first band with k d > 6.
 *
 * Output spectrum: the source term is a wave action density source term
 * S_bt(k,theta) on the same (k,theta) spectral grid as the action density
 * spectrum N(k,theta) = F(k,theta) / sigma used internally by WW3, where
 * F(k,theta) is the variance (energy) density spectrum and sigma the radian
 * (intrinsic) frequency. The spectra are stored as 1-D (flattened) arrays in
 * Fortran order, so that the index of direction ith and wavenumber ik (both
 * 0-based) is is = ith + nth * ik, with nth the number of directions. As the
 * diagonal term D_bt depends on wavenumber only, it is the same for all
 * directions of a wavenumber band.
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

#ifndef WW4_BT1_WW3_H
#define WW4_BT1_WW3_H

// System
#include <cstddef>
#include <span>

namespace ww4_bt1 {

// --- compute_w3sbt1 ---------------------------------------------------------
/**
 * @brief Computes the BT1 JONSWAP bottom friction source term (W3SBT1).
 * @details Bottom friction S_bt = D_bt N with the empirical JONSWAP diagonal
 * term D_bt(k) = (C_bt / d) max[0, c_g k / sigma - 0.5], see file description.
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
void compute_w3sbt1(std::span<const float> action_density,
                    std::span<const float> group_velocity,
                    std::span<const float> wavenumber, float depth,
                    std::span<float> source_term,
                    std::span<float> diagonal_term, std::size_t num_directions,
                    std::span<const float> radian_frequency,
                    float friction_coefficient) noexcept;

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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
void w3sbt1_cpp(const float *A, const float *CG, const float *WN, float depth,
                float *S, float *D, int nth, int nk, const float *sig,
                float sbtc1) noexcept;

} // extern "C"

#endif // WW4_BT1_WW3_H
