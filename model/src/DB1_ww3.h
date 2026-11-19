/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file DB1_ww3.h
 * @brief Declarations of the single-point Battjes-Janssen depth-induced
 * breaking (DB1) source term routine for WAVEWATCH IV (WW4).
 * @details Declares the C++20 DB1 routine for the depth-induced breaking
 * source term (W3SDB1) for one spatial point, and its C-compatible interface
 * called from the WAVEWATCH III (WW3) Fortran module w3sdb1md.F90.
 * Implementations are in DB1_ww3.cpp. The ww4_db1 routine uses descriptive
 * variable names. The extern "C" wrapper keeps the WW3 Fortran names (A, DEPTH,
 * EMEAN, FMEAN, WNMEAN, S, D, ...) for easy cross-referencing with
 * w3sdb1md.F90.
 *
 * Depth-induced breaking according to the bore model of Battjes and Janssen
 * (1978). The source term is linear in the spectrum,
 * S_db(k,theta) = D_db N(k,theta), with a diagonal term D_db that is the same
 * for all spectral components. The maximum wave height H_max is either
 * H_max = gamma d (McCowan-type limit, FDONLY true), or
 * H_max = (gamma / k_m) tanh(k_m d) (Miche-type limit, FDONLY false), with
 * gamma the breaker index (SDBC2, BJGAM), d the water depth and k_m the mean
 * wavenumber. With H_rms = sqrt(8 E) and b^2 = H_rms^2 / H_max^2, the fraction
 * of breaking waves Q_b follows from (1 - Q_b) / (-ln Q_b) = b^2, and is
 * obtained by iteration. The diagonal term is
 * D_db = -2 alpha f_m Q_b / b^2 for b^2 < 1, and D_db = -2 alpha f_m b^2 in the
 * degenerate regime b^2 >= 1, where all waves are broken and the excess
 * energy is removed smoothly. Here alpha is the proportionality constant
 * (SDBC1, BJALFA), E the total wave energy and f_m the mean frequency. Q_b and
 * D_db are computed in double precision, as in WW3.
 *
 * Differences from the WW3 Fortran code: the alternative breaking formulation
 * of Thornton and Guza (1983) in W3SDB1 (IWB = 2) cannot be selected in WW3
 * (IWB is fixed to 1), and is not converted. The spectral integrals ETOT and
 * FMEAN2 of W3SDB1 are only used by that formulation and by debug output
 * (W3_DEBUGRUN), and are not computed.
 *
 * Input and output spectrum: the source term is a wave action density source
 * term S_db(k,theta) on the same (k,theta) spectral grid as the action density
 * spectrum N(k,theta) = F(k,theta) / sigma used internally by WW3, where
 * F(k,theta) is the variance (energy) density spectrum and sigma the radian
 * (intrinsic) frequency. As D_db does not depend on k or theta, the same D_db
 * applies to energy and action density spectra. The spectra are
 * stored as 1-D (flattened) arrays in Fortran order, so that the index of
 * direction ith and wavenumber ik (both 0-based) is is = ith + nth * ik, with
 * nth the number of directions.
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

#ifndef WW4_DB1_WW3_H
#define WW4_DB1_WW3_H

// System
#include <optional>
#include <span>

namespace ww4_db1 {

// --- compute_w3sdb1 ---------------------------------------------------------
/**
 * @brief Computes the DB1 Battjes-Janssen depth-induced breaking source term
 * (W3SDB1).
 * @details Source term S_db = D_db N with a diagonal term D_db that is the same
 * for all spectral components, see file description. If the mean wave energy
 * is zero, nothing is computed and the output arrays are not changed.
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
[[nodiscard]] std::optional<bool>
compute_w3sdb1(std::span<const float> action_density, float depth,
               float mean_energy, float mean_frequency, float mean_wavenumber,
               std::span<float> source_term, std::span<float> diagonal_term,
               float proportionality_constant, float breaker_index,
               bool depth_only) noexcept;

} // namespace ww4_db1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sdb1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_db1::compute_w3sdb1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_db1::compute_w3sdb1. Bound in w3sdb1md.F90 as
 * BIND(C, name="w3sdb1_cpp").
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-19
 * @date Last update: 2026-11-19
 */
int w3sdb1_cpp(const float *A, float depth, float emean, float fmean,
               float wnmean, float *S, float *D, int nspec, float sdbc1,
               float sdbc2, int fdonly) noexcept;

} // extern "C"

#endif // WW4_DB1_WW3_H
