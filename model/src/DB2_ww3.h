/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file DB2_ww3.h
 * @brief Declarations of the single-point Thornton-Guza type depth-induced
 * breaking (DB2) source term routine for WAVEWATCH IV (WW4).
 * @details Declares the C++20 DB2 routine for the depth-induced breaking
 * source term for one spatial point, and its C-compatible interface for the
 * WAVEWATCH III (WW3) Fortran code. Implementations are in DB2_ww3.cpp. The
 * ww4_db2 routine uses descriptive variable names. The extern "C" wrapper uses
 * the WW3 Fortran names (A, CG, DEPTH, EMEAN, FMEAN, S, D, ...).
 *
 * DB2 is the second breaking formulation of the WW3 routine W3SDB1 (option
 * IWB = 2, based on Thornton and Guza, 1983). In WW3 this option could not be
 * selected, as IWB is fixed to 1 in W3SDB1. In WW4 it is a separate source
 * term, fully independent of DB1 (Battjes and Janssen, 1978), so that each can
 * be selected on its own. DB2 does not use any DB1 code.
 *
 * The source term is linear in the spectrum,
 * S_db(k,theta) = D_db N(k,theta), with a diagonal term D_db that is the same
 * for all spectral components:
 * D_db = alpha (sqrt(pi) / 16) f_m H_rms^3 / (d E), with
 * H_rms = sqrt(8 E_m), where alpha is the proportionality constant (SDBC1),
 * f_m the mean frequency, d the water depth, E_m the mean wave energy (EMEAN)
 * and E the total wave energy obtained by integrating the spectrum. D_db is
 * zero if E does not exceed the single precision threshold TINY(1.0).
 *
 * Sign of the source term: as in W3SDB1 (IWB = 2), D_db is positive, so that
 * this source term adds energy instead of dissipating it, and the breaking
 * flag (LBREAK) is always false. This is reproduced exactly here, and needs to
 * be reviewed before DB2 is used.
 *
 * Input and output spectrum: the source term is a wave action density source
 * term S_db(k,theta) on the same (k,theta) spectral grid as the action density
 * spectrum N(k,theta) = F(k,theta) / sigma used internally by WW3, where
 * F(k,theta) is the variance (energy) density spectrum and sigma the radian
 * (intrinsic) frequency. The total energy is
 * E = sum_k [ (sum_theta N(k,theta)) DDEN(k) / c_g(k) ], with
 * DDEN = d_sigma d_theta sigma and c_g the group velocity. The spectra are
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
 * @author Code Heritage: Converted and adapted from the IWB = 2 option of
 * WAVEWATCH III (WW3) w3sdb1md.F90 (W3SDB1) by J. H. G. M. Alves and Hendrik
 * L. Tolman, with contributions by A. Roland.
 * @date Initial: 2026-11-20
 * @date Last update: 2026-11-20
 */

#ifndef WW4_DB2_WW3_H
#define WW4_DB2_WW3_H

// System
#include <cstddef>
#include <optional>
#include <span>

namespace ww4_db2 {

// --- compute_w3sdb2 ---------------------------------------------------------
/**
 * @brief Computes the DB2 Thornton-Guza type depth-induced breaking source
 * term (W3SDB1 option IWB = 2).
 * @details Source term S_db = D_db N with
 * D_db = alpha (sqrt(pi) / 16) f_m H_rms^3 / (d E), see file description for
 * the sign. If the mean wave energy is zero, nothing is computed and the output
 * arrays are not changed.
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-20
 * @date Last update: 2026-11-20
 */
[[nodiscard]] std::optional<bool>
compute_w3sdb2(std::span<const float> action_density,
               std::span<const float> group_velocity, float depth,
               float mean_energy, float mean_frequency,
               std::span<float> source_term, std::span<float> diagonal_term,
               std::size_t num_directions,
               std::span<const float> integration_factor,
               float proportionality_constant) noexcept;

} // namespace ww4_db2

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sdb2_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_db2::compute_w3sdb2.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_db2::compute_w3sdb2. For binding as BIND(C, name="w3sdb2_cpp").
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-20
 * @date Last update: 2026-11-20
 */
int w3sdb2_cpp(const float *A, const float *CG, float depth, float emean,
               float fmean, float *S, float *D, int nth, int nk,
               const float *dden, float sdbc1) noexcept;

} // extern "C"

#endif // WW4_DB2_WW3_H
