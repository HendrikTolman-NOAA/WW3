/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file BT4_ww3.h
 * @brief Declarations of the single-point SHOWEX movable-bed bottom friction
 * (BT4) source term routines for WAVEWATCH IV (WW4).
 * @details Declares the C++20 BT4 routines for the SHOWEX bottom friction
 * source term (W3SBT4) for one spatial point and for its initialization
 * (TABU_ERF, called from INSBT4), and their C-compatible interfaces called
 * from the WAVEWATCH III (WW3) Fortran module w3sbt4md.F90. Implementations
 * are in BT4_ww3.cpp. The ww4_bt4 routines use descriptive variable names. The
 * extern "C" wrappers keep the WW3 Fortran names (A, CG, WN, DEPTH, D50, PSIC,
 * TAUBBL, BEDFORM, S, D, ...) for easy cross-referencing with w3sbt4md.F90.
 *
 * Bottom friction over a movable sandy bed, based on the eddy viscosity model
 * of Grant and Madsen (1979) with the roughness parameterization of Tolman
 * (1994) as adjusted by Ardhuin et al. (2003), and a sub-grid treatment of
 * the water depth following Tolman (1995). The source term is linear in the
 * spectrum,
 * S_bt(k,theta) = -f_e u_b sigma^2 / (2 g sinh^2(k d)) N(k,theta),
 * where f_e is the dissipation factor (from a table of the friction factor as
 * a function of a_b / k_N), u_b the r.m.s. bottom orbital velocity, a_b the
 * r.m.s. bottom orbital displacement and k_N the Nikuradse roughness. The
 * roughness is the sum of a ripple roughness k_r = a_b A_1 (psi / psi_c)^A_2
 * and a sheet flow roughness k_s = 0.0655 a_b (u_b^2 / ((s - 1) g a_b))^1.4
 * (Wilson, 1989) when the normalized Shields number psi / psi_c exceeds A_3,
 * and a relict roughness max(A_5, A_6 D_50, A_4 a_b) otherwise. Here psi is
 * the Shields number from the skin friction factor (roughness D_50), psi_c
 * the critical Shields number, s the sediment specific gravity and A_1 to A_6
 * the empirical coefficients of the BT4 namelist (SBTCX).
 *
 * Sub-grid depth variability: the Shields number is computed for the depths
 * d (1 - 0.001), d (1 + 0.001) and d. If its variation with depth, scaled
 * with the relative depth standard deviation sigma_d (SBTCX(5)), is large
 * enough, the roughness is the probability-weighted sum of the relict and
 * the active (ripple and sheet flow) roughness, with the probability of
 * active ripples taken from a tabulated function (erf table). The bulk
 * quantities (u_b, a_b, ripple wavelength) are those of the central depth d.
 *
 * Other outputs: the bedform parameters (roughness k_N and the two components
 * of the potential ripple wavelength, 0.6 times the orbital diameter along
 * the mean direction of the orbital motion, negative for relict ripples), and
 * the momentum flux to the bottom boundary layer TAUBBL.
 *
 * Input and output spectrum: the source term is a wave action density source
 * term S_bt(k,theta) on the same (k,theta) spectral grid as the action density
 * spectrum N(k,theta) = F(k,theta) / sigma used internally by WW3, where
 * F(k,theta) is the variance (energy) density spectrum and sigma the radian
 * (intrinsic) frequency. The diagonal term depends on wavenumber only. The
 * spectra are stored as 1-D (flattened) arrays in Fortran order, so that the
 * index of direction ith and wavenumber ik (both 0-based) is is = ith + nth *
 * ik, with nth the number of directions.
 *
 * Possible bugs in the WW3 code, reproduced here to keep the results
 * identical (to be reviewed before use in WW4):
 * - The momentum flux TAUBBL uses ECOS(IS) and ESIN(IS) with the spectral index
 *   IS instead of the direction index ITH. In WW3 only the first NTH elements
 *   of ECOS and ESIN (the directions) are set and the others are zero, so that
 *   only the first frequency band contributes to TAUBBL. The arrays are passed
 *   in full for this reason. The intended expression is probably the sum over
 *   all bands with ECOS(ITH) and ESIN(ITH).
 * - The erf table is (1 + erf(x)) / 2 for x from -4 to 4. The factor
 *   2**(1/2) in TABU_ERF is an integer expression equal to 1 (1/2 is integer
 *   division). The cumulative normal distribution (1 + erf(x / sqrt(2))) / 2,
 *   with x in standard deviations (XERFMAX = 4 "number of stdev"), may have
 *   been intended.
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
 * w3sbt4md.F90 (W3SBT4, INSBT4 and TABU_ERF) by F. Ardhuin and J. Lepesqueur.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */

#ifndef WW4_BT4_WW3_H
#define WW4_BT4_WW3_H

// System
#include <cstddef>
#include <span>

namespace ww4_bt4 {

// --- compute_tabu_erf -------------------------------------------------------
/**
 * @brief Tabulates the function (1 + erf(x)) / 2 used in the sub-grid
 * treatment (TABU_ERF).
 * @details The table has equally spaced values of x from -x_max to x_max,
 * see file description.
 * @param[out] erf_table Tabulated values (size n + 1 for n intervals,
 * ERFTABLE(0:SIZEERFTABLE)).
 * @param[in] max_stdev Range x_max of the table (XERFMAX).
 * @return The spacing of the table in x (DELXERF).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
[[nodiscard]] float compute_tabu_erf(std::span<float> erf_table,
                                     float max_stdev) noexcept;

// --- compute_w3sbt4 ---------------------------------------------------------
/**
 * @brief Computes the BT4 SHOWEX bottom friction source term (W3SBT4).
 * @details Movable-bed bottom friction with ripple, sheet flow and relict
 * roughness and a sub-grid treatment of the depth, see file description.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nspec = nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] wavenumber Wavenumber k of each frequency band (size nk).
 * @param[in] depth Mean water depth d.
 * @param[in] grain_size Median grain size D_50 (D50).
 * @param[in] critical_shields Critical Shields number psi_c (PSIC).
 * @param[out] bottom_stress Momentum flux to the bottom boundary layer
 * (TAUBBL, size 2).
 * @param[out] bedform Bedform roughness and ripple wavelength components
 * (BEDFORM, size 3).
 * @param[out] source_term Bottom friction source term S_bt as 1-D array (size
 * nspec).
 * @param[out] diagonal_term Diagonal term of the source term derivative as 1-D
 * array (size nspec).
 * @param[in] num_directions Number of directions nth.
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] integration_factor Factor DDEN = d_sigma d_theta sigma of each
 * frequency band (size nk).
 * @param[in] cos_table WW3 array ECOS, cosines of the directions in the first
 * nth elements (size nspec), see file description.
 * @param[in] sin_table WW3 array ESIN, sines of the directions in the first nth
 * elements (size nspec).
 * @param[in] parameters Empirical coefficients SBTCX(1:7): A_1, A_2, A_3, A_4,
 * sigma_d, A_5, A_6 (size 7).
 * @param[in] friction_table Friction factor table FWTABLE(0:SIZEFWTABLE).
 * @param[in] friction_log_min Lower limit of log10(a_b / k_N) in the friction
 * factor table (ABMIN).
 * @param[in] friction_log_step Step of log10(a_b / k_N) in the friction factor
 * table (DELAB).
 * @param[in] erf_table Table of (1 + erf(x)) / 2 (ERFTABLE(0:SIZEERFTABLE)).
 * @param[in] erf_max Range x_max of the erf table (XERFMAX).
 * @param[in] erf_step Spacing of the erf table (DELXERF).
 * @param[in] gravity Acceleration of gravity g.
 * @param[in] sediment_gravity Specific gravity s of the sediment (SED_SG).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void compute_w3sbt4(
    std::span<const float> action_density,
    std::span<const float> group_velocity, std::span<const float> wavenumber,
    float depth, float grain_size, float critical_shields,
    std::span<float> bottom_stress, std::span<float> bedform,
    std::span<float> source_term, std::span<float> diagonal_term,
    std::size_t num_directions, std::span<const float> radian_frequency,
    std::span<const float> integration_factor, std::span<const float> cos_table,
    std::span<const float> sin_table, std::span<const float> parameters,
    std::span<const float> friction_table, float friction_log_min,
    float friction_log_step, std::span<const float> erf_table, float erf_max,
    float erf_step, float gravity, float sediment_gravity) noexcept;

} // namespace ww4_bt4

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- tabu_erf_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to ww4_bt4::compute_tabu_erf.
 * @details Bound in w3sbt4md.F90 as BIND(C, name="tabu_erf_cpp").
 * @param[out] erftable Tabulated values ERFTABLE(0:sizeerftable) (size
 * sizeerftable + 1).
 * @param[in] sizeerftable Number of intervals of the table (SIZEERFTABLE).
 * @param[in] xerfmax Range of the table (XERFMAX).
 * @return The spacing of the table (DELXERF).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
float tabu_erf_cpp(float *erftable, int sizeerftable, float xerfmax) noexcept;

// --- w3sbt4_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_bt4::compute_w3sbt4.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_bt4::compute_w3sbt4. Bound in w3sbt4md.F90 as
 * BIND(C, name="w3sbt4_cpp").
 * @param[in] A Action density spectrum as 1-D array (size nth * nk).
 * @param[in] CG Group velocities (size nk).
 * @param[in] WN Wavenumbers (size nk).
 * @param[in] depth Water depth.
 * @param[in] d50 Median grain size.
 * @param[in] psic Critical Shields number.
 * @param[out] TAUBBL Momentum flux to the bottom boundary layer (size 2).
 * @param[out] BEDFORM Bedform parameters (size 3).
 * @param[out] S Bottom friction source term (size nth * nk).
 * @param[out] D Diagonal term of the derivative (size nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Integration factor DDEN (size nk).
 * @param[in] ecos Array ECOS (size nth * nk).
 * @param[in] esin Array ESIN (size nth * nk).
 * @param[in] sbtcx Coefficients SBTCX(1:7) (size 7).
 * @param[in] fwtable Friction factor table (size sizefwtable + 1).
 * @param[in] sizefwtable Number of intervals of the friction factor table.
 * @param[in] abmin Lower limit of the friction factor table (ABMIN).
 * @param[in] delab Step of the friction factor table (DELAB).
 * @param[in] erftable Erf table (size sizeerftable + 1).
 * @param[in] sizeerftable Number of intervals of the erf table.
 * @param[in] xerfmax Range of the erf table (XERFMAX).
 * @param[in] delxerf Spacing of the erf table (DELXERF).
 * @param[in] grav Acceleration of gravity.
 * @param[in] sed_sg Specific gravity of the sediment (SED_SG).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void w3sbt4_cpp(const float *A, const float *CG, const float *WN, float depth,
                float d50, float psic, float *TAUBBL, float *BEDFORM, float *S,
                float *D, int nth, int nk, const float *sig, const float *dden,
                const float *ecos, const float *esin, const float *sbtcx,
                const float *fwtable, int sizefwtable, float abmin, float delab,
                const float *erftable, int sizeerftable, float xerfmax,
                float delxerf, float grav, float sed_sg) noexcept;

} // extern "C"

#endif // WW4_BT4_WW3_H
