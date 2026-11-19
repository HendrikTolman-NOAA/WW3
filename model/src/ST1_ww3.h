/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file ST1_ww3.h
 * @brief Declarations of the single-point WAM-3 (ST1) source term routines for
 * WAVEWATCH IV (WW4).
 * @details Declares the C++20 ST1 routines for the mean wave parameters
 * (W3SPR1), wind input (W3SIN1) and whitecapping dissipation (W3SDS1) for one
 * spatial point, and their C-compatible interfaces called from the WAVEWATCH
 * III (WW3) Fortran module w3src1md.F90. Implementations are in ST1_ww3.cpp.
 * The ww4_st1 routines use descriptive variable names. The extern "C"
 * wrappers keep the WW3 Fortran names (A, K, SIG2, EMEAN, ...) for easy
 * cross-referencing with w3src1md.F90.
 *
 * Input spectrum: all routines take the wave action density spectrum
 * N(k,theta) = F(k,theta) / sigma, as a function of wavenumber k and wave
 * direction theta, where F(k,theta) is the variance (energy) density spectrum
 * and sigma the radian (intrinsic) frequency. This is the spectrum used
 * internally by WW3. The 2-D spectrum is stored in Fortran (column-major)
 * order A(ith, ik), so that the 1-D (flattened) index of direction ith and
 * wavenumber ik (both 0-based) is is = ith + nth * ik.
 *
 * References:
 * - Snyder, R. L., F. W. Dobson, J. A. Elliott and R. B. Long, 1981: Array
 *   measurements of atmospheric pressure fluctuations above surface gravity
 *   waves. J. Fluid Mech., 102, 1-59.
 * - Komen, G. J., S. Hasselmann and K. Hasselmann, 1984: On the existence of
 *   a fully developed wind-sea spectrum. J. Phys. Oceanogr., 14, 1271-1285.
 * - The WAMDI Group, 1988: The WAM model - A third generation ocean wave
 *   prediction model. J. Phys. Oceanogr., 18, 1775-1810.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section ST1.
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
 * @date Initial: 2026-11-18
 * @date Last update: 2026-11-19
 */

#ifndef WW4_ST1_WW3_H
#define WW4_ST1_WW3_H

// System
#include <cstddef>
#include <span>

namespace ww4_st1 {

// --- compute_w3sin1 ---------------------------------------------------------
/**
 * @brief Computes the ST1 wind input source term (W3SIN1).
 * @details WAM-3 wind input of Snyder et al. (1981) as adapted by Komen et al.
 * (1984), S_in = D * N with diagonal term D.
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
void compute_w3sin1(std::span<const float> action_density,
                    std::span<const float> wavenumber, float friction_velocity,
                    float cos_wind_dir, float sin_wind_dir,
                    std::span<float> source_term,
                    std::span<float> source_diagonal, float input_coefficient,
                    std::span<const float> radian_frequency,
                    std::span<const float> cos_direction,
                    std::span<const float> sin_direction) noexcept;

// --- compute_w3sds1 ---------------------------------------------------------
/**
 * @brief Computes the ST1 whitecapping dissipation source term (W3SDS1).
 * @details WAM-3 whitecapping of Komen et al. (1984), S_ds = D * N with
 * diagonal term D.
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
void compute_w3sds1(std::span<const float> action_density,
                    std::span<const float> wavenumber, float total_variance,
                    float mean_frequency, float mean_wavenumber,
                    std::span<float> source_term,
                    std::span<float> source_diagonal,
                    float dissipation_coefficient) noexcept;

// --- compute_w3spr1 ---------------------------------------------------------
/**
 * @brief Computes the mean wave parameters used by the ST1 source terms
 * (W3SPR1).
 * @details Calculates the total variance, mean frequency, mean wavenumber and
 * maximum action density of the WAM-3 method (WAMDI Group, 1988).
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-17
 * @date Last update: 2026-11-19
 */
void compute_w3spr1(std::span<const float> action_density,
                    std::span<const float> group_velocity,
                    std::span<const float> wavenumber,
                    std::span<const float> radian_frequency,
                    std::span<const float> band_integration_factor,
                    std::size_t num_directions, float tail_factor_energy,
                    float tail_factor_frequency, float tail_factor_wavenumber,
                    float inv_two_pi, float &total_variance,
                    float &mean_frequency, float &mean_wavenumber,
                    float &max_action_density) noexcept;

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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
void w3sin1_cpp(const float *A, const float *K, float ustar, float cosu,
                float sinu, float *S, float *D, int nspec, float sinc1,
                const float *sig2, const float *ecos,
                const float *esin) noexcept;

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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-19
 */
void w3sds1_cpp(const float *A, const float *K, float emean, float fmean,
                float wnmean, float *S, float *D, int nspec,
                float sdsc1) noexcept;

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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-17
 * @date Last update: 2026-11-19
 */
void w3spr1_cpp(const float *A, const float *cg, const float *wn,
                const float *sig, const float *dden, int nth, int nk, float fte,
                float ftf, float ftwn, float tpiinv, float *emean, float *fmean,
                float *wnmean, float *amax) noexcept;

} // extern "C"

#endif // WW4_ST1_WW3_H
