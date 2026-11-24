/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file TR1_ww3.h
 * @brief Declarations of the single-point Lumped Triad Approximation (TR1)
 * triad interaction source term routine for WAVEWATCH IV (WW4).
 * @details Declares the C++20 TR1 routine for the triad nonlinear interaction
 * source term (W3STR1) for one spatial point, and its C-compatible interface
 * called from the WAVEWATCH III (WW3) Fortran module w3str1md.F90.
 * Implementations are in TR1_ww3.cpp. The ww4_tr1 routine uses descriptive
 * variable names. The extern "C" wrapper keeps the WW3 Fortran names (A, CG,
 * WN, DEPTH, S, D, ...) for easy cross-referencing with w3str1md.F90.
 *
 * Triad interactions according to the Lumped Triad Approximation (LTA) of
 * Eldeberky (1996), as implemented in SWAN. Only self-sum interactions are
 * considered, in which a component of frequency sigma / 2 interacts with
 * itself to exchange energy with the component at sigma. The interactions are
 * computed per direction (directionally decoupled). For the energy density
 * E(sigma,theta), the source term is
 * S_tr(sigma,theta) = S+(sigma,theta) - 2 S+(2 sigma,theta), with
 * S+(sigma) = max[0, alpha c c_g J^2 |sin beta| (E(sigma/2)^2
 * - 2 E(sigma/2) E(sigma))], where alpha is a proportionality constant, c and
 * c_g the phase and group velocity at sigma, and J the interaction
 * coefficient of Madsen and Sorensen (1993),
 * J = k'^2 (g d + 2 c'^2) / (k d (g d + (2/15) g d^3 k^2 - (2/5) sigma^2 d^2)),
 * with primed quantities at sigma / 2. The biphase is
 * beta = (pi / 2) (tanh(0.2 / Ur) - 1), with the Ursell number
 * Ur = g H_s / (2 sqrt(2) sigma_m01^2 d^2), H_s = 4 sqrt(E_tot) and sigma_m01
 * the mean radian frequency from the first spectral moment. As the frequency
 * grid is logarithmic, the components at sigma / 2 and 2 sigma are obtained by
 * linear interpolation between the two nearest discrete frequencies.
 *
 * Interactions are computed only if Ur >= 0.01, and only for frequencies below
 * 10 sigma_m01 (and at least up to the band at 2 sigma_1). The constants are
 * fixed in the WW3 code: alpha = 1, frequency limit 10 sigma_m01, biphase
 * constant 0.2 and minimum Ursell number 0.01. Note that the WW3 manual (TR1)
 * gives alpha = 0.05 and a frequency limit of 2.5 f_m01.
 *
 * Input and output spectrum: the input is the wave action density spectrum
 * N(k,theta) = F(k,theta) / sigma used internally by WW3, where F(k,theta) is
 * the variance (energy) density spectrum and sigma the radian (intrinsic)
 * frequency. The interactions are computed for the energy density
 * E(sigma,theta) = 2 pi sigma N(k,theta) / c_g, and converted back to an
 * action density source term S(k,theta) on the same (k,theta) grid, with
 * diagonal term D = S / N (zero where N <= 1e-12). The spectra are stored as
 * 1-D (flattened) arrays in Fortran order, so that the index of direction ith
 * and wavenumber ik (both 0-based) is is = ith + nth * ik, with nth the number
 * of directions. As in WW3, the source term and diagonal term are not set for
 * the highest wavenumber band, and are not set at all if the spectrum is zero.
 * The caller initializes them (to zero in W3SRCE).
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
 * w3str1md.F90 (W3STR1) by A. J. van der Westhuysen, with contributions by A.
 * Roland, based on the SWAN code by M. Zijlema.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */

#ifndef WW4_TR1_WW3_H
#define WW4_TR1_WW3_H

// System
#include <cstddef>
#include <span>

namespace ww4_tr1 {

// --- compute_w3str1 ---------------------------------------------------------
/**
 * @brief Computes the TR1 Lumped Triad Approximation source term (W3STR1).
 * @details Self-sum triad interactions of Eldeberky (1996) per direction, see
 * file description. If the spectrum is zero, nothing is computed and the
 * output arrays are not changed. The highest wavenumber band of the outputs is
 * never set.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nspec = nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] wavenumber Wavenumber k of each frequency band (size nk).
 * @param[in] depth Mean water depth d.
 * @param[in,out] source_term Triad source term S in action density form as
 * 1-D array (size nspec).
 * @param[in,out] diagonal_term Diagonal term D of the source term derivative
 * as 1-D array (size nspec).
 * @param[in] num_directions Number of directions nth.
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] integration_factor Factor DDEN = d_sigma d_theta sigma of each
 * frequency band for the spectral integrals (size nk).
 * @param[in] tail_energy_factor Tail factor for the energy integral (FTE).
 * @param[in] tail_frequency_factor Tail factor for the first moment (FTF).
 * @param[in] gravity Acceleration of gravity g.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void compute_w3str1(std::span<const float> action_density,
                    std::span<const float> group_velocity,
                    std::span<const float> wavenumber, float depth,
                    std::span<float> source_term,
                    std::span<float> diagonal_term, std::size_t num_directions,
                    std::span<const float> radian_frequency,
                    std::span<const float> integration_factor,
                    float tail_energy_factor, float tail_frequency_factor,
                    float gravity) noexcept;

} // namespace ww4_tr1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3str1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_tr1::compute_w3str1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_tr1::compute_w3str1. Bound in w3str1md.F90 as
 * BIND(C, name="w3str1_cpp").
 * @param[in] A Action density spectrum as 1-D array (size nth * nk).
 * @param[in] CG Group velocities (size nk).
 * @param[in] WN Wavenumbers (size nk).
 * @param[in] depth Mean water depth.
 * @param[in,out] S Triad source term as 1-D array (size nth * nk).
 * @param[in,out] D Diagonal term of the derivative as 1-D array (size
 * nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Integration factor DDEN (size nk).
 * @param[in] fte Tail factor for the energy integral (FTE).
 * @param[in] ftf Tail factor for the first moment (FTF).
 * @param[in] grav Acceleration of gravity.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void w3str1_cpp(const float *A, const float *CG, const float *WN, float depth,
                float *S, float *D, int nth, int nk, const float *sig,
                const float *dden, float fte, float ftf, float grav) noexcept;

} // extern "C"

#endif // WW4_TR1_WW3_H
