/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file LN1_ww3.h
 * @brief Declarations of the single-point linear wind input (LN1) source term
 * routine for WAVEWATCH IV (WW4).
 * @details Declares the C++20 LN1 routine for the linear wind input source
 * term (W3SLN1) for one spatial point, and its C-compatible interface called
 * from the WAVEWATCH III (WW3) Fortran module w3sln1md.F90. Implementations
 * are in LN1_ww3.cpp. The ww4_ln1 routine uses descriptive variable names. The
 * extern "C" wrapper keeps the WW3 Fortran names (K, FHIGH, USTAR, USDIR, ...)
 * for easy cross-referencing with w3sln1md.F90.
 *
 * Linear wind input according to Cavaleri and Malanotte-Rizzoli (1981),
 * filtered for low frequencies according to Tolman (1992). Converted to
 * action density as a function of wavenumber k and direction theta, the
 * source term is
 * S_ln(k,theta) = C_ln k^-1 max[0, u_* cos(theta - theta_w)]^4 G,
 * with C_ln = 80 (rho_a / rho_w)^2 g^-2 (SLNC1), and low-frequency filter
 * G = exp[-(sigma / sigma_filt)^-4], which is set to zero for
 * sigma / sigma_filt < 0.5. The filter frequency is
 * sigma_filt = min[max(alpha_PM g / (28 u_*), alpha_hf min(sigma_nk,
 * sigma_hf)), 2 sigma_nk], where sigma_nk is the highest discrete radian
 * frequency and sigma_hf the cut-off frequency of the prognostic part of the
 * spectrum.
 *
 * Output spectrum: the source term is a wave action density source term
 * S_ln(k,theta) on the same (k,theta) spectral grid as the action density
 * spectrum N(k,theta) = F(k,theta) / sigma used internally by WW3, where
 * F(k,theta) is the variance (energy) density spectrum and sigma the radian
 * (intrinsic) frequency. The 2-D spectrum is stored in Fortran (column-major)
 * order S(ith, ik), so that the 1-D (flattened) index of direction ith and
 * wavenumber ik (both 0-based) is is = ith + nth * ik.
 *
 * References:
 * - Cavaleri, L. and P. Malanotte-Rizzoli, 1981: Wind-wave prediction in
 *   shallow water: Theory and applications. J. Geophys. Res., 86,
 *   10961-10973.
 * - Tolman, H. L., 1992: Effects of numerics on the physics in a
 *   third-generation wind-wave model. J. Phys. Oceanogr., 22, 1095-1111.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section LN1.
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
 * w3sln1md.F90 (W3SLN1) by Hendrik L. Tolman.
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */

#ifndef WW4_LN1_WW3_H
#define WW4_LN1_WW3_H

// System
#include <cstddef>
#include <span>

namespace ww4_ln1 {

// --- compute_w3sln1 ---------------------------------------------------------
/**
 * @brief Computes the LN1 linear wind input source term (W3SLN1).
 * @details Linear wind input of Cavaleri and Malanotte-Rizzoli (1981) filtered
 * for low frequencies as in Tolman (1992), see file description. The source
 * term is the product of a directional factor max[0, cos(theta - theta_w)]^4
 * and a frequency dependent factor C_ln u_*^4 k^-1 G. The order of operations
 * matches the Fortran original so results are reproducible.
 * @param[in] wavenumber Wavenumber k of each frequency band (size nk).
 * @param[in] cutoff_frequency Cut-off frequency sigma_hf in the integration
 * (rad/s) (FHIGH).
 * @param[in] friction_velocity Wind friction velocity u_*.
 * @param[in] wind_direction Direction theta_w of the friction velocity
 * (radians, USDIR).
 * @param[out] source_term Linear input source term in action density form as
 * 2-D array S(ith, ik) in Fortran (column-major) order (size nth * nk), see
 * file description.
 * @param[in] num_directions Number of directions nth.
 * @param[in] cos_direction Cosine of each spectral direction (size nth).
 * @param[in] sin_direction Sine of each spectral direction (size nth).
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] linear_input_coefficient Combined input constant
 * 80 (rho_a / rho_w)^2 g^-2 (SLNC1).
 * @param[in] pm_frequency_factor Factor alpha_PM on the Pierson-Moskowitz
 * frequency in the filter frequency (FSPM).
 * @param[in] high_frequency_factor Factor alpha_hf on the cut-off frequency in
 * the filter frequency (FSHF).
 * @param[in] gravity Acceleration of gravity g.
 * @pre wavenumber and radian_frequency have the same size nk > 0,
 * cos_direction and sin_direction have size nth = num_directions, and
 * source_term has size nth * nk.
 * @post The array source_term is populated with the calculated linear input
 * source term.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SLN1 by Hendrik L. Tolman.
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */
void compute_w3sln1(std::span<const float> wavenumber, float cutoff_frequency,
                    float friction_velocity, float wind_direction,
                    std::span<float> source_term, std::size_t num_directions,
                    std::span<const float> cos_direction,
                    std::span<const float> sin_direction,
                    std::span<const float> radian_frequency,
                    float linear_input_coefficient, float pm_frequency_factor,
                    float high_frequency_factor, float gravity) noexcept;

} // namespace ww4_ln1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sln1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_ln1::compute_w3sln1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_ln1::compute_w3sln1. S is the Fortran array S(nth, nk). Bound in
 * w3sln1md.F90 as BIND(C, name="w3sln1_cpp").
 * @param[in] K Wavenumbers (size nk).
 * @param[in] fhigh Cut-off frequency in the integration (rad/s).
 * @param[in] ustar Wind friction velocity u_*.
 * @param[in] usdir Direction of the friction velocity (radians).
 * @param[out] S Linear input source term as 2-D array S(nth, nk) (size
 * nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] ecos Cosine of each spectral direction (size nth).
 * @param[in] esin Sine of each spectral direction (size nth).
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] slnc1 Combined input constant (SLNC1).
 * @param[in] fspm Factor on the Pierson-Moskowitz frequency in the filter.
 * @param[in] fshf Factor on the cut-off frequency in the filter.
 * @param[in] grav Acceleration of gravity.
 * @pre nk > 0 and all arrays are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */
void w3sln1_cpp(const float *K, float fhigh, float ustar, float usdir, float *S,
                int nth, int nk, const float *ecos, const float *esin,
                const float *sig, float slnc1, float fspm, float fshf,
                float grav) noexcept;

} // extern "C"

#endif // WW4_LN1_WW3_H
