/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file NL1_ww3.h
 * @brief Declarations of the single-point Discrete Interaction Approximation
 * (NL1) nonlinear interaction source term routine for WAVEWATCH IV (WW4).
 * @details Declares the C++20 NL1 routine for the nonlinear
 * wave-wave interactions (W3SNL1) for one spatial point, and its C-compatible
 * interface called from the WAVEWATCH III (WW3) Fortran module w3snl1md.F90.
 * Implementations are in NL1_ww3.cpp. The ww4_nl1 routine uses descriptive
 * variable names. The extern "C" wrapper keeps the WW3 Fortran names (A, CG,
 * KDMEAN, IP11, AWG1, ...) for easy cross-referencing with w3snl1md.F90.
 *
 * Discrete Interaction Approximation (DIA): for each spectral component k
 * only two mirror-image quadruplets are used, with k_1 = k, f_r2 = (1 +
 * lambda) f_r and f_r3 = (1 - lambda) f_r (lambda = 0.25 by default), and the
 * angles given by the geometry of the "figure of 8" (Hasselmann et al.,
 * 1985). The source term is calculated for the variance density spectrum
 * F(f_r,theta), with f_r the relative frequency, using the elementary
 * contributions
 * dS(l,m,n) = C / g^4 f_r,l^11 [F_l^2 (F_m / (1 + lambda)^4 +
 * F_n / (1 - lambda)^4) - 2 F_l F_m F_n / (1 - lambda^2)^4]
 * for the quadruplets in which k takes each of the roles. Spectral densities
 * at the interacting components are obtained with bilinear interpolation in
 * frequency and direction. For shallow water the source term is scaled by
 * D = 1 + c_1 / (k d) [1 - c_2 k d] exp(-c_3 k d), with the mean relative
 * depth k d = max(KDCON k_hat d, KDMN) (Hasselmann et al., 1985).
 *
 * Input and output spectrum: the routine takes the wave action density
 * spectrum N(k,theta) = F(k,theta) / sigma, as a function of wavenumber k and
 * wave direction theta, where F(k,theta) is the variance (energy) density
 * spectrum and sigma the radian (intrinsic) frequency. This is the spectrum
 * used internally by WW3. Internally it is converted to F(f_r,theta) =
 * 2 pi sigma N(k,theta) / c_g, and extended with a parametric tail to higher
 * frequencies and with empty bins at low frequencies. The source term and its
 * diagonal are returned in action density form on the original grid. The 2-D
 * spectrum is stored in Fortran (column-major) order A(ith, ik), so that the
 * 1-D (flattened) index of direction ith and wavenumber ik (both 0-based) is
 * is = ith + nth * ik. The address arrays (IP11, ..., IC82) hold the 1-based
 * Fortran addresses set in INSNL1. In the extended arrays, the addresses
 * 1 - nth to 0 are the empty bins at low frequencies.
 *
 * References:
 * - Hasselmann, S. and K. Hasselmann, 1985: Computations and
 *   parameterizations of the nonlinear energy transfer in a gravity-wave
 *   spectrum. Part I: A new method for efficient computations of the exact
 *   nonlinear transfer integral. J. Phys. Oceanogr., 15, 1369-1377.
 * - Hasselmann, S., K. Hasselmann, J. H. Allender and T. P. Barnett, 1985:
 *   Computations and parameterizations of the nonlinear energy transfer in a
 *   gravity-wave spectrum. Part II: Parameterizations of the nonlinear energy
 *   transfer for application in wave models. J. Phys. Oceanogr., 15,
 *   1378-1391.
 * - The WAMDI Group, 1988: The WAM model - A third generation ocean wave
 *   prediction model. J. Phys. Oceanogr., 18, 1775-1810.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section NL1.
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
 * w3snl1md.F90 (W3SNL1) by Hendrik L. Tolman.
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */

#ifndef WW4_NL1_WW3_H
#define WW4_NL1_WW3_H

// System
#include <cstddef>
#include <span>

namespace ww4_nl1 {

// --- compute_w3snl1 ---------------------------------------------------------
/**
 * @brief Computes the NL1 nonlinear wave-wave interaction source term (W3SNL1).
 * @details Discrete Interaction Approximation (DIA) of Hasselmann et al.
 * (1985), see file description. Calculates the nonlinear interactions and the
 * diagonal term of their derivative in four steps: the proportionality
 * constant including the shallow water factor, the auxiliary spectrum
 * F(f_r,theta) with tail and empty low frequency bins, the interactions for
 * the extended spectrum, and the source and diagonal terms obtained by
 * folding the interactions back to the discrete spectrum. The order of
 * operations matches the Fortran original so results are reproducible.
 * @param[in] action_density Wave action density spectrum N(k,theta) in 1-D
 * (flattened) form (size nspec = nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] mean_relative_depth Mean relative depth k_hat d (KDMEAN).
 * @param[out] source_term Nonlinear interaction source term in action density
 * form (size nspec).
 * @param[out] source_diagonal Diagonal term of the derivative of source_term
 * with respect to action_density (size nspec).
 * @param[in] num_directions Number of directions nth.
 * @param[in] num_extended_frequencies Number of frequency bands of the extended
 * auxiliary spectrum, including the tail (NFRHGH).
 * @param[in] num_auxiliary_bins Number of bins of the extended auxiliary
 * spectrum (NSPECY = NFRHGH * nth).
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] tail_factor Factor applied per frequency band to extend the
 * spectrum beyond nk (FACHFE).
 * @param[in] dia_coefficient Scaled proportionality constant C / g^4 of the DIA
 * (SNLC1).
 * @param[in] depth_coefficient_1 Constant c_1 in the shallow water factor
 * (SNLS1).
 * @param[in] depth_coefficient_2 Constant c_2 in the shallow water factor
 * (SNLS2).
 * @param[in] depth_coefficient_3 Constant c_3 in the shallow water factor
 * (SNLS3).
 * @param[in] depth_conversion_factor Conversion factor from k_hat d to the mean
 * relative depth k d (KDCON).
 * @param[in] min_relative_depth Minimum mean relative depth k d (KDMN).
 * @param[in] inv_two_pi Inverse of 2*Pi (TPIINV), passed from Fortran so that
 * results stay bit-identical to WW3.
 * @param[in] plus_address_q1_c1 Address in the auxiliary spectrum of
 * interpolation corner 1 of the component at (1+lambda) f_r for quadruplet 1,
 * weight interp_weight_1 (IP11, size nspecx).
 * @param[in] plus_address_q1_c2 Address in the auxiliary spectrum of
 * interpolation corner 2 of the component at (1+lambda) f_r for quadruplet 1,
 * weight interp_weight_2 (IP12, size nspecx).
 * @param[in] plus_address_q1_c3 Address in the auxiliary spectrum of
 * interpolation corner 3 of the component at (1+lambda) f_r for quadruplet 1,
 * weight interp_weight_3 (IP13, size nspecx).
 * @param[in] plus_address_q1_c4 Address in the auxiliary spectrum of
 * interpolation corner 4 of the component at (1+lambda) f_r for quadruplet 1,
 * weight interp_weight_4 (IP14, size nspecx).
 * @param[in] minus_address_q1_c1 Address in the auxiliary spectrum of
 * interpolation corner 1 of the component at (1-lambda) f_r for quadruplet 1,
 * weight interp_weight_5 (IM11, size nspecx).
 * @param[in] minus_address_q1_c2 Address in the auxiliary spectrum of
 * interpolation corner 2 of the component at (1-lambda) f_r for quadruplet 1,
 * weight interp_weight_6 (IM12, size nspecx).
 * @param[in] minus_address_q1_c3 Address in the auxiliary spectrum of
 * interpolation corner 3 of the component at (1-lambda) f_r for quadruplet 1,
 * weight interp_weight_7 (IM13, size nspecx).
 * @param[in] minus_address_q1_c4 Address in the auxiliary spectrum of
 * interpolation corner 4 of the component at (1-lambda) f_r for quadruplet 1,
 * weight interp_weight_8 (IM14, size nspecx).
 * @param[in] plus_address_q2_c1 Address in the auxiliary spectrum of
 * interpolation corner 1 of the component at (1+lambda) f_r for quadruplet 2,
 * weight interp_weight_1 (IP21, size nspecx).
 * @param[in] plus_address_q2_c2 Address in the auxiliary spectrum of
 * interpolation corner 2 of the component at (1+lambda) f_r for quadruplet 2,
 * weight interp_weight_2 (IP22, size nspecx).
 * @param[in] plus_address_q2_c3 Address in the auxiliary spectrum of
 * interpolation corner 3 of the component at (1+lambda) f_r for quadruplet 2,
 * weight interp_weight_3 (IP23, size nspecx).
 * @param[in] plus_address_q2_c4 Address in the auxiliary spectrum of
 * interpolation corner 4 of the component at (1+lambda) f_r for quadruplet 2,
 * weight interp_weight_4 (IP24, size nspecx).
 * @param[in] minus_address_q2_c1 Address in the auxiliary spectrum of
 * interpolation corner 1 of the component at (1-lambda) f_r for quadruplet 2,
 * weight interp_weight_5 (IM21, size nspecx).
 * @param[in] minus_address_q2_c2 Address in the auxiliary spectrum of
 * interpolation corner 2 of the component at (1-lambda) f_r for quadruplet 2,
 * weight interp_weight_6 (IM22, size nspecx).
 * @param[in] minus_address_q2_c3 Address in the auxiliary spectrum of
 * interpolation corner 3 of the component at (1-lambda) f_r for quadruplet 2,
 * weight interp_weight_7 (IM23, size nspecx).
 * @param[in] minus_address_q2_c4 Address in the auxiliary spectrum of
 * interpolation corner 4 of the component at (1-lambda) f_r for quadruplet 2,
 * weight interp_weight_8 (IM24, size nspecx).
 * @param[in] fold_address_c1_q1 Address in the interaction arrays of the
 * contribution of weight 1 from quadruplet 1 to each spectral bin (IC11, size
 * nspec).
 * @param[in] fold_address_c1_q2 Address in the interaction arrays of the
 * contribution of weight 1 from quadruplet 2 to each spectral bin (IC12, size
 * nspec).
 * @param[in] fold_address_c2_q1 Address in the interaction arrays of the
 * contribution of weight 2 from quadruplet 1 to each spectral bin (IC21, size
 * nspec).
 * @param[in] fold_address_c2_q2 Address in the interaction arrays of the
 * contribution of weight 2 from quadruplet 2 to each spectral bin (IC22, size
 * nspec).
 * @param[in] fold_address_c3_q1 Address in the interaction arrays of the
 * contribution of weight 3 from quadruplet 1 to each spectral bin (IC31, size
 * nspec).
 * @param[in] fold_address_c3_q2 Address in the interaction arrays of the
 * contribution of weight 3 from quadruplet 2 to each spectral bin (IC32, size
 * nspec).
 * @param[in] fold_address_c4_q1 Address in the interaction arrays of the
 * contribution of weight 4 from quadruplet 1 to each spectral bin (IC41, size
 * nspec).
 * @param[in] fold_address_c4_q2 Address in the interaction arrays of the
 * contribution of weight 4 from quadruplet 2 to each spectral bin (IC42, size
 * nspec).
 * @param[in] fold_address_c5_q1 Address in the interaction arrays of the
 * contribution of weight 5 from quadruplet 1 to each spectral bin (IC51, size
 * nspec).
 * @param[in] fold_address_c5_q2 Address in the interaction arrays of the
 * contribution of weight 5 from quadruplet 2 to each spectral bin (IC52, size
 * nspec).
 * @param[in] fold_address_c6_q1 Address in the interaction arrays of the
 * contribution of weight 6 from quadruplet 1 to each spectral bin (IC61, size
 * nspec).
 * @param[in] fold_address_c6_q2 Address in the interaction arrays of the
 * contribution of weight 6 from quadruplet 2 to each spectral bin (IC62, size
 * nspec).
 * @param[in] fold_address_c7_q1 Address in the interaction arrays of the
 * contribution of weight 7 from quadruplet 1 to each spectral bin (IC71, size
 * nspec).
 * @param[in] fold_address_c7_q2 Address in the interaction arrays of the
 * contribution of weight 7 from quadruplet 2 to each spectral bin (IC72, size
 * nspec).
 * @param[in] fold_address_c8_q1 Address in the interaction arrays of the
 * contribution of weight 8 from quadruplet 1 to each spectral bin (IC81, size
 * nspec).
 * @param[in] fold_address_c8_q2 Address in the interaction arrays of the
 * contribution of weight 8 from quadruplet 2 to each spectral bin (IC82, size
 * nspec).
 * @param[in] lambda_weight_1 Weight 1 / (1 + lambda)^4 (DAL1).
 * @param[in] lambda_weight_2 Weight 1 / (1 - lambda)^4 (DAL2).
 * @param[in] lambda_weight_3 Weight 2 * lambda_weight_1 * lambda_weight_2
 * (DAL3).
 * @param[in] frequency_scaling Scaling array f_r^11 of the extended spectrum
 * (AF11, size nspecx).
 * @param[in] interp_weight_1 Bilinear interpolation weight 1 (AWG1).
 * @param[in] interp_weight_2 Bilinear interpolation weight 2 (AWG2).
 * @param[in] interp_weight_3 Bilinear interpolation weight 3 (AWG3).
 * @param[in] interp_weight_4 Bilinear interpolation weight 4 (AWG4).
 * @param[in] interp_weight_5 Bilinear interpolation weight 5 (AWG5).
 * @param[in] interp_weight_6 Bilinear interpolation weight 6 (AWG6).
 * @param[in] interp_weight_7 Bilinear interpolation weight 7 (AWG7).
 * @param[in] interp_weight_8 Bilinear interpolation weight 8 (AWG8).
 * @param[in] squared_weight_1 Squared interpolation weight 1 for the diagonal
 * term (SWG1).
 * @param[in] squared_weight_2 Squared interpolation weight 2 for the diagonal
 * term (SWG2).
 * @param[in] squared_weight_3 Squared interpolation weight 3 for the diagonal
 * term (SWG3).
 * @param[in] squared_weight_4 Squared interpolation weight 4 for the diagonal
 * term (SWG4).
 * @param[in] squared_weight_5 Squared interpolation weight 5 for the diagonal
 * term (SWG5).
 * @param[in] squared_weight_6 Squared interpolation weight 6 for the diagonal
 * term (SWG6).
 * @param[in] squared_weight_7 Squared interpolation weight 7 for the diagonal
 * term (SWG7).
 * @param[in] squared_weight_8 Squared interpolation weight 8 for the diagonal
 * term (SWG8).
 * @pre action_density, source_term and source_diagonal have size nspec =
 * nth * nk, group_velocity and radian_frequency have size nk > 0, the plus
 * and minus address arrays and frequency_scaling have size nspecx, the fold
 * address arrays have size nspec, and num_auxiliary_bins =
 * num_extended_frequencies * nth. All addresses are valid for the auxiliary
 * arrays (1 - nth to num_auxiliary_bins, and 1 - nth to nspecx for the
 * interaction arrays).
 * @post The arrays source_term and source_diagonal are populated with the
 * calculated nonlinear interaction source term and diagonal term.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SNL1 by Hendrik L. Tolman.
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */
void compute_w3snl1(
    std::span<const float> action_density,
    std::span<const float> group_velocity, float mean_relative_depth,
    std::span<float> source_term, std::span<float> source_diagonal,
    std::size_t num_directions, std::size_t num_extended_frequencies,
    std::size_t num_auxiliary_bins, std::span<const float> radian_frequency,
    float tail_factor, float dia_coefficient, float depth_coefficient_1,
    float depth_coefficient_2, float depth_coefficient_3,
    float depth_conversion_factor, float min_relative_depth, float inv_two_pi,
    std::span<const int> plus_address_q1_c1,
    std::span<const int> plus_address_q1_c2,
    std::span<const int> plus_address_q1_c3,
    std::span<const int> plus_address_q1_c4,
    std::span<const int> minus_address_q1_c1,
    std::span<const int> minus_address_q1_c2,
    std::span<const int> minus_address_q1_c3,
    std::span<const int> minus_address_q1_c4,
    std::span<const int> plus_address_q2_c1,
    std::span<const int> plus_address_q2_c2,
    std::span<const int> plus_address_q2_c3,
    std::span<const int> plus_address_q2_c4,
    std::span<const int> minus_address_q2_c1,
    std::span<const int> minus_address_q2_c2,
    std::span<const int> minus_address_q2_c3,
    std::span<const int> minus_address_q2_c4,
    std::span<const int> fold_address_c1_q1,
    std::span<const int> fold_address_c1_q2,
    std::span<const int> fold_address_c2_q1,
    std::span<const int> fold_address_c2_q2,
    std::span<const int> fold_address_c3_q1,
    std::span<const int> fold_address_c3_q2,
    std::span<const int> fold_address_c4_q1,
    std::span<const int> fold_address_c4_q2,
    std::span<const int> fold_address_c5_q1,
    std::span<const int> fold_address_c5_q2,
    std::span<const int> fold_address_c6_q1,
    std::span<const int> fold_address_c6_q2,
    std::span<const int> fold_address_c7_q1,
    std::span<const int> fold_address_c7_q2,
    std::span<const int> fold_address_c8_q1,
    std::span<const int> fold_address_c8_q2, float lambda_weight_1,
    float lambda_weight_2, float lambda_weight_3,
    std::span<const float> frequency_scaling, float interp_weight_1,
    float interp_weight_2, float interp_weight_3, float interp_weight_4,
    float interp_weight_5, float interp_weight_6, float interp_weight_7,
    float interp_weight_8, float squared_weight_1, float squared_weight_2,
    float squared_weight_3, float squared_weight_4, float squared_weight_5,
    float squared_weight_6, float squared_weight_7,
    float squared_weight_8) noexcept;

} // namespace ww4_nl1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3snl1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_nl1::compute_w3snl1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_nl1::compute_w3snl1. Bound in w3snl1md.F90 as
 * BIND(C, name="w3snl1_cpp").
 * @param[in] a Wave action density spectrum N(k,theta) in 1-D form (size
 * nspec).
 * @param[in] cg Group velocities (size nk).
 * @param[in] kdmean Mean relative depth k_hat d.
 * @param[out] s Nonlinear interaction source term (size nspec).
 * @param[out] d Diagonal term of the derivative (size nspec).
 * @param[in] nk Number of wavenumbers.
 * @param[in] nth Number of directions.
 * @param[in] nspec Number of spectral bins (nth * nk).
 * @param[in] nfrhgh Number of frequency bands of the extended spectrum.
 * @param[in] nspecx Size of the interaction arrays.
 * @param[in] nspecy Size of the extended auxiliary spectrum.
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] fachfe Tail extension factor.
 * @param[in] snlc1 Scaled proportionality constant.
 * @param[in] snls1 Constant c_1 in the shallow water factor.
 * @param[in] snls2 Constant c_2 in the shallow water factor.
 * @param[in] snls3 Constant c_3 in the shallow water factor.
 * @param[in] kdcon Conversion factor for the relative depth.
 * @param[in] kdmn Minimum relative depth.
 * @param[in] tpiinv Inverse of 2*Pi.
 * @param[in] ip11 Plus address array, see plus_address_q1_c1 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip12 Plus address array, see plus_address_q1_c2 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip13 Plus address array, see plus_address_q1_c3 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip14 Plus address array, see plus_address_q1_c4 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im11 Minus address array, see minus_address_q1_c1 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im12 Minus address array, see minus_address_q1_c2 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im13 Minus address array, see minus_address_q1_c3 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im14 Minus address array, see minus_address_q1_c4 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip21 Plus address array, see plus_address_q2_c1 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip22 Plus address array, see plus_address_q2_c2 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip23 Plus address array, see plus_address_q2_c3 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ip24 Plus address array, see plus_address_q2_c4 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im21 Minus address array, see minus_address_q2_c1 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im22 Minus address array, see minus_address_q2_c2 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im23 Minus address array, see minus_address_q2_c3 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] im24 Minus address array, see minus_address_q2_c4 in
 * ww4_nl1::compute_w3snl1 (size nspecx).
 * @param[in] ic11 Fold address array, see fold_address_c1_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic12 Fold address array, see fold_address_c1_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic21 Fold address array, see fold_address_c2_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic22 Fold address array, see fold_address_c2_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic31 Fold address array, see fold_address_c3_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic32 Fold address array, see fold_address_c3_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic41 Fold address array, see fold_address_c4_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic42 Fold address array, see fold_address_c4_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic51 Fold address array, see fold_address_c5_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic52 Fold address array, see fold_address_c5_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic61 Fold address array, see fold_address_c6_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic62 Fold address array, see fold_address_c6_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic71 Fold address array, see fold_address_c7_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic72 Fold address array, see fold_address_c7_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic81 Fold address array, see fold_address_c8_q1 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] ic82 Fold address array, see fold_address_c8_q2 in
 * ww4_nl1::compute_w3snl1 (size nspec).
 * @param[in] dal1 Weight 1 / (1 + lambda)^4.
 * @param[in] dal2 Weight 1 / (1 - lambda)^4.
 * @param[in] dal3 Weight 2 * DAL1 * DAL2.
 * @param[in] af11 Scaling array f_r^11 (size nspecx).
 * @param[in] awg1 Interpolation weight 1.
 * @param[in] awg2 Interpolation weight 2.
 * @param[in] awg3 Interpolation weight 3.
 * @param[in] awg4 Interpolation weight 4.
 * @param[in] awg5 Interpolation weight 5.
 * @param[in] awg6 Interpolation weight 6.
 * @param[in] awg7 Interpolation weight 7.
 * @param[in] awg8 Interpolation weight 8.
 * @param[in] swg1 Squared interpolation weight 1.
 * @param[in] swg2 Squared interpolation weight 2.
 * @param[in] swg3 Squared interpolation weight 3.
 * @param[in] swg4 Squared interpolation weight 4.
 * @param[in] swg5 Squared interpolation weight 5.
 * @param[in] swg6 Squared interpolation weight 6.
 * @param[in] swg7 Squared interpolation weight 7.
 * @param[in] swg8 Squared interpolation weight 8.
 * @pre nk > 0 and all arrays are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-14
 * @date Last update: 2026-11-19
 */
void w3snl1_cpp(const float *a, const float *cg, float kdmean, float *s,
                float *d, int nk, int nth, int nspec, int nfrhgh, int nspecx,
                int nspecy, const float *sig, float fachfe, float snlc1,
                float snls1, float snls2, float snls3, float kdcon, float kdmn,
                float tpiinv, const int *ip11, const int *ip12, const int *ip13,
                const int *ip14, const int *im11, const int *im12,
                const int *im13, const int *im14, const int *ip21,
                const int *ip22, const int *ip23, const int *ip24,
                const int *im21, const int *im22, const int *im23,
                const int *im24, const int *ic11, const int *ic12,
                const int *ic21, const int *ic22, const int *ic31,
                const int *ic32, const int *ic41, const int *ic42,
                const int *ic51, const int *ic52, const int *ic61,
                const int *ic62, const int *ic71, const int *ic72,
                const int *ic81, const int *ic82, float dal1, float dal2,
                float dal3, const float *af11, float awg1, float awg2,
                float awg3, float awg4, float awg5, float awg6, float awg7,
                float awg8, float swg1, float swg2, float swg3, float swg4,
                float swg5, float swg6, float swg7, float swg8) noexcept;

} // extern "C"

#endif // WW4_NL1_WW3_H
