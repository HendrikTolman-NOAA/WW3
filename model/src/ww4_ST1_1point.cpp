// WAVEWATCH IV (WW4) Source Code
// Copyright 2026 National Weather Service (NWS), National Oceanic and Atmospheric Administration. All rights reserved.
// NWS often uses Generative AI (GenAI) for code development and refactoring. Whenever GenAI is used, NWS requires a full human review of code before it is added to its repositories.

/**
 * @file ww4_ST1_1point.cpp
 * @brief Pure C++ implementations of the ST1 source term computations (mean parameters, input and dissipation).
 * @details This file contains C++20 optimized implementations of WAVEWATCH III (WW3) ST1 source terms
 *          intended as the starting point of the WAVEWATCH IV (WW4) conversion.
 * @author Main Authors: Aldgisl (AI Persona), Hendrik Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI)
 * @author Original WW3 Author: H. L. Tolman
 * @date 2026-08-12 (Initial Date)
 * @date 2026-11-17 (Most Recent Update Date)
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

namespace ww4_st1 {

/**
 * @brief Performs the actual computation for the ST1 input source term (W3SIN1).
 * @details Calculates the diagonal of the input source term and the source term itself using the WAM-3 method.
 * @param[in] A Action density spectrum (1-D array of size nspec).
 * @param[in] K Wavenumber for the entire spectrum (1-D array of size nspec).
 * @param[in] ustar Friction velocity.
 * @param[in] cosu Cosine of the direction of ustar.
 * @param[in] sinu Sine of the direction of ustar.
 * @param[out] S Computed input source term (1-D array of size nspec).
 * @param[out] D Computed diagonal term of derivative (1-D array of size nspec).
 * @param[in] nspec Number of spectral bins.
 * @param[in] sinc1 Input calibration constant.
 * @param[in] sig2 Group of frequencies squared (1-D array of size nspec).
 * @param[in] ecos Cosine of direction bins (1-D array of size nspec).
 * @param[in] esin Sine of direction bins (1-D array of size nspec).
 * @pre All input and output arrays must be allocated and have at least `nspec` elements.
 * @post The arrays S and D are populated with the calculated input source terms.
 */
void compute_w3sin1(
    const float* const A,
    const float* const K,
    const float ustar,
    const float cosu,
    const float sinu,
    float* const S,
    float* const D,
    const int nspec,
    const float sinc1,
    const float* const sig2,
    const float* const ecos,
    const float* const esin
) noexcept {
    // 2.  Diagonal and Source via C++
    for (int is = 0; is < nspec; ++is) {
        const float term = ustar * (ecos[is] * cosu + esin[is] * sinu) * K[is] / sig2[is] - 0.035714f;
        D[is] = sinc1 * sig2[is] * std::max(0.0f, term);
        S[is] = D[is] * A[is];
    }
}

/**
 * @brief Performs the actual computation for the ST1 dissipation source term (W3SDS1).
 * @details Calculates the whitecapping source term and diagonal term of its derivative using the WAM-3 method.
 * @param[in] A Action density spectrum (1-D array of size nspec).
 * @param[in] K Wavenumber for the entire spectrum (1-D array of size nspec).
 * @param[in] emean Mean wave energy.
 * @param[in] fmean Mean wave frequency.
 * @param[in] wnmean Mean wavenumber.
 * @param[out] S Computed dissipation source term (1-D array of size nspec).
 * @param[out] D Computed diagonal term of derivative (1-D array of size nspec).
 * @param[in] nspec Number of spectral bins.
 * @param[in] sdsc1 Dissipation calibration constant.
 * @pre All input and output arrays must be allocated and have at least `nspec` elements.
 * @post The arrays S and D are populated with the calculated dissipation source terms.
 */
void compute_w3sds1(
    const float* const A,
    const float* const K,
    const float emean,
    const float fmean,
    const float wnmean,
    float* const S,
    float* const D,
    const int nspec,
    const float sdsc1
) noexcept {
    // 1.  Common factor (calculated in C++)
    const float factor = sdsc1 * fmean * (wnmean * wnmean * wnmean) * (emean * emean);

    // 3.  Source term via C++
    for (int is = 0; is < nspec; ++is) {
        D[is] = factor * K[is];
        S[is] = D[is] * A[is];
    }
}

/**
 * @brief Computes the mean wave parameters used by the ST1 source terms (W3SPR1).
 * @details Calculates the mean wave energy, mean frequency, mean wavenumber and maximum action density
 *          of a spectrum (WAM-3 method), including the parametric tail beyond the discrete spectrum.
 *          The summation order matches the Fortran original so results are reproducible.
 *          The 2-D spectrum uses Fortran (column-major) layout, A(ith, ik) = A[ith + nth * ik].
 * @param[in] A Action density as a function of direction and wavenumber (size nth * nk).
 * @param[in] cg Group velocities (size nk).
 * @param[in] wn Wavenumbers (size nk).
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Frequency band width factors, DDEN (size nk).
 * @param[in] nth Number of directions.
 * @param[in] fte Factor in tail integration of energy.
 * @param[in] ftf Factor in tail integration of frequency.
 * @param[in] ftwn Factor in tail integration of wavenumber.
 * @param[in] tpiinv Inverse of 2*Pi.
 * @param[out] emean Mean wave energy.
 * @param[out] fmean Mean wave frequency.
 * @param[out] wnmean Mean wavenumber.
 * @param[out] amax Maximum action density in the spectrum.
 * @pre cg, wn, sig and dden have the same size nk > 0, and A has size nth * nk.
 * @post emean, fmean, wnmean and amax are set.
 */
inline void compute_w3spr1(
    const std::span<const float> A,
    const std::span<const float> cg,
    const std::span<const float> wn,
    const std::span<const float> sig,
    const std::span<const float> dden,
    const int nth,
    const float fte,
    const float ftf,
    const float ftwn,
    const float tpiinv,
    float& emean,
    float& fmean,
    float& wnmean,
    float& amax
) noexcept {
    constexpr float kMinSum = 1.0e-7f;

    const std::size_t nk = cg.size();
    const std::size_t nthu = static_cast<std::size_t>(nth);

    float e = 0.0f;
    float f = 0.0f;
    float w = 0.0f;
    float amx = 0.0f;
    float eb = 0.0f;

    for (std::size_t ik = 0; ik < nk; ++ik) {
        // 1.  Integral over directions
        eb = 0.0f;
        for (std::size_t ith = 0; ith < nthu; ++ith) {
            const float a = A[ith + nthu * ik];
            eb += a;
            amx = std::max(amx, a);
        }
        // 2.  Integrate over frequencies
        eb = eb * dden[ik] / cg[ik];
        e += eb;
        f += eb / sig[ik];
        w += eb / std::sqrt(wn[ik]);
    }

    // 3.  Add tail beyond discrete spectrum (eb holds the last band)
    const float eband = eb / dden[nk - 1];
    e += eband * fte;
    f += eband * ftf;
    w += eband * ftwn;

    // 4.  Final processing
    fmean = tpiinv * e / std::max(kMinSum, f);
    const float wnm = e / std::max(kMinSum, w);
    wnmean = wnm * wnm;
    emean = e;
    amax = amx;
}

} // namespace ww4_st1

// C-compatible ABI interfaces for Fortran interoperability
extern "C" {

/**
 * @brief Fortran interface to compute_w3sin1.
 */
void w3sin1_cpp(
    const float* const A,
    const float* const K,
    const float ustar,
    const float cosu,
    const float sinu,
    float* const S,
    float* const D,
    const int nspec,
    const float sinc1,
    const float* const sig2,
    const float* const ecos,
    const float* const esin
) noexcept {
    ww4_st1::compute_w3sin1(A, K, ustar, cosu, sinu, S, D, nspec, sinc1, sig2, ecos, esin);
}

/**
 * @brief Fortran interface to compute_w3sds1.
 */
void w3sds1_cpp(
    const float* const A,
    const float* const K,
    const float emean,
    const float fmean,
    const float wnmean,
    float* const S,
    float* const D,
    const int nspec,
    const float sdsc1
) noexcept {
    ww4_st1::compute_w3sds1(A, K, emean, fmean, wnmean, S, D, nspec, sdsc1);
}

/**
 * @brief Fortran interface to compute_w3spr1.
 * @details A is the Fortran array A(nth, nk). Output scalars are passed by reference (Fortran default).
 */
void w3spr1_cpp(
    const float* const A,
    const float* const cg,
    const float* const wn,
    const float* const sig,
    const float* const dden,
    const int nth,
    const int nk,
    const float fte,
    const float ftf,
    const float ftwn,
    const float tpiinv,
    float* const emean,
    float* const fmean,
    float* const wnmean,
    float* const amax
) noexcept {
    const auto n = static_cast<std::size_t>(nk);
    const auto nspec = static_cast<std::size_t>(nth) * n;
    ww4_st1::compute_w3spr1(
        std::span<const float>(A, nspec), std::span<const float>(cg, n),
        std::span<const float>(wn, n), std::span<const float>(sig, n),
        std::span<const float>(dden, n), nth, fte, ftf, ftwn, tpiinv,
        *emean, *fmean, *wnmean, *amax);
}

}
