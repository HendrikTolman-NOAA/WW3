// WAVEWATCH IV (WW4) Source Code
// Copyright 2024 National Weather Service (NWS), National Oceanic and Atmospheric Administration. All rights reserved.
// NWS often uses Generative AI (GenAI) for code development and refactoring. Whenever GenAI is used, NWS requires a full human review of code before it is added to its repositories.

/**
 * @file w3src1_cpp.cpp
 * @brief Pure C++ implementations of the ST1 source term computations (input and dissipation).
 * @details This file contains C++20 optimized implementations of WAVEWATCH III (WW3) ST1 source terms
 *          intended as the starting point of the WAVEWATCH IV (WW4) conversion.
 * @author Main Authors: Aldgisl (AI Persona), Jules (Developer)
 * @author Original WW3 Author: H. L. Tolman
 * @date 2024-08-03 (Initial Date)
 * @date 2024-08-03 (Most Recent Update Date)
 */

#include <cmath>
#include <algorithm>

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
    const float factor = sdsc1 * fmean * (wnmean * wnmean * wnmean) * (emean * emean);
    for (int is = 0; is < nspec; ++is) {
        D[is] = factor * K[is];
        S[is] = D[is] * A[is];
    }
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

}
