// WAVEWATCH IV (WW4) Source Code
// Copyright 2026 National Weather Service (NWS), National Oceanic and Atmospheric Administration. All rights reserved.
// NWS often uses Generative AI (GenAI) for code development and refactoring. Whenever GenAI is used, NWS requires a full human review of code before it is added to its repositories.

/**
 * @file ww4_LN1_1point.cpp
 * @brief Pure C++ implementations of the LN1 linear wind input source term computation.
 * @details This file contains C++20 optimized implementations of WAVEWATCH III (WW3) LN1 source term
 *          intended as part of the WAVEWATCH IV (WW4) conversion.
 *          Linear wind input according to Cavaleri and Melanotte-Rizzoli (1982) filtered for low
 *          frequencies according to Tolman (1992).
 * @author Main Authors: Aldgisl (AI Persona), Hendrik Tolman
 * @author Contributors: Jules (Agentic AI)
 * @author Original WW3 Author: H. L. Tolman
 * @date 2026-08-14 (Initial Date)
 * @date 2026-08-14 (Most Recent Update Date)
 */

#include <cmath>
#include <algorithm>
#include <vector>

namespace ww4_ln1 {

/**
 * @brief Performs the actual computation for the LN1 linear wind input source term (W3SLN1).
 * @details Calculates the linear wind input source term spectrum according to Cavaleri and Melanotte-Rizzoli (1982)
 *          with low frequency filtering (Tolman 1992).
 * @param[in] K Wavenumber for the entire spectrum (1-D array of size nk).
 * @param[in] fhigh Cut-off frequency in integration (rad/s).
 * @param[in] ustar Friction velocity.
 * @param[in] usdir Wind direction (radians).
 * @param[out] S Computed linear input source term (1-D array of size nth * nk, indexed as S[ith + ik * nth]).
 * @param[in] nth Number of direction bins.
 * @param[in] nk Number of frequency bins.
 * @param[in] ecos Cosine of direction bins (1-D array of size nth).
 * @param[in] esin Sine of direction bins (1-D array of size nth).
 * @param[in] sig Frequency array in rad/s (1-D array of size nk).
 * @param[in] slnc1 Linear input constant.
 * @param[in] fspm PM frequency factor.
 * @param[in] fshf High-frequency cut-off factor.
 * @param[in] grav Acceleration of gravity.
 * @pre All input and output arrays must be allocated with appropriate size (nth, nk, or nth * nk).
 * @post The array S is populated with the calculated linear input source terms.
 */
void compute_w3sln1(
    const float* const K,
    const float fhigh,
    const float ustar,
    const float usdir,
    float* const S,
    const int nth,
    const int nk,
    const float* const ecos,
    const float* const esin,
    const float* const sig,
    const float slnc1,
    const float fspm,
    const float fshf,
    const float grav
) noexcept {
    const int nspec = nth * nk;
    if (nspec <= 0) {
        return;
    }

    // Handle zero or negative ustar gracefully by zeroing out S
    if (ustar <= 0.0f) {
        std::fill_n(S, nspec, 0.0f);
        return;
    }

    const float cosu = std::cos(usdir);
    const float sinu = std::sin(usdir);

    // 1. Directional factor DIRF(ITH)
    std::vector<float> dirf(nth);
    for (int ith = 0; ith < nth; ++ith) {
        const float proj = std::max(0.0f, ecos[ith] * cosu + esin[ith] * sinu);
        const float proj2 = proj * proj;
        dirf[ith] = proj2 * proj2;
    }

    // 2. Frequency filtering factors
    const float fac = slnc1 * (ustar * ustar * ustar * ustar);
    const float ff1 = fspm * grav / (28.0f * ustar);
    const float ff2 = fshf * std::min(sig[nk - 1], fhigh);
    const float ffilt = std::min(std::max(ff1, ff2), 2.0f * sig[nk - 1]);

    std::vector<float> wnf(nk);
    for (int ik = 0; ik < nk; ++ik) {
        const float rfr = sig[ik] / ffilt;
        if (rfr < 0.5f) {
            wnf[ik] = 0.0f;
        } else {
            const float rfr2 = rfr * rfr;
            const float rfr4 = rfr2 * rfr2;
            wnf[ik] = (fac / K[ik]) * std::exp(-1.0f / rfr4);
        }
    }

    // 3. Compose source term S(ITH, IK) -> S[ith + ik * nth]
    for (int ik = 0; ik < nk; ++ik) {
        const float wnf_val = wnf[ik];
        const int ik_offset = ik * nth;
        for (int ith = 0; ith < nth; ++ith) {
            S[ith + ik_offset] = wnf_val * dirf[ith];
        }
    }
}

} // namespace ww4_ln1

// C-compatible ABI interfaces for Fortran interoperability
extern "C" {

/**
 * @brief Fortran interface to compute_w3sln1.
 */
void w3sln1_cpp(
    const float* const K,
    const float fhigh,
    const float ustar,
    const float usdir,
    float* const S,
    const int nth,
    const int nk,
    const float* const ecos,
    const float* const esin,
    const float* const sig,
    const float slnc1,
    const float fspm,
    const float fshf,
    const float grav
) noexcept {
    ww4_ln1::compute_w3sln1(K, fhigh, ustar, usdir, S, nth, nk, ecos, esin, sig, slnc1, fspm, fshf, grav);
}

}
