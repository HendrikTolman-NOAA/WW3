// WAVEWATCH IV (WW4) Source Code
// Copyright 2026 National Weather Service (NWS), National Oceanic and Atmospheric Administration. All rights reserved.
// NWS often uses Generative AI (GenAI) for code development and refactoring. Whenever GenAI is used, NWS requires a full human review of code before it is added to its repositories.

/**
 * @file ww4_NL1_1point.cpp
 * @brief Pure C++ implementation of the NL1 nonlinear interaction source term computation (W3SNL1).
 * @details This file contains a C++20 optimized implementation of the WAVEWATCH III (WW3) NL1 Discrete
 *          Interaction Approximation (DIA) source term intended for WAVEWATCH IV (WW4).
 * @author Main Authors: Aldgisl (AI Persona), Hendrik Tolman
 * @author Contributors: Jules (Agentic AI)
 * @author Original WW3 Author: H. L. Tolman
 * @date 2026-08-14 (Initial Date)
 * @date 2026-08-14 (Most Recent Update Date)
 */

#include <cmath>
#include <algorithm>
#include <vector>

namespace ww4_nl1 {

/**
 * @brief Performs the actual computation for the NL1 nonlinear interaction source term (W3SNL1).
 * @details Calculates the nonlinear wave-wave interactions and diagonal term of its derivative
 *          according to the Discrete Interaction Approximation (DIA) method.
 * @param[in] A Action density spectrum (1-D array of size nspec).
 * @param[in] CG Group velocities (1-D array of size nk).
 * @param[in] kdmean Mean relative depth.
 * @param[out] S Computed source term (1-D array of size nspec).
 * @param[out] D Computed diagonal term of derivative (1-D array of size nspec).
 * @param[in] nk Number of spectral frequencies.
 * @param[in] nth Number of spectral directions.
 * @param[in] nspec Total number of spectral bins (nk * nth).
 * @param[in] nfrhgh High frequency spectral limit.
 * @param[in] nspecx Address space size for quadruplet interaction calculation.
 * @param[in] nspecy Extended spectral size including high-frequency tail.
 * @param[in] sig Frequencies array (1-D array of size nk).
 * @param[in] fachfe Tail extension factor.
 * @param[in] snlc1 DIA constant 1.
 * @param[in] snls1 DIA constant 2.
 * @param[in] snls2 DIA constant 3.
 * @param[in] snls3 DIA constant 4.
 * @param[in] kdcon Shallow water correction constant 1.
 * @param[in] kdmn Shallow water correction constant 2.
 * @param[in] IP11 Direction/frequency index mapping arrays (each size nspecx).
 * @param[in] IP12 Index mapping array.
 * @param[in] IP13 Index mapping array.
 * @param[in] IP14 Index mapping array.
 * @param[in] IM11 Index mapping array.
 * @param[in] IM12 Index mapping array.
 * @param[in] IM13 Index mapping array.
 * @param[in] IM14 Index mapping array.
 * @param[in] IP21 Index mapping array.
 * @param[in] IP22 Index mapping array.
 * @param[in] IP23 Index mapping array.
 * @param[in] IP24 Index mapping array.
 * @param[in] IM21 Index mapping array.
 * @param[in] IM22 Index mapping array.
 * @param[in] IM23 Index mapping array.
 * @param[in] IM24 Index mapping array.
 * @param[in] IC11 Output source/derivative index mapping arrays (each size nspec).
 * @param[in] IC12 Index mapping array.
 * @param[in] IC21 Index mapping array.
 * @param[in] IC22 Index mapping array.
 * @param[in] IC31 Index mapping array.
 * @param[in] IC32 Index mapping array.
 * @param[in] IC41 Index mapping array.
 * @param[in] IC42 Index mapping array.
 * @param[in] IC51 Index mapping array.
 * @param[in] IC52 Index mapping array.
 * @param[in] IC61 Index mapping array.
 * @param[in] IC62 Index mapping array.
 * @param[in] IC71 Index mapping array.
 * @param[in] IC72 Index mapping array.
 * @param[in] IC81 Index mapping array.
 * @param[in] IC82 Index mapping array.
 * @param[in] dal1 Weight factor 1.
 * @param[in] dal2 Weight factor 2.
 * @param[in] dal3 Weight factor 3.
 * @param[in] af11 Scaling array (1-D array of size nspecx).
 * @param[in] awg1 Angular interpolation weight 1.
 * @param[in] awg2 Angular interpolation weight 2.
 * @param[in] awg3 Angular interpolation weight 3.
 * @param[in] awg4 Angular interpolation weight 4.
 * @param[in] awg5 Angular interpolation weight 5.
 * @param[in] awg6 Angular interpolation weight 6.
 * @param[in] awg7 Angular interpolation weight 7.
 * @param[in] awg8 Angular interpolation weight 8.
 * @param[in] swg1 Derivative weight 1.
 * @param[in] swg2 Derivative weight 2.
 * @param[in] swg3 Derivative weight 3.
 * @param[in] swg4 Derivative weight 4.
 * @param[in] swg5 Derivative weight 5.
 * @param[in] swg6 Derivative weight 6.
 * @param[in] swg7 Derivative weight 7.
 * @param[in] swg8 Derivative weight 8.
 * @pre All input arrays must be properly populated and allocated to their respective dimensions.
 * @post S and D arrays are populated with calculated DIA nonlinear source and diagonal derivative terms.
 */
void compute_w3snl1(
    const float* const A,
    const float* const CG,
    const float kdmean,
    float* const S,
    float* const D,
    const int nk,
    const int nth,
    const int nspec,
    const int nfrhgh,
    const int nspecx,
    const int nspecy,
    const float* const sig,
    const float fachfe,
    const float snlc1,
    const float snls1,
    const float snls2,
    const float snls3,
    const float kdcon,
    const float kdmn,
    const int* const IP11,
    const int* const IP12,
    const int* const IP13,
    const int* const IP14,
    const int* const IM11,
    const int* const IM12,
    const int* const IM13,
    const int* const IM14,
    const int* const IP21,
    const int* const IP22,
    const int* const IP23,
    const int* const IP24,
    const int* const IM21,
    const int* const IM22,
    const int* const IM23,
    const int* const IM24,
    const int* const IC11,
    const int* const IC12,
    const int* const IC21,
    const int* const IC22,
    const int* const IC31,
    const int* const IC32,
    const int* const IC41,
    const int* const IC42,
    const int* const IC51,
    const int* const IC52,
    const int* const IC61,
    const int* const IC62,
    const int* const IC71,
    const int* const IC72,
    const int* const IC81,
    const int* const IC82,
    const float dal1,
    const float dal2,
    const float dal3,
    const float* const af11,
    const float awg1,
    const float awg2,
    const float awg3,
    const float awg4,
    const float awg5,
    const float awg6,
    const float awg7,
    const float awg8,
    const float swg1,
    const float swg2,
    const float swg3,
    const float swg4,
    const float swg5,
    const float swg6,
    const float swg7,
    const float swg8
) noexcept {
    constexpr float tpiinv = 0.5f / 3.14159265358979323846f;

    // 1. Calculate proportionality constant
    const float x = std::max(kdcon * kdmean, kdmn);
    const float x2 = std::max(-1.0e15f, snls3 * x);
    const float cons = snlc1 * (1.0f + (snls1 / x) * (1.0f - snls2 * x) * std::exp(x2));

    // 2. Prepare auxiliary spectrum and arrays
    // Buffer for UE (size nth + nspecy to support 1-based Fortran indices from 1-nth to nspecy)
    std::vector<float> ue_buf(nth + nspecy, 0.0f);
    float* const UE = ue_buf.data() + nth - 1;

    std::vector<float> CON(nspec, 0.0f);

    for (int ifr = 0; ifr < nk; ++ifr) {
        const float conx = tpiinv / sig[ifr] * CG[ifr];
        for (int ith = 0; ith < nth; ++ith) {
            const int isp = ith + ifr * nth;
            UE[isp + 1] = A[isp] / conx;
            CON[isp] = conx;
        }
    }

    for (int ifr = nk; ifr < nfrhgh; ++ifr) {
        for (int ith = 0; ith < nth; ++ith) {
            const int isp = ith + ifr * nth;
            UE[isp + 1] = UE[isp + 1 - nth] * fachfe;
        }
    }

    // Buffers for DIA interaction arrays (size nth + nspecx for 1-based Fortran indexing)
    std::vector<float> sa1_buf(nth + nspecx, 0.0f);
    std::vector<float> sa2_buf(nth + nspecx, 0.0f);
    std::vector<float> da1c_buf(nth + nspecx, 0.0f);
    std::vector<float> da1p_buf(nth + nspecx, 0.0f);
    std::vector<float> da1m_buf(nth + nspecx, 0.0f);
    std::vector<float> da2c_buf(nth + nspecx, 0.0f);
    std::vector<float> da2p_buf(nth + nspecx, 0.0f);
    std::vector<float> da2m_buf(nth + nspecx, 0.0f);

    float* const SA1  = sa1_buf.data()  + nth - 1;
    float* const SA2  = sa2_buf.data()  + nth - 1;
    float* const DA1C = da1c_buf.data() + nth - 1;
    float* const DA1P = da1p_buf.data() + nth - 1;
    float* const DA1M = da1m_buf.data() + nth - 1;
    float* const DA2C = da2c_buf.data() + nth - 1;
    float* const DA2P = da2p_buf.data() + nth - 1;
    float* const DA2M = da2m_buf.data() + nth - 1;

    // Zero out boundary region from 1-nth to 0
    for (int k = 1 - nth; k <= 0; ++k) {
        UE[k]   = 0.0f;
        SA1[k]  = 0.0f;
        SA2[k]  = 0.0f;
        DA1C[k] = 0.0f;
        DA1P[k] = 0.0f;
        DA1M[k] = 0.0f;
        DA2C[k] = 0.0f;
        DA2P[k] = 0.0f;
        DA2M[k] = 0.0f;
    }

    // 3. Calculate interactions for extended spectrum
    for (int isp = 0; isp < nspecx; ++isp) {
        const int isp_f = isp + 1;
        const float e00 = UE[isp_f];

        const float ep1 = awg1 * UE[IP11[isp]] + awg2 * UE[IP12[isp]]
                        + awg3 * UE[IP13[isp]] + awg4 * UE[IP14[isp]];
        const float em1 = awg5 * UE[IM11[isp]] + awg6 * UE[IM12[isp]]
                        + awg7 * UE[IM13[isp]] + awg8 * UE[IM14[isp]];
        const float ep2 = awg1 * UE[IP21[isp]] + awg2 * UE[IP22[isp]]
                        + awg3 * UE[IP23[isp]] + awg4 * UE[IP24[isp]];
        const float em2 = awg5 * UE[IM21[isp]] + awg6 * UE[IM22[isp]]
                        + awg7 * UE[IM23[isp]] + awg8 * UE[IM24[isp]];

        const float factor = cons * af11[isp] * e00;

        const float sa1a = e00 * (ep1 * dal1 + em1 * dal2);
        const float sa1b = sa1a - ep1 * em1 * dal3;
        const float sa2a = e00 * (ep2 * dal1 + em2 * dal2);
        const float sa2b = sa2a - ep2 * em2 * dal3;

        SA1[isp_f] = factor * sa1b;
        SA2[isp_f] = factor * sa2b;

        DA1C[isp_f] = cons * af11[isp] * (sa1a + sa1b);
        DA1P[isp_f] = factor * (dal1 * e00 - dal3 * em1);
        DA1M[isp_f] = factor * (dal2 * e00 - dal3 * ep1);

        DA2C[isp_f] = cons * af11[isp] * (sa2a + sa2b);
        DA2P[isp_f] = factor * (dal1 * e00 - dal3 * em2);
        DA2M[isp_f] = factor * (dal2 * e00 - dal3 * ep2);
    }

    // 4. Put source and diagonal term together
    for (int isp = 0; isp < nspec; ++isp) {
        const int isp_f = isp + 1;

        S[isp] = CON[isp] * ( -2.0f * (SA1[isp_f] + SA2[isp_f])
            + awg1 * (SA1[IC11[isp]] + SA2[IC12[isp]])
            + awg2 * (SA1[IC21[isp]] + SA2[IC22[isp]])
            + awg3 * (SA1[IC31[isp]] + SA2[IC32[isp]])
            + awg4 * (SA1[IC41[isp]] + SA2[IC42[isp]])
            + awg5 * (SA1[IC51[isp]] + SA2[IC52[isp]])
            + awg6 * (SA1[IC61[isp]] + SA2[IC62[isp]])
            + awg7 * (SA1[IC71[isp]] + SA2[IC72[isp]])
            + awg8 * (SA1[IC81[isp]] + SA2[IC82[isp]]) );

        D[isp] = -2.0f * (DA1C[isp_f] + DA2C[isp_f])
            + swg1 * (DA1P[IC11[isp]] + DA2P[IC12[isp]])
            + swg2 * (DA1P[IC21[isp]] + DA2P[IC22[isp]])
            + swg3 * (DA1P[IC31[isp]] + DA2P[IC32[isp]])
            + swg4 * (DA1P[IC41[isp]] + DA2P[IC42[isp]])
            + swg5 * (DA1M[IC51[isp]] + DA2M[IC52[isp]])
            + swg6 * (DA1M[IC61[isp]] + DA2M[IC62[isp]])
            + swg7 * (DA1M[IC71[isp]] + DA2M[IC72[isp]])
            + swg8 * (DA1M[IC81[isp]] + DA2M[IC82[isp]]);
    }
}

} // namespace ww4_nl1

// C-compatible ABI interface for Fortran interoperability
extern "C" {

/**
 * @brief Fortran interface to compute_w3snl1.
 */
void w3snl1_cpp(
    const float* const A,
    const float* const CG,
    const float kdmean,
    float* const S,
    float* const D,
    const int nk,
    const int nth,
    const int nspec,
    const int nfrhgh,
    const int nspecx,
    const int nspecy,
    const float* const sig,
    const float fachfe,
    const float snlc1,
    const float snls1,
    const float snls2,
    const float snls3,
    const float kdcon,
    const float kdmn,
    const int* const IP11,
    const int* const IP12,
    const int* const IP13,
    const int* const IP14,
    const int* const IM11,
    const int* const IM12,
    const int* const IM13,
    const int* const IM14,
    const int* const IP21,
    const int* const IP22,
    const int* const IP23,
    const int* const IP24,
    const int* const IM21,
    const int* const IM22,
    const int* const IM23,
    const int* const IM24,
    const int* const IC11,
    const int* const IC12,
    const int* const IC21,
    const int* const IC22,
    const int* const IC31,
    const int* const IC32,
    const int* const IC41,
    const int* const IC42,
    const int* const IC51,
    const int* const IC52,
    const int* const IC61,
    const int* const IC62,
    const int* const IC71,
    const int* const IC72,
    const int* const IC81,
    const int* const IC82,
    const float dal1,
    const float dal2,
    const float dal3,
    const float* const af11,
    const float awg1,
    const float awg2,
    const float awg3,
    const float awg4,
    const float awg5,
    const float awg6,
    const float awg7,
    const float awg8,
    const float swg1,
    const float swg2,
    const float swg3,
    const float swg4,
    const float swg5,
    const float swg6,
    const float swg7,
    const float swg8
) noexcept {
    ww4_nl1::compute_w3snl1(
        A, CG, kdmean, S, D, nk, nth, nspec, nfrhgh, nspecx, nspecy,
        sig, fachfe, snlc1, snls1, snls2, snls3, kdcon, kdmn,
        IP11, IP12, IP13, IP14, IM11, IM12, IM13, IM14,
        IP21, IP22, IP23, IP24, IM21, IM22, IM23, IM24,
        IC11, IC12, IC21, IC22, IC31, IC32, IC41, IC42,
        IC51, IC52, IC61, IC62, IC71, IC72, IC81, IC82,
        dal1, dal2, dal3, af11,
        awg1, awg2, awg3, awg4, awg5, awg6, awg7, awg8,
        swg1, swg2, swg3, swg4, swg5, swg6, swg7, swg8
    );
}

}
