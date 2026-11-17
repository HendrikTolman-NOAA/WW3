/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file ww4_ST1_1point.cpp
 * @brief Single-point WAM-3 (ST1) source term routines for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementations of the ST1 mean wave parameters (W3SPR1),
 * wind input (W3SIN1) and whitecapping dissipation (W3SDS1) for one spatial
 * point, with C-compatible interfaces called from the WAVEWATCH III (WW3)
 * Fortran module w3src1md.F90. Starting point of the WW4 conversion of the
 * ST1 source terms.
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
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-17
 */

// Library
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

namespace ww4_st1 {

// --- Constants --------------------------------------------------------------
/// Offset in the wind input diagonal term of the WAM-3 method.
constexpr float kSin1Offset = 0.035714F;

/// Lower limit on the frequency and wavenumber integrals in the mean parameters (avoids division by zero).
constexpr float kMinMeanSum = 1.0e-7F;

// --- compute_w3sin1 ---------------------------------------------------------
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
 * @param[in] sinc1 Input calibration constant.
 * @param[in] sig2 Group of frequencies squared (1-D array of size nspec).
 * @param[in] ecos Cosine of direction bins (1-D array of size nspec).
 * @param[in] esin Sine of direction bins (1-D array of size nspec).
 * @pre All spans have the same size nspec.
 * @post The arrays S and D are populated with the calculated input source terms.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SIN1 by Hendrik L. Tolman.
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-17
 */
constexpr void compute_w3sin1(const std::span<const float> A, const std::span<const float> K, const float ustar,
                              const float cosu, const float sinu, const std::span<float> S, const std::span<float> D,
                              const float sinc1, const std::span<const float> sig2, const std::span<const float> ecos,
                              const std::span<const float> esin) noexcept {
    // === Section: Diagonal and source term ==================================
    for (std::size_t is = 0; is < A.size(); ++is) {
        const float term = (((ustar * ((ecos[is] * cosu) + (esin[is] * sinu))) * K[is]) / sig2[is]) - kSin1Offset;
        D[is] = (sinc1 * sig2[is]) * std::max(0.0F, term);
        S[is] = D[is] * A[is];
    }
}

// --- compute_w3sds1 ---------------------------------------------------------
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
 * @param[in] sdsc1 Dissipation calibration constant.
 * @pre All spans have the same size nspec.
 * @post The arrays S and D are populated with the calculated dissipation source terms.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SDS1 by Hendrik L. Tolman.
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-17
 */
constexpr void compute_w3sds1(const std::span<const float> A, const std::span<const float> K, const float emean,
                              const float fmean, const float wnmean, const std::span<float> S, const std::span<float> D,
                              const float sdsc1) noexcept {
    // === Section: Common factor =============================================
    const float factor = ((sdsc1 * fmean) * (wnmean * wnmean * wnmean)) * (emean * emean);

    // === Section: Diagonal and source term ==================================
    for (std::size_t is = 0; is < A.size(); ++is) {
        D[is] = factor * K[is];
        S[is] = D[is] * A[is];
    }
}

// --- compute_w3spr1 ---------------------------------------------------------
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
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SPR1 by Hendrik L. Tolman.
 * @date Initial: 2026-11-17
 * @date Last update: 2026-11-17
 */
void compute_w3spr1(const std::span<const float> A, const std::span<const float> cg, const std::span<const float> wn,
                    const std::span<const float> sig, const std::span<const float> dden, const std::size_t nth,
                    const float fte, const float ftf, const float ftwn, const float tpiinv, float &emean, float &fmean,
                    float &wnmean, float &amax) noexcept {
    const std::size_t nk = cg.size();

    float esum = 0.0F;
    float fsum = 0.0F;
    float wsum = 0.0F;
    float amx = 0.0F;
    float eband = 0.0F;

    for (std::size_t ik = 0; ik < nk; ++ik) {
        // === Section: Integral over directions ==============================
        eband = 0.0F;
        for (std::size_t ith = 0; ith < nth; ++ith) {
            const float act = A[ith + (nth * ik)];
            eband += act;
            amx = std::max(amx, act);
        }
        // === Section: Integral over frequencies =============================
        eband = (eband * dden[ik]) / cg[ik];
        esum += eband;
        fsum += eband / sig[ik];
        wsum += eband / std::sqrt(wn[ik]);
    }

    // === Section: Tail beyond discrete spectrum =============================
    // eband holds the energy of the last discrete band
    const float etail = eband / dden[nk - 1];
    esum += etail * fte;
    fsum += etail * ftf;
    wsum += etail * ftwn;

    // === Section: Final processing ==========================================
    fmean = (tpiinv * esum) / std::max(kMinMeanSum, fsum);
    const float wnm = esum / std::max(kMinMeanSum, wsum);
    wnmean = wnm * wnm;
    emean = esum;
    amax = amx;
}

} // namespace ww4_st1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3sin1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to compute_w3sin1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to ww4_st1::compute_w3sin1.
 * @param[in] A Action density spectrum (size nspec).
 * @param[in] K Wavenumber for the entire spectrum (size nspec).
 * @param[in] ustar Friction velocity.
 * @param[in] cosu Cosine of the direction of ustar.
 * @param[in] sinu Sine of the direction of ustar.
 * @param[out] S Computed input source term (size nspec).
 * @param[out] D Computed diagonal term of derivative (size nspec).
 * @param[in] nspec Number of spectral bins.
 * @param[in] sinc1 Input calibration constant.
 * @param[in] sig2 Group of frequencies squared (size nspec).
 * @param[in] ecos Cosine of direction bins (size nspec).
 * @param[in] esin Sine of direction bins (size nspec).
 * @pre All arrays are allocated with at least nspec elements.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-17
 */
void w3sin1_cpp(const float *const A, const float *const K, const float ustar, const float cosu, const float sinu,
                float *const S, float *const D, const int nspec, const float sinc1, const float *const sig2,
                const float *const ecos, const float *const esin) noexcept {
    const auto n = static_cast<std::size_t>(nspec);
    ww4_st1::compute_w3sin1(std::span<const float>(A, n), std::span<const float>(K, n), ustar, cosu, sinu,
                            std::span<float>(S, n), std::span<float>(D, n), sinc1, std::span<const float>(sig2, n),
                            std::span<const float>(ecos, n), std::span<const float>(esin, n));
}

// --- w3sds1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to compute_w3sds1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to ww4_st1::compute_w3sds1.
 * @param[in] A Action density spectrum (size nspec).
 * @param[in] K Wavenumber for the entire spectrum (size nspec).
 * @param[in] emean Mean wave energy.
 * @param[in] fmean Mean wave frequency.
 * @param[in] wnmean Mean wavenumber.
 * @param[out] S Computed dissipation source term (size nspec).
 * @param[out] D Computed diagonal term of derivative (size nspec).
 * @param[in] nspec Number of spectral bins.
 * @param[in] sdsc1 Dissipation calibration constant.
 * @pre All arrays are allocated with at least nspec elements.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Jules (Agentic AI), Kiro (Agentic AI/SI)
 * @date Initial: 2026-08-12
 * @date Last update: 2026-11-17
 */
void w3sds1_cpp(const float *const A, const float *const K, const float emean, const float fmean, const float wnmean,
                float *const S, float *const D, const int nspec, const float sdsc1) noexcept {
    const auto n = static_cast<std::size_t>(nspec);
    ww4_st1::compute_w3sds1(std::span<const float>(A, n), std::span<const float>(K, n), emean, fmean, wnmean,
                            std::span<float>(S, n), std::span<float>(D, n), sdsc1);
}

// --- w3spr1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to compute_w3spr1.
 * @details A is the Fortran array A(nth, nk). Output scalars are passed by reference (Fortran default).
 * @param[in] A Action density as a function of direction and wavenumber (size nth * nk).
 * @param[in] cg Group velocities (size nk).
 * @param[in] wn Wavenumbers (size nk).
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Frequency band width factors (size nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] fte Factor in tail integration of energy.
 * @param[in] ftf Factor in tail integration of frequency.
 * @param[in] ftwn Factor in tail integration of wavenumber.
 * @param[in] tpiinv Inverse of 2*Pi.
 * @param[out] emean Mean wave energy.
 * @param[out] fmean Mean wave frequency.
 * @param[out] wnmean Mean wavenumber.
 * @param[out] amax Maximum action density in the spectrum.
 * @pre nk > 0 and all arrays are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-17
 * @date Last update: 2026-11-17
 */
void w3spr1_cpp(const float *const A, const float *const cg, const float *const wn, const float *const sig,
                const float *const dden, const int nth, const int nk, const float fte, const float ftf,
                const float ftwn, const float tpiinv, float *const emean, float *const fmean, float *const wnmean,
                float *const amax) noexcept {
    const auto nkz = static_cast<std::size_t>(nk);
    const auto nthz = static_cast<std::size_t>(nth);
    ww4_st1::compute_w3spr1(std::span<const float>(A, nthz * nkz), std::span<const float>(cg, nkz),
                            std::span<const float>(wn, nkz), std::span<const float>(sig, nkz),
                            std::span<const float>(dden, nkz), nthz, fte, ftf, ftwn, tpiinv, *emean, *fmean, *wnmean,
                            *amax);
}

} // extern "C"
