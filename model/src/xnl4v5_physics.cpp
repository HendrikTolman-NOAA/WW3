/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file xnl4v5_physics.cpp
 * @brief Implementation of the physics functions of the Xnl library for the
 * WRT (exact) nonlinear four-wave interactions.
 * @details Pure C++20 translation of the physics functions of the Xnl library
 * by G. Ph. van Vledder (X_DISPER, X_COSK, X_FLOCUS, X_LOCUS1, X_LOCUS2,
 * X_JACOBIAN, X_CPLE, XC_WEBB, XC_HH, Z_WNUMB and Z_CMPCG), with C-compatible
 * interfaces called from the corresponding Fortran functions in
 * mod_xnl4v5.f90. Second stage of the translation of the Xnl library, which
 * is kept as a separate library (namespace xnl4v5, CMake target xnl4v5_lib).
 * The functions, their use, and the parts of the WW3 code that are kept for
 * review are described in xnl4v5_physics.h. Each routine has its full
 * documentation (details, contracts and heritage) in front of its definition,
 * and a shorter version with its declaration in xnl4v5_physics.h.
 *
 * References:
 * - Webb, D. J., 1978: Non-linear transfers between sea waves. Deep-Sea
 *   Res., 25, 279-298.
 * - Herterich, K. and K. Hasselmann, 1980: A similarity relation for the
 *   nonlinear energy transfer in a finite-depth gravity-wave spectrum. J.
 *   Fluid Mech., 97, 215-224.
 * - Tracy, B. and D. T. Resio, 1982: Theory and calculation of the nonlinear
 *   energy transfer between sea waves in deep water. WES Report 11, US Army
 *   Corps of Engineers.
 * - van Vledder, G. Ph., 2006: The WRT method for the computation of
 *   non-linear four-wave interactions in discrete spectral wave models.
 *   Coastal Eng., 53, 223-242.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section NL2.
 *
 * Attribution: original Xnl library (GurboQuad version 5.03), 1999-2003,
 * Gerbrant Ph. van Vledder, ALKYON Hydraulic Consultancy & Research. No
 * separate copyright statement was included in the original source. The
 * author shared the code through WAVEWATCH III (WW3) before WW3 was licensed;
 * the history of the code is documented in the WW3 repository.
 *
 * The C++ code is a straight translation of the Fortran code, without further
 * development or refactoring of the method: the computations and their order
 * are unchanged, and the results are bitwise identical to those of the
 * Fortran code. The code is only organized in C++ functions and structures
 * as required by the C++ coding rules, and variables that were computed but
 * not used are omitted.
 *
 * The National Weather Service (NWS) often uses Generative SI (GenSI, formerly
 * identified as GenAI) for code development and refactoring. Whenever GenSI is
 * used, NWS requires a full human review of code before it is added to its
 * repositories.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from the Xnl library by Gerbrant Ph. van
 * Vledder: mod_xnl4v5.f90 (X_DISPER, X_COSK, X_FLOCUS, X_LOCUS1, X_LOCUS2,
 * X_JACOBIAN, X_CPLE, XC_WEBB, XC_HH) and serv_xnl4v5.f90 (Z_WNUMB,
 * Z_CMPCG), as included in WAVEWATCH III (WW3) for the NL2 source term.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */

// Local
#include "xnl4v5_physics.h"

// System
#include <algorithm>
#include <cmath>

namespace xnl4v5 {

namespace {

// --- Constants --------------------------------------------------------------
/// Value returned by Z_WNUMB and Z_CMPCG for zero depth or frequency.
constexpr float kInvalidValue = -10.0F;

/// k d above which Z_CMPCG uses the deep water group velocity.
constexpr float kCmpcgDeepLimit = 30.0F;

/// k d above which X_DISPER and X_JACOBIAN use the deep water relation.
constexpr float kDeepWaterLimit = 20.0F;

/// Argument above which TANZ and COSZ use the argument 25 (XC_HH).
constexpr float kHyperbolicLimit = 20.0F;

/// Argument used by TANZ and COSZ above kHyperbolicLimit (XC_HH).
constexpr float kHyperbolicCap = 25.0F;

/// Acceleration of gravity used inside XC_HH.
constexpr float kGravityHh = 9.81F;

/// Constant pi / 4 of XC_HH (DATA statement).
constexpr float kQuarterPiHh = 0.785398163F;

/// Small value in the denominators of XC_HH.
constexpr float kEpsilonHh = 1.0e-20F;

/// Small value in the denominators of XC_WEBB.
constexpr float kEpsilonWebb = 1.0e-30F;

/// Value of the locus functions for an invalid dispersion type.
constexpr float kInvalidLocus = -1.0F;

// --- tanz -------------------------------------------------------------------
/**
 * @brief tanh with the argument limited to 25 above 20 (TANZ in XC_HH).
 * @param value The argument x.
 * @return tanh(x), or tanh(25) for x > 20.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl TANZ by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float tanz(const float value) noexcept {
  return std::tanh(value > kHyperbolicLimit ? kHyperbolicCap : value);
}

// --- cosz -------------------------------------------------------------------
/**
 * @brief cosh with the argument limited to 25 above 20 (COSZ in XC_HH).
 * @param value The argument x.
 * @return cosh(x), or cosh(25) for x > 20.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl COSZ by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float cosz(const float value) noexcept {
  return std::cosh(value > kHyperbolicLimit ? kHyperbolicCap : value);
}

// --- HhTerm -----------------------------------------------------------------
/**
 * @brief Partial coupling coefficient of one permutation in XC_HH.
 */
struct HhTerm {
  float coefficient; ///< Partial coefficient SCPL.
  float omega1;      ///< Radian frequency OM1 of the permutation.
  float omega2;      ///< Radian frequency OM2 of the permutation.
  float omega3;      ///< Radian frequency OM3 of the permutation.
};

// --- HhInput ----------------------------------------------------------------
/**
 * @brief Wavenumbers and signs of one permutation in XC_HH.
 */
struct HhInput {
  float k1x0; ///< x-component of the first wavenumber.
  float k1y0; ///< y-component of the first wavenumber.
  float k2x0; ///< x-component of the second wavenumber.
  float k2y0; ///< y-component of the second wavenumber.
  float k3x0; ///< x-component of the third wavenumber.
  float k3y0; ///< y-component of the third wavenumber.
  float s1;   ///< Sign of the first wavenumber.
  float s2;   ///< Sign of the second wavenumber.
  float s3;   ///< Sign of the third wavenumber.
};

// --- hh_term ----------------------------------------------------------------
/**
 * @brief Partial coupling coefficient of one permutation (XC_HH).
 * @details One of the three blocks of XC_HH, with the order of operations of
 * the Fortran code. The term eps is added in the denominator of T1 only in
 * the first block, as in the Fortran code. Powers x**2 with the real
 * exponent 2. are computed as x * x, as gfortran does at -O2 and -O3 (the
 * WW3 release optimization). At -O0, gfortran calls powf for these, which
 * can differ in the last bit.
 * @param in Wavenumbers and signs of the permutation.
 * @param depth Water depth h.
 * @param add_epsilon Add eps to the denominator of T1 (first block).
 * @return The partial coefficient and the radian frequencies.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl XC_HH by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
HhTerm hh_term(const HhInput &in, const float depth,
               const bool add_epsilon) noexcept {
  const float g = kGravityHh;
  const float gsq = g * g;
  const float k1x = in.s1 * in.k1x0;
  const float k1y = in.s1 * in.k1y0;
  const float k2x = in.s2 * in.k2x0;
  const float k2y = in.s2 * in.k2y0;
  const float k3x = in.s3 * in.k3x0;
  const float k3y = in.s3 * in.k3y0;
  const float k1 = std::sqrt((k1x * k1x) + (k1y * k1y));
  const float k2 = std::sqrt((k2x * k2x) + (k2y * k2y));
  const float k3 = std::sqrt((k3x * k3x) + (k3y * k3y));
  const float tnz1 = tanz(k1 * depth);
  const float tnz2 = tanz(k2 * depth);
  const float tnz3 = tanz(k3 * depth);
  const float csz1 = cosz(k1 * depth);
  const float csz2 = cosz(k2 * depth);
  const float csz3 = cosz(k3 * depth);
  const float om1 = std::sqrt((g * k1) * tnz1);
  const float om2 = std::sqrt((g * k2) * tnz2);
  const float om3 = std::sqrt((g * k3) * tnz3);
  const float som1 = in.s1 * om1;
  const float som2 = in.s2 * om2;
  const float som3 = in.s3 * om3;
  const float dot23 = (k2x * k3x) + (k2y * k3y);
  const float k23x = k2x + k3x;
  const float k23y = k2y + k3y;
  const float k23 = std::sqrt((k23x * k23x) + (k23y * k23y));
  const float tnz23 = tanz(k23 * depth);
  const float csz23 = cosz(k23 * depth);
  const float omsq23 = (g * k23) * tnz23;
  const float dot123 = (k1x * k23x) + (k1y * k23y);
  const float s23 = som2 + som3;
  const float s123 = (som1 + som2) + som3;

  const float di = -(s23 * ((((k2 * k3) * tnz2) * tnz3) - dot23)) +
                   (0.5F * (((som2 * (k3 * k3)) / (csz3 * csz3)) +
                            ((som3 * (k2 * k2)) / (csz2 * csz2))));
  const float e =
      (0.5F / g) * (dot23 - (((som2 * som3) / gsq) *
                             (((om2 * om2) + (om3 * om3)) + (som2 * som3))));
  const float p1 = (2.0F * s123) * ((((om1 * om1) * omsq23) / gsq) - dot123);
  const float p2 = -((som1 * (k23 * k23)) / (csz23 * csz23));
  const float p3 = -((s23 * (k1 * k1)) / (csz1 * csz1));
  float denominator = omsq23 - (s23 * s23);
  if (add_epsilon) {
    denominator = denominator + kEpsilonHh;
  }
  const float t1 = (di / denominator) * ((p1 + p2) + p3);
  const float t2 = -(((di * som1) / gsq) * ((om1 * om1) + omsq23));
  const float p4 = (g * (k1 * k1)) / (csz1 * csz1);
  const float som1_cubed = (som1 * som1) * som1;
  const float t3 = e * ((((som1_cubed * s23) / g) - (g * dot123)) - p4);
  const float t4 =
      (((0.5F * som1) / gsq) * dot23) *
      ((s123 * ((om2 * om2) + (om3 * om3))) + ((som2 * som3) * s23));
  const float t5 = -(((((0.5F * som1) * (om2 * om2)) * (k3 * k3)) / gsq) *
                     ((som1 + som2) + (2.0F * som3))) -
                   (((((0.5F * som1) * (om3 * om3)) * (k2 * k2)) / gsq) *
                    ((som1 + (2.0F * som2)) + som3));
  return HhTerm{
      .coefficient = (((t1 + t2) + t3) + t4) + t5,
      .omega1 = om1,
      .omega2 = om2,
      .omega3 = om3,
  };
}

} // namespace

// --- compute_z_wnumb --------------------------------------------------------
/**
 * @brief Wavenumber from radian frequency and depth (Z_WNUMB).
 * @details Explicit approximation of the linear dispersion relation with a
 * polynomial in y = sigma^2 d / g, evaluated with Horner's scheme in the
 * order of the Fortran expression. Returns -10 for zero depth or frequency.
 * Moved here from xnl4v5_kernel.cpp, unchanged.
 * @param radian_frequency Radian frequency sigma.
 * @param depth Water depth d.
 * @param gravity Acceleration of gravity g.
 * @return The wavenumber.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl Z_WNUMB by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */
float compute_z_wnumb(const float radian_frequency, const float depth,
                      const float gravity) noexcept {
  if (depth <= 0.0F || radian_frequency <= 0.0F) {
    return kInvalidValue;
  }
  const float omega = (radian_frequency * radian_frequency) / gravity;
  const float y = omega * depth;
  float poly = 0.00039F + (y * 0.00011F);
  poly = 0.00171F + (y * poly);
  poly = 0.00654F + (y * poly);
  poly = 0.02174F + (y * poly);
  poly = 0.06320F + (y * poly);
  poly = 0.16084F + (y * poly);
  poly = 0.35550F + (y * poly);
  poly = 0.66667F + (y * poly);
  const float xx = y * (y + (1.0F / (1.0F + (y * poly))));
  const float x = std::sqrt(xx);
  return x / depth;
}

// --- compute_z_cmpcg --------------------------------------------------------
/**
 * @brief Group velocity from radian frequency and depth (Z_CMPCG).
 * @details Linear theory with the wavenumber of Z_WNUMB, and the deep water
 * value for k d > 30. Returns -10 for zero depth or frequency. Moved here
 * from xnl4v5_kernel.cpp, unchanged.
 * @param radian_frequency Radian frequency sigma.
 * @param depth Water depth d.
 * @param gravity Acceleration of gravity g.
 * @return The group velocity.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl Z_CMPCG by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */
float compute_z_cmpcg(const float radian_frequency, const float depth,
                      const float gravity) noexcept {
  const float k = compute_z_wnumb(radian_frequency, depth, gravity);
  if (depth <= 0.0F || radian_frequency <= 0.0F) {
    return kInvalidValue;
  }
  if (depth * k > kCmpcgDeepLimit) {
    return gravity / (2.0F * radian_frequency);
  }
  return (radian_frequency / k) *
         (0.5F + ((depth * k) / std::sinh((2.0F * depth) * k)));
}

// --- compute_x_disper -------------------------------------------------------
/**
 * @brief Radian frequency from wavenumber and depth (X_DISPER).
 * @details sigma = sqrt(g k) for deep water or k d > 20, and
 * sigma = sqrt(g k tanh(k d)) for finite depth.
 * @param wavenumber Wavenumber k.
 * @param depth Water depth d.
 * @param dispersion 1 deep water, 2 finite depth (IQ_DISP).
 * @param gravity Acceleration of gravity g (Q_GRAV).
 * @return The radian frequency, or -1 for an invalid dispersion type.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_DISPER by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_disper(const float wavenumber, const float depth,
                       const int dispersion, const float gravity) noexcept {
  const float kd = wavenumber * depth;
  const int type = kd > kDeepWaterLimit ? 1 : dispersion;
  if (type == 1) {
    return std::sqrt(gravity * wavenumber);
  }
  if (type == 2) {
    return std::sqrt((gravity * wavenumber) * std::tanh(wavenumber * depth));
  }
  return kInvalidLocus;
}

// --- compute_x_cosk ---------------------------------------------------------
/**
 * @brief Cosine of the polar angle of the locus for a wavenumber (X_COSK).
 * @details Deep water: ((q / sqrt(g) + sqrt(k))^4 - k^2 - |p|^2) / (2 k |p|).
 * Finite depth: (k_z^2 - k^2 - |p|^2) / (2 k |p|), with k_z the wavenumber
 * of the radian frequency q + sigma(k). Limited to [-1, 1]. The power x**4
 * gives the same result as gfortran at all optimization levels, (x^2)^2.
 * @param wavenumber Wavenumber k along the symmetry axis of the locus.
 * @param locus Parameters of the locus.
 * @return The cosine.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_COSK by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_cosk(const float wavenumber,
                     const LocusParameters &locus) noexcept {
  const float k = wavenumber;
  const float pmag = locus.difference_magnitude;
  float value = 0.0F;
  if (locus.dispersion == 1) {
    const float qq = locus.frequency_difference / std::sqrt(locus.gravity);
    const float base = qq + std::sqrt(k);
    const float base_sq = base * base;
    value =
        (((base_sq * base_sq) - (k * k)) - (pmag * pmag)) / ((2.0F * k) * pmag);
  } else if (locus.dispersion == 2) {
    const float wk =
        locus.frequency_difference +
        compute_x_disper(k, locus.depth, locus.dispersion, locus.gravity);
    const float kz = compute_z_wnumb(wk, locus.depth, locus.gravity);
    value = (((kz * kz) - (k * k)) - (pmag * pmag)) / ((2.0F * k) * pmag);
  }
  value = std::max(-1.0F, value);
  value = std::min(1.0F, value);
  return value;
}

// --- compute_x_flocus -------------------------------------------------------
/**
 * @brief Locus function in the wavenumber plane (X_FLOCUS).
 * @details z = q + sigma(k2) - sigma(k2 + p). The radian frequencies are
 * stored in double precision and the sum is computed in double precision,
 * as in the Fortran code. Deep water uses sigma = sqrt(g) (|k|^2)^(1/4) via
 * powf(|k|^2, 0.25), as gfortran does.
 * @param kx x-component of the wavenumber k2.
 * @param ky y-component of the wavenumber k2.
 * @param locus Parameters of the locus.
 * @return The locus function, or -1 for an invalid dispersion type.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_FLOCUS by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_flocus(const float kx, const float ky,
                       const LocusParameters &locus) noexcept {
  const float px = locus.difference_x;
  const float py = locus.difference_y;
  const float k4x = kx + px;
  const float k4y = ky + py;
  double w2 = 0.0;
  double w4 = 0.0;
  if (locus.dispersion == 1) {
    w2 = static_cast<double>(locus.sqrt_gravity *
                             std::pow((kx * kx) + (ky * ky), 0.25F));
    w4 = static_cast<double>(locus.sqrt_gravity *
                             std::pow((k4x * k4x) + (k4y * k4y), 0.25F));
  } else if (locus.dispersion == 2) {
    const float k2m = std::sqrt((kx * kx) + (ky * ky));
    const float k4m = std::sqrt((k4x * k4x) + (k4y * k4y));
    w2 = static_cast<double>(
        compute_x_disper(k2m, locus.depth, locus.dispersion, locus.gravity));
    w4 = static_cast<double>(
        compute_x_disper(k4m, locus.depth, locus.dispersion, locus.gravity));
  } else {
    return kInvalidLocus;
  }
  return static_cast<float>(
      (static_cast<double>(locus.frequency_difference) + w2) - w4);
}

// --- compute_x_locus1 -------------------------------------------------------
/**
 * @brief Locus function along the symmetry axis (X_LOCUS1).
 * @details z = q + sigma(k2) - sigma(|k2 - |p||).
 * @param wavenumber Magnitude of the wavenumber k2 on the symmetry axis.
 * @param locus Parameters of the locus.
 * @return The locus function, or -1 for an invalid dispersion type.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_LOCUS1 by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_locus1(const float wavenumber,
                       const LocusParameters &locus) noexcept {
  const float k2 = wavenumber;
  const float k4 = std::abs(-locus.difference_magnitude + k2);
  float w2 = 0.0F;
  float w4 = 0.0F;
  if (locus.dispersion == 1) {
    w2 = locus.sqrt_gravity * std::sqrt(k2);
    w4 = locus.sqrt_gravity * std::sqrt(k4);
  } else if (locus.dispersion == 2) {
    w2 = compute_x_disper(k2, locus.depth, locus.dispersion, locus.gravity);
    w4 = compute_x_disper(k4, locus.depth, locus.dispersion, locus.gravity);
  } else {
    return kInvalidLocus;
  }
  return (locus.frequency_difference + w2) - w4;
}

// --- compute_x_locus2 -------------------------------------------------------
/**
 * @brief Locus function perpendicular to the symmetry axis (X_LOCUS2).
 * @details k2 = k_mid + lambda (-p_y, p_x), k4 = k2 + p, and
 * z = q + sigma(k2) - sigma(k4).
 * @param lambda Distance from the midpoint, in units of |p|.
 * @param locus Parameters of the locus.
 * @return The locus function, or -1 for an invalid dispersion type.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_LOCUS2 by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_locus2(const float lambda,
                       const LocusParameters &locus) noexcept {
  const float kk2x = locus.midpoint_x - (lambda * locus.difference_y);
  const float kk2y = locus.midpoint_y + (lambda * locus.difference_x);
  const float kk2m = std::sqrt((kk2x * kk2x) + (kk2y * kk2y));
  const float kk4x = kk2x + locus.difference_x;
  const float kk4y = kk2y + locus.difference_y;
  const float kk4m = std::sqrt((kk4x * kk4x) + (kk4y * kk4y));
  float w2 = 0.0F;
  float w4 = 0.0F;
  if (locus.dispersion == 1) {
    w2 = locus.sqrt_gravity * std::sqrt(kk2m);
    w4 = locus.sqrt_gravity * std::sqrt(kk4m);
  } else if (locus.dispersion == 2) {
    w2 = compute_x_disper(kk2m, locus.depth, locus.dispersion, locus.gravity);
    w4 = compute_x_disper(kk4m, locus.depth, locus.dispersion, locus.gravity);
  } else {
    return kInvalidLocus;
  }
  return (locus.frequency_difference + w2) - w4;
}

// --- compute_x_jacobian -----------------------------------------------------
/**
 * @brief Jacobian term |c_g2 - c_g4| of the locus integral (X_JACOBIAN).
 * @details Group velocities of linear theory (deep water value for
 * k d > 20), with the directions from atan2(x, y) as in the Fortran code.
 * @param k2x x-component of the wavenumber k2.
 * @param k2y y-component of the wavenumber k2.
 * @param k4x x-component of the wavenumber k4.
 * @param k4y y-component of the wavenumber k4.
 * @param depth Water depth (Q_DEPTH).
 * @param gravity Acceleration of gravity (Q_GRAV).
 * @return The Jacobian term.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_JACOBIAN by G. Ph. van
 * Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_jacobian(const float k2x, const float k2y, const float k4x,
                         const float k4y, const float depth,
                         const float gravity) noexcept {
  const float k2m = std::sqrt((k2x * k2x) + (k2y * k2y));
  const float k4m = std::sqrt((k4x * k4x) + (k4y * k4y));
  const float ang2 = std::atan2(k2x, k2y);
  const float ang4 = std::atan2(k4x, k4y);
  const float sig2 = std::sqrt((gravity * k2m) * std::tanh(k2m * depth));
  const float sig4 = std::sqrt((gravity * k4m) * std::tanh(k4m * depth));
  const float k2md = k2m * depth;
  const float k4md = k4m * depth;
  const float cg2 =
      k2md > kDeepWaterLimit
          ? (0.5F * gravity) / sig2
          : (sig2 / k2m) * (0.5F + (k2md / std::sinh(2.0F * k2md)));
  const float cg4 =
      k4md > kDeepWaterLimit
          ? (0.5F * gravity) / sig4
          : (sig4 / k4m) * (0.5F + (k4md / std::sinh(2.0F * k4md)));
  return std::sqrt(((cg2 * cg2) + (cg4 * cg4)) -
                   (((2.0F * cg2) * cg4) * std::cos(ang2 - ang4)));
}

// --- compute_xc_webb --------------------------------------------------------
/**
 * @brief Deep water coupling coefficient of Webb (1978) (XC_WEBB).
 * @details The mixed precision of the Fortran code is kept: the squared sums
 * of frequencies WSQ* are double precision variables computed from single
 * precision expressions, and the expressions with them are evaluated in
 * double precision and rounded to single precision on assignment.
 * @param k1x x-component of the wavenumber k1.
 * @param k1y y-component of the wavenumber k1.
 * @param k2x x-component of the wavenumber k2.
 * @param k2y y-component of the wavenumber k2.
 * @param k3x x-component of the wavenumber k3.
 * @param k3y y-component of the wavenumber k3.
 * @param k4x x-component of the wavenumber k4.
 * @param k4y y-component of the wavenumber k4.
 * @param gravity Acceleration of gravity.
 * @return The coupling coefficient.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl XC_WEBB by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_xc_webb(const float k1x, const float k1y, const float k2x,
                      const float k2y, const float k3x, const float k3y,
                      const float k4x, const float k4y,
                      const float gravity) noexcept {
  const float pi = 4.0F * std::atan(1.0F);
  const float pi4 = 0.25F * pi;
  const float k1 = std::sqrt((k1x * k1x) + (k1y * k1y));
  const float k2 = std::sqrt((k2x * k2x) + (k2y * k2y));
  const float k3 = std::sqrt((k3x * k3x) + (k3y * k3y));
  const float k4 = std::sqrt((k4x * k4x) + (k4y * k4y));
  const float w1 = std::sqrt(k1);
  const float w2 = std::sqrt(k2);
  const float w3 = std::sqrt(k3);
  const float w4 = std::sqrt(k4);
  const float dot12 = (k1x * k2x) + (k1y * k2y);
  const float dot13 = (k1x * k3x) + (k1y * k3y);
  const float dot14 = (k1x * k4x) + (k1y * k4y);
  const float dot23 = (k2x * k3x) + (k2y * k3y);
  const float dot24 = (k2x * k4x) + (k2y * k4y);
  const float dot34 = (k3x * k4x) + (k3y * k4y);
  const auto magnitude = [](const float x, const float y) {
    return static_cast<double>(std::sqrt((x * x) + (y * y)));
  };
  const auto square = [](const float x) { return static_cast<double>(x * x); };
  const double wsqp12 = magnitude(k1x + k2x, k1y + k2y);
  const double wsq12 = square(w1 + w2);
  const double wsqm13 = magnitude(k1x - k3x, k1y - k3y);
  const double wsq13 = square(w1 - w3);
  const double wsqm14 = magnitude(k1x - k4x, k1y - k4y);
  const double wsq14 = square(w1 - w4);
  const auto z12 = static_cast<float>(wsqp12 - wsq12);
  const auto z13 = static_cast<float>(wsqm13 - wsq13);
  const auto z14 = static_cast<float>(wsqm14 - wsq14);
  auto z = static_cast<float>(
      ((2.0 * wsq12) * static_cast<double>((k1 * k2) - dot12)) *
      static_cast<double>((k3 * k4) - dot34));
  const float p1 = z / (z12 + kEpsilonWebb);
  z = static_cast<float>(
      ((2.0 * wsq13) * static_cast<double>((k1 * k3) + dot13)) *
      static_cast<double>((k2 * k4) + dot24));
  const float p2 = z / (z13 + kEpsilonWebb);
  z = static_cast<float>(
      ((2.0 * wsq14) * static_cast<double>((k1 * k4) + dot14)) *
      static_cast<double>((k2 * k3) + dot23));
  const float p3 = z / (z14 + kEpsilonWebb);
  const float p4 =
      0.5F * (((dot12 * dot34) + (dot13 * dot24)) + (dot14 * dot23));
  const auto p5 = static_cast<float>(
      (static_cast<double>(0.25F * (dot13 + dot24)) * wsq13) * wsq13);
  const auto p6 = static_cast<float>(
      -((static_cast<double>(0.25F * (dot12 + dot34)) * wsq12) * wsq12));
  const auto p7 = static_cast<float>(
      (static_cast<double>(0.25F * (dot14 + dot23)) * wsq14) * wsq14);
  const float p8 = (((2.5F * k1) * k2) * k3) * k4;
  const auto p9 = static_cast<float>(
      ((wsq12 * wsq13) * wsq14) * static_cast<double>(((k1 + k2) + k3) + k4));
  const float dwebb = (((((((p1 + p2) + p3) + p4) + p5) + p6) + p7) + p8) + p9;
  return ((((gravity * gravity) * pi4) * dwebb) * dwebb) /
         ((((w1 * w2) * w3) * w4) + kEpsilonWebb);
}

// --- compute_xc_hh ----------------------------------------------------------
/**
 * @brief Finite depth coupling coefficient of Herterich and Hasselmann (1980)
 * (XC_HH).
 * @details Sum of three partial coefficients for cyclic permutations of the
 * three wavenumbers with signs (+, +, -), (+, -, +) and (-, +, +) (hh_term),
 * squared and scaled with the frequencies of the last permutation, as in the
 * Fortran code. XC_HH uses g = 9.81 internally. The fourth wavenumber of the
 * Fortran argument list is not used and is omitted.
 * @param w1x x-component of the first wavenumber.
 * @param w1y y-component of the first wavenumber.
 * @param w2x x-component of the second wavenumber.
 * @param w2y y-component of the second wavenumber.
 * @param w3x x-component of the third wavenumber.
 * @param w3y y-component of the third wavenumber.
 * @param depth Water depth h.
 * @return The coupling coefficient.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl XC_HH by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_xc_hh(const float w1x, const float w1y, const float w2x,
                    const float w2y, const float w3x, const float w3y,
                    const float depth) noexcept {
  const float gsq = kGravityHh * kGravityHh;
  const HhTerm term1 = hh_term(
      HhInput{
          .k1x0 = w1x,
          .k1y0 = w1y,
          .k2x0 = w2x,
          .k2y0 = w2y,
          .k3x0 = w3x,
          .k3y0 = w3y,
          .s1 = 1.0F,
          .s2 = 1.0F,
          .s3 = -1.0F,
      },
      depth, true);
  const HhTerm term2 = hh_term(
      HhInput{
          .k1x0 = w2x,
          .k1y0 = w2y,
          .k2x0 = w3x,
          .k2y0 = w3y,
          .k3x0 = w1x,
          .k3y0 = w1y,
          .s1 = 1.0F,
          .s2 = -1.0F,
          .s3 = 1.0F,
      },
      depth, false);
  const HhTerm term3 = hh_term(
      HhInput{
          .k1x0 = w3x,
          .k1y0 = w3y,
          .k2x0 = w1x,
          .k2y0 = w1y,
          .k3x0 = w2x,
          .k3y0 = w2y,
          .s1 = -1.0F,
          .s2 = 1.0F,
          .s3 = 1.0F,
      },
      depth, false);
  const float stot =
      (term1.coefficient + term2.coefficient) + term3.coefficient;
  const float om1 = term3.omega1;
  const float om2 = term3.omega2;
  const float om3 = term3.omega3;
  const float om4 = (om2 + om3) - om1;
  const float dsq = (((stot * stot) * kQuarterPiHh) * gsq) /
                    ((((om1 * om2) * om3) * om4) + kEpsilonHh);
  return dsq * gsq;
}

// --- compute_x_cple ---------------------------------------------------------
/**
 * @brief Coupling coefficient for a quadruplet (X_CPLE).
 * @details IQ_CPLE = 1: XC_WEBB with (k1, k2, k3, k4). IQ_CPLE = 2: XC_HH
 * with (k4, k3, k2), as in the Fortran code. Zero for IQ_CPLE outside
 * [1, 4]. For IQ_CPLE = 3 and 4 the Fortran code returns an undefined
 * value; the C++ code returns zero.
 * @param k1x x-component of the wavenumber k1.
 * @param k1y y-component of the wavenumber k1.
 * @param k2x x-component of the wavenumber k2.
 * @param k2y y-component of the wavenumber k2.
 * @param k3x x-component of the wavenumber k3.
 * @param k3y y-component of the wavenumber k3.
 * @param k4x x-component of the wavenumber k4.
 * @param k4y y-component of the wavenumber k4.
 * @param coupling 1 Webb, 2 Herterich and Hasselmann (IQ_CPLE).
 * @param depth Water depth.
 * @param gravity Acceleration of gravity.
 * @return The coupling coefficient.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Translated from Xnl X_CPLE by G. Ph. van Vledder.
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float compute_x_cple(const float k1x, const float k1y, const float k2x,
                     const float k2y, const float k3x, const float k3y,
                     const float k4x, const float k4y, const int coupling,
                     const float depth, const float gravity) noexcept {
  if (coupling == 1) {
    return compute_xc_webb(k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, gravity);
  }
  if (coupling == 2) {
    // XC_HH is called with (k4, k3, k2), as in the Fortran code
    return compute_xc_hh(k4x, k4y, k3x, k3y, k2x, k2y, // NOLINT
                         depth);
  }
  return 0.0F;
}

} // namespace xnl4v5

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- x_disper_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_disper (X_DISPER).
 * @param k Wavenumber.
 * @param d Water depth.
 * @param iq_disp Dispersion type.
 * @param q_grav Acceleration of gravity.
 * @return The radian frequency.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_disper_cpp(const float k, const float d, const int iq_disp,
                   const float q_grav) noexcept {
  return xnl4v5::compute_x_disper(k, d, iq_disp, q_grav);
}

// --- x_cosk_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_cosk (X_COSK).
 * @param k Wavenumber along the symmetry axis.
 * @param q Difference of radian frequencies.
 * @param pmag Magnitude of k1 - k3.
 * @param iq_disp Dispersion type.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The cosine of the polar angle.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_cosk_cpp(const float k, const float q, const float pmag,
                 const int iq_disp, const float q_depth,
                 const float q_grav) noexcept {
  const xnl4v5::LocusParameters locus{
      .frequency_difference = q,
      .difference_x = 0.0F,
      .difference_y = 0.0F,
      .difference_magnitude = pmag,
      .midpoint_x = 0.0F,
      .midpoint_y = 0.0F,
      .depth = q_depth,
      .gravity = q_grav,
      .sqrt_gravity = 0.0F,
      .dispersion = iq_disp,
  };
  return xnl4v5::compute_x_cosk(k, locus);
}

// --- x_flocus_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_flocus (X_FLOCUS).
 * @param kxx x-component of k2.
 * @param kyy y-component of k2.
 * @param q Difference of radian frequencies.
 * @param px x-component of k1 - k3.
 * @param py y-component of k1 - k3.
 * @param sqrtg sqrt(g) of module m_constants.
 * @param iq_disp Dispersion type.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The locus function.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_flocus_cpp(const float kxx, const float kyy, const float q,
                   const float px, const float py, const float sqrtg,
                   const int iq_disp, const float q_depth,
                   const float q_grav) noexcept {
  const xnl4v5::LocusParameters locus{
      .frequency_difference = q,
      .difference_x = px,
      .difference_y = py,
      .difference_magnitude = 0.0F,
      .midpoint_x = 0.0F,
      .midpoint_y = 0.0F,
      .depth = q_depth,
      .gravity = q_grav,
      .sqrt_gravity = sqrtg,
      .dispersion = iq_disp,
  };
  return xnl4v5::compute_x_flocus(kxx, kyy, locus);
}

// --- x_locus1_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_locus1 (X_LOCUS1).
 * @param k2 Magnitude of k2.
 * @param q Difference of radian frequencies.
 * @param pmag Magnitude of k1 - k3.
 * @param sqrtg sqrt(g) of module m_constants.
 * @param iq_disp Dispersion type.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The locus function.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_locus1_cpp(const float k2, const float q, const float pmag,
                   const float sqrtg, const int iq_disp, const float q_depth,
                   const float q_grav) noexcept {
  const xnl4v5::LocusParameters locus{
      .frequency_difference = q,
      .difference_x = 0.0F,
      .difference_y = 0.0F,
      .difference_magnitude = pmag,
      .midpoint_x = 0.0F,
      .midpoint_y = 0.0F,
      .depth = q_depth,
      .gravity = q_grav,
      .sqrt_gravity = sqrtg,
      .dispersion = iq_disp,
  };
  return xnl4v5::compute_x_locus1(k2, locus);
}

// --- x_locus2_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_locus2 (X_LOCUS2).
 * @param lambda Distance from the midpoint, in units of |p|.
 * @param q Difference of radian frequencies.
 * @param px x-component of k1 - k3.
 * @param py y-component of k1 - k3.
 * @param kmidx x-component of the locus midpoint.
 * @param kmidy y-component of the locus midpoint.
 * @param sqrtg sqrt(g) of module m_constants.
 * @param iq_disp Dispersion type.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The locus function.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_locus2_cpp(const float lambda, const float q, const float px,
                   const float py, const float kmidx, const float kmidy,
                   const float sqrtg, const int iq_disp, const float q_depth,
                   const float q_grav) noexcept {
  const xnl4v5::LocusParameters locus{
      .frequency_difference = q,
      .difference_x = px,
      .difference_y = py,
      .difference_magnitude = 0.0F,
      .midpoint_x = kmidx,
      .midpoint_y = kmidy,
      .depth = q_depth,
      .gravity = q_grav,
      .sqrt_gravity = sqrtg,
      .dispersion = iq_disp,
  };
  return xnl4v5::compute_x_locus2(lambda, locus);
}

// --- x_jacobian_cpp ---------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_jacobian (X_JACOBIAN).
 * @param x2 x-component of k2.
 * @param y2 y-component of k2.
 * @param x4 x-component of k4.
 * @param y4 y-component of k4.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The Jacobian term.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_jacobian_cpp(const float x2, const float y2, const float x4,
                     const float y4, const float q_depth,
                     const float q_grav) noexcept {
  return xnl4v5::compute_x_jacobian(x2, y2, x4, y4, q_depth, q_grav);
}

// --- x_cple_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_x_cple (X_CPLE).
 * @param k1x x-component of k1.
 * @param k1y y-component of k1.
 * @param k2x x-component of k2.
 * @param k2y y-component of k2.
 * @param k3x x-component of k3.
 * @param k3y y-component of k3.
 * @param k4x x-component of k4.
 * @param k4y y-component of k4.
 * @param iq_cple Type of coupling coefficient.
 * @param depth Water depth.
 * @param grav Acceleration of gravity.
 * @return The coupling coefficient.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
float x_cple_cpp(const float k1x, const float k1y, const float k2x,
                 const float k2y, const float k3x, const float k3y,
                 const float k4x, const float k4y, const int iq_cple,
                 const float depth, const float grav) noexcept {
  return xnl4v5::compute_x_cple(k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, iq_cple,
                                depth, grav);
}

} // extern "C"
