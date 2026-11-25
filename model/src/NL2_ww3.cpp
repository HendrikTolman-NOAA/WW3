/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file NL2_ww3.cpp
 * @brief Implementation of the single-point interaction kernel of the NL2
 * (WRT, exact) nonlinear four-wave interaction source term for WAVEWATCH IV
 * (WW4).
 * @details Pure C++20 implementation of the kernel of the Webb-Resio-Tracy
 * method of the Xnl library (Q_XNL4V4, Q_T13V4, Q_GETLOCUS, Z_CMPCG and
 * Z_WNUMB), with a C-compatible interface called from Q_XNL4V4 in
 * mod_xnl4v5.f90. First stage of the WW4 conversion of the NL2 source term.
 * The method, the input and output spectrum and what remains in Fortran are
 * described in NL2_ww3.h. Each routine has its full documentation (details,
 * contracts and heritage) in front of its definition, and a shorter version
 * with its declaration in NL2_ww3.h.
 *
 * References:
 * - Hasselmann, K., 1962: On the non-linear energy transfer in a gravity wave
 *   spectrum, Part 1. General theory. J. Fluid Mech., 12, 481-500.
 * - Hasselmann, K., 1963: On the non-linear energy transfer in a gravity wave
 *   spectrum, Part 2. Conservation theory, wave-particle correspondence,
 *   irreversibility. J. Fluid Mech., 15, 273-281.
 * - Hasselmann, K., 1963: On the non-linear energy transfer in a gravity wave
 *   spectrum, Part 3. Evaluation of energy flux and sea-swell interactions
 *   for a Neuman spectrum. J. Fluid Mech., 15, 385-398.
 * - Webb, D. J., 1978: Non-linear transfers between sea waves. Deep-Sea
 *   Res., 25, 279-298.
 * - Tracy, B. and D. T. Resio, 1982: Theory and calculation of the nonlinear
 *   energy transfer between sea waves in deep water. WES Report 11, US Army
 *   Corps of Engineers.
 * - Resio, D. T. and W. Perrie, 1991: A numerical study of nonlinear energy
 *   fluxes due to wave-wave interactions. Part 1: Methodology and basic
 *   results. J. Fluid Mech., 223, 603-629.
 * - van Vledder, G. Ph., 2006: The WRT method for the computation of
 *   non-linear four-wave interactions in discrete spectral wave models.
 *   Coastal Eng., 53, 223-242.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section NL2.
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
 * @author Code Heritage: Converted and adapted from the Xnl library
 * mod_xnl4v5.f90 (Q_XNL4V4, Q_T13V4, Q_GETLOCUS) and serv_xnl4v5.f90
 * (Z_CMPCG, Z_WNUMB) by Gerbrant Ph. van Vledder, as used in WAVEWATCH III
 * (WW3) w3snl2md.F90 by Hendrik L. Tolman and G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */

// Local
#include "NL2_ww3.h"

// System
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <span>
#include <vector>

namespace ww4_nl2 {

namespace {

// --- Constants --------------------------------------------------------------
/// Value returned by Z_WNUMB and Z_CMPCG for zero depth or frequency.
constexpr float kInvalidValue = -10.0F;

/// k d above which Z_CMPCG uses the deep water group velocity.
constexpr float kDeepWaterLimit = 30.0F;

/// Mirror offset in the direction indices of transformed loci (IMIRROR).
constexpr int kMirrorOffset = 1;

/// Offset added before truncating the reference direction index.
constexpr float kDirectionRounding = 0.01F;

/// Half circle in degrees, used for the directional difference of k1 and k3.
constexpr float kHalfCircleDeg = 180.0F;

// --- fortran_maxval ---------------------------------------------------------
/**
 * @brief Maximum of an array as Fortran MAXVAL in gfortran.
 * @details Starts from -infinity, takes the first value that is not NaN, and
 * then the maximum of the remaining values. Returns NaN if all values are
 * NaN, and -infinity for an empty array.
 * @param values The array.
 * @return The maximum value.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
float fortran_maxval(const std::span<const float> values) noexcept {
  float result = -std::numeric_limits<float>::infinity();
  std::size_t i = 0;
  for (; i < values.size(); ++i) {
    if (values[i] >= result) {
      result = values[i];
      break;
    }
  }
  if (i == values.size()) {
    return values.empty() ? result : std::numeric_limits<float>::quiet_NaN();
  }
  for (++i; i < values.size(); ++i) {
    // std::max(a, b) returns b only if a < b, as the test b > a
    result = std::max(result, values[i]);
  }
  return result;
}

// --- sixth_power ------------------------------------------------------------
/**
 * @brief Sixth power x**6 as evaluated by gfortran at -O3.
 * @details Computed as x**3 * x**3 with x**3 = x**2 * x, the multiplication
 * chain gfortran generates for a constant exponent at -O2 and -O3 (the WW3
 * release optimization). At -O0, gfortran calls __powisf2, which computes
 * x**2 * x**4 and can differ in the last bit.
 * @param value The value x.
 * @return The sixth power x**6.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
constexpr float sixth_power(const float value) noexcept {
  const float power2 = value * value;
  const float power3 = power2 * value;
  return power3 * power3;
}

// --- compute_z_wnumb --------------------------------------------------------
/**
 * @brief Wavenumber from radian frequency and depth (Z_WNUMB).
 * @details Explicit approximation of the linear dispersion relation with a
 * polynomial in y = sigma^2 d / g, evaluated with Horner's scheme in the
 * order of the Fortran expression. Returns -10 for zero depth or frequency.
 * @param radian_frequency Radian frequency sigma.
 * @param depth Water depth d.
 * @param gravity Acceleration of gravity g.
 * @return The wavenumber.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Z_WNUMB by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
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
 * value for k d > 30. Returns -10 for zero depth or frequency.
 * @param radian_frequency Radian frequency sigma.
 * @param depth Water depth d.
 * @param gravity Acceleration of gravity g.
 * @return The group velocity.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Z_CMPCG by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
float compute_z_cmpcg(const float radian_frequency, const float depth,
                      const float gravity) noexcept {
  const float k = compute_z_wnumb(radian_frequency, depth, gravity);
  if (depth <= 0.0F || radian_frequency <= 0.0F) {
    return kInvalidValue;
  }
  if (depth * k > kDeepWaterLimit) {
    return gravity / (2.0F * radian_frequency);
  }
  return (radian_frequency / k) *
         (0.5F + ((depth * k) / std::sinh((2.0F * depth) * k)));
}

// --- TransformedLocus -------------------------------------------------------
/**
 * @brief Locus of one (k1, k3) combination, transformed from the reference
 * locus (T_* arrays of the Xnl library).
 * @details Indices are 1-based, as in the Fortran code.
 */
struct TransformedLocus {
  int num_points = 0;            ///< Number of points (NLOCUSX).
  std::vector<int> k2_index;     ///< k index of k2 (T_IK2).
  std::vector<int> a2_index;     ///< Direction index of k2 (T_IA2).
  std::vector<int> k4_index;     ///< k index of k4 (T_IK4).
  std::vector<int> a4_index;     ///< Direction index of k4 (T_IA4).
  std::vector<float> k2_weight1; ///< Weight 1 of k2 (T_W1K2).
  std::vector<float> k2_weight2; ///< Weight 2 of k2 (T_W2K2).
  std::vector<float> k2_weight3; ///< Weight 3 of k2 (T_W3K2).
  std::vector<float> k2_weight4; ///< Weight 4 of k2 (T_W4K2).
  std::vector<float> k4_weight1; ///< Weight 1 of k4 (T_W1K4).
  std::vector<float> k4_weight2; ///< Weight 2 of k4 (T_W2K4).
  std::vector<float> k4_weight3; ///< Weight 3 of k4 (T_W3K4).
  std::vector<float> k4_weight4; ///< Weight 4 of k4 (T_W4K4).
  std::vector<float> factor;     ///< Factor Z along the locus (T_ZZ).

  /**
   * @brief Allocates the arrays for a maximum number of points.
   * @param size Maximum number of points on a locus (KLOCUS).
   */
  explicit TransformedLocus(const std::size_t size)
      : k2_index(size), a2_index(size), k4_index(size), a4_index(size),
        k2_weight1(size), k2_weight2(size), k2_weight3(size), k2_weight4(size),
        k4_weight1(size), k4_weight2(size), k4_weight3(size), k4_weight4(size),
        factor(size) {}
};

// --- LocusResult ------------------------------------------------------------
/**
 * @brief Result of the locus lookup (Q_GETLOCUS).
 */
enum class LocusResult : std::uint8_t {
  kFound = 0,          ///< Locus found and transformed.
  kEmpty = 1,          ///< Locus has no points (T13 is zero).
  kInvalidAddress = 2, ///< Direction entry above IAMAX.
};

// --- KernelContext ----------------------------------------------------------
/**
 * @brief Grid, loci and spectrum shared by the kernel routines.
 * @details Holds copies of the (non-owning) grid, locus and settings views.
 */
struct KernelContext {
  std::size_t num_wavenumbers;    ///< Number of wavenumbers (NKQ).
  std::size_t num_directions;     ///< Number of directions (NAQ).
  WavenumberGrid grid;            ///< Wavenumber and direction grid.
  InteractionGrid interactions;   ///< Precomputed loci.
  KernelSettings settings;        ///< Settings.
  std::span<const float> density; ///< Spectrum N(k,theta) (NSPEC).

  /**
   * @brief 0-based index in a (k, theta) array of 1-based indices.
   * @param ik 1-based wavenumber index.
   * @param ia 1-based direction index.
   * @return The index ik - 1 + nk (ia - 1).
   */
  [[nodiscard]] std::size_t at(const int ik, const int ia) const noexcept {
    return static_cast<std::size_t>(ik - 1) +
           (num_wavenumbers * static_cast<std::size_t>(ia - 1));
  }

  /**
   * @brief Wavenumber of a 1-based index.
   * @param ik 1-based wavenumber index.
   * @return Q_K(ik).
   */
  [[nodiscard]] float k_at(const int ik) const noexcept {
    return grid.wavenumber[static_cast<std::size_t>(ik - 1)];
  }

  /**
   * @brief Wavenumber step of a 1-based index.
   * @param ik 1-based wavenumber index.
   * @return Q_DK(ik).
   */
  [[nodiscard]] float dk_at(const int ik) const noexcept {
    return grid.wavenumber_step[static_cast<std::size_t>(ik - 1)];
  }

  /**
   * @brief Direction in degrees of a 1-based index.
   * @param ia 1-based direction index.
   * @return Q_AD(ia).
   */
  [[nodiscard]] float angle_at(const int ia) const noexcept {
    return grid.direction_deg[static_cast<std::size_t>(ia - 1)];
  }
};

// --- direction_index --------------------------------------------------------
/**
 * @brief Direction index of a mirrored locus point.
 * @details ibeta + 2 iaref - ia - imirror, evaluated in single precision and
 * truncated as in Q_GETLOCUS (exact for the small integers involved).
 * @param ibeta Rotation index.
 * @param reference_direction Reference direction IAREF.
 * @param index Direction index of the reference locus point.
 * @return The transformed direction index.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
int direction_index(const int ibeta, const int reference_direction,
                    const int index) noexcept {
  const float value = ((static_cast<float>(ibeta) +
                        (2.0F * static_cast<float>(reference_direction))) -
                       static_cast<float>(index)) -
                      static_cast<float>(kMirrorOffset);
  return static_cast<int>(value);
}

// --- reference_rotation -----------------------------------------------------
/**
 * @brief Rotation index int(theta / dtheta + 0.01) of a direction.
 * @param context Kernel context.
 * @param direction 1-based direction index.
 * @return The rotation index.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
int reference_rotation(const KernelContext &context,
                       const int direction) noexcept {
  return static_cast<int>(
      (context.angle_at(direction) / context.grid.direction_step_deg) +
      kDirectionRounding);
}

// --- LocusEntry -------------------------------------------------------------
/**
 * @brief Entry of the reference locus of a (k1, k3) combination.
 */
struct LocusEntry {
  int k_entry; ///< 1-based wavenumber entry KMEM.
  int a_entry; ///< 1-based direction entry AMEM.
  int it1;     ///< Direction index of k1, shifted for circle grids.
  int it3;     ///< Direction index of k3, shifted for circle grids.
  int itmin;   ///< Smallest of it1 and it3.
  int ikmin;   ///< Smallest wavenumber index of k1 and k3.
};

// --- locus_entry ------------------------------------------------------------
/**
 * @brief Determines the reference locus entry of a (k1, k3) combination.
 * @details First part of Q_GETLOCUS: the entries KMEM and AMEM from the
 * differences of the wavenumber and direction indices of k1 and k3, with the
 * shortest directional difference on circle grids (IQ_GRID = 3).
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param ik3 1-based wavenumber index of k3.
 * @param ia3 1-based direction index of k3.
 * @return The locus entry and the (shifted) direction indices.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_GETLOCUS by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
LocusEntry locus_entry(const KernelContext &context, const int ik1,
                       const int ia1, const int ik3, const int ia3) noexcept {
  const KernelSettings &settings = context.settings;
  const auto nkq = static_cast<int>(context.num_wavenumbers);
  const auto naq = static_cast<int>(context.num_directions);
  LocusEntry entry{
      .k_entry = 0,
      .a_entry = 0,
      .it1 = ia1,
      .it3 = ia3,
      .itmin = std::min(ia1, ia3),
      .ikmin = std::min(ik1, ik3),
  };
  const int ikdif = std::abs(ik1 - ik3);
  int iadif = std::abs(ia1 - ia3);
  if (settings.grid_type == 3) {
    const int nhalf = naq / 2;
    if (iadif > nhalf) {
      if (entry.it1 > nhalf) {
        entry.it1 = entry.it1 - naq;
      }
      if (entry.it3 > nhalf) {
        entry.it3 = entry.it3 - naq;
      }
    }
    entry.itmin = std::min(entry.it1, entry.it3);
    iadif = (naq - std::abs(naq - (2 * std::abs(entry.it1 - entry.it3)))) / 2;
  }
  entry.a_entry = iadif + 1;
  if (settings.geometric_grid == 0) {
    const int jk1 = std::min(ik1, ik3);
    const int jk3 = std::max(ik1, ik3);
    entry.k_entry = (jk3 - jk1 + 1) - (((jk1 - (2 * nkq) - 2) * (jk1 - 1)) / 2);
  } else {
    entry.k_entry = ikdif + 1;
  }
  return entry;
}

// --- LocusTransform ---------------------------------------------------------
/**
 * @brief Transformation of a reference locus to an actual (k1, k3).
 */
struct LocusTransform {
  int k_offset;        ///< Shift of the wavenumber indices (KDIF).
  int ibeta;           ///< Rotation of the direction indices.
  bool swap_k;         ///< Wavenumber indices of k2 and k4 exchanged.
  bool mirror;         ///< Directions mirrored.
  bool swap_direction; ///< Direction indices of k2 and k4 exchanged.
  bool swap_weights;   ///< Weights exchanged for mirrored directions.
};

// --- select_transform -------------------------------------------------------
/**
 * @brief Selects the transformation of the reference locus (six cases of
 * Q_GETLOCUS).
 * @details Case 1 (k3 > k1, theta3 >= theta1) and case 5 (k3 = k1,
 * theta3 > theta1) shift and rotate. Case 4 (k1 > k3, theta1 > theta3) also
 * exchanges k2 and k4. Cases 2, 3 and 6 mirror the directions and exchange
 * the weights. Case 3 (k1 > k3, theta3 >= theta1) exchanges k2 and k4 in the
 * wavenumber and direction indices, and case 6 (k1 = k3, theta1 > theta3)
 * only in the wavenumber indices, as in the Fortran code.
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param ik3 1-based wavenumber index of k3.
 * @param ia3 1-based direction index of k3.
 * @param entry Locus entry and shifted direction indices.
 * @param k_offset Shift of the wavenumber indices.
 * @param[out] transform The selected transformation.
 * @return False if none of the six cases applies (k1 = k3).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_GETLOCUS by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
bool select_transform(const KernelContext &context, const int ik1,
                      const int ia1, const int ik3, const int ia3,
                      const LocusEntry &entry, const int k_offset,
                      LocusTransform &transform) noexcept {
  const int it1 = entry.it1;
  const int it3 = entry.it3;
  const int shift = entry.itmin - context.settings.reference_direction;
  // Transformation with flags (swap_k, mirror, swap_direction, swap_weights)
  const auto make = [k_offset](const int ibeta, const bool swap_k,
                               const bool mirror, const bool swap_direction,
                               const bool swap_weights) {
    return LocusTransform{
        .k_offset = k_offset,
        .ibeta = ibeta,
        .swap_k = swap_k,
        .mirror = mirror,
        .swap_direction = swap_direction,
        .swap_weights = swap_weights,
    };
  };
  if ((ik3 > ik1 && it3 >= it1) || (ik1 == ik3 && it3 > it1)) { // Case 1, 5
    transform = make(shift, false, false, false, false);
  } else if (ik3 > ik1 && it3 < it1) { // Case 2
    transform =
        make(reference_rotation(context, ia1), false, true, false, true);
  } else if (ik1 > ik3 && it3 >= it1) { // Case 3
    transform = make(reference_rotation(context, ia3), true, true, true, true);
  } else if (ik1 > ik3 && it1 > it3) { // Case 4
    transform = make(shift, true, false, true, false);
  } else if (ik1 == ik3 && it1 > it3) { // Case 6
    transform = make(reference_rotation(context, ia1), true, true, false, true);
  } else {
    return false;
  }
  return true;
}

// --- copy_weights -----------------------------------------------------------
/**
 * @brief Copies the interpolation weights of one locus point.
 * @details With swap, the weights of the two direction neighbours are
 * exchanged (1 <-> 3, 2 <-> 4) for both k2 and k4, as for mirrored loci in
 * Q_GETLOCUS.
 * @param quad Precomputed loci.
 * @param i Index of the reference locus point in the QUAD arrays.
 * @param ui Index of the point in the transformed locus.
 * @param swap Exchange the weights of the direction neighbours.
 * @param[out] locus Transformed locus.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_GETLOCUS by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
void copy_weights(const InteractionGrid &quad, const std::size_t i,
                  const std::size_t ui, const bool swap,
                  TransformedLocus &locus) noexcept {
  if (swap) {
    locus.k2_weight1[ui] = quad.k2_weight3[i];
    locus.k2_weight2[ui] = quad.k2_weight4[i];
    locus.k2_weight3[ui] = quad.k2_weight1[i];
    locus.k2_weight4[ui] = quad.k2_weight2[i];
    locus.k4_weight1[ui] = quad.k4_weight3[i];
    locus.k4_weight2[ui] = quad.k4_weight4[i];
    locus.k4_weight3[ui] = quad.k4_weight1[i];
    locus.k4_weight4[ui] = quad.k4_weight2[i];
  } else {
    locus.k2_weight1[ui] = quad.k2_weight1[i];
    locus.k2_weight2[ui] = quad.k2_weight2[i];
    locus.k2_weight3[ui] = quad.k2_weight3[i];
    locus.k2_weight4[ui] = quad.k2_weight4[i];
    locus.k4_weight1[ui] = quad.k4_weight1[i];
    locus.k4_weight2[ui] = quad.k4_weight2[i];
    locus.k4_weight3[ui] = quad.k4_weight3[i];
    locus.k4_weight4[ui] = quad.k4_weight4[i];
  }
}

// --- transform_locus --------------------------------------------------------
/**
 * @brief Copies a reference locus with a transformation applied.
 * @param context Kernel context.
 * @param entry 0-based index of the locus entry (KMEM, AMEM) in QUAD_NLOC.
 * @param num_points Number of points on the locus.
 * @param transform The transformation.
 * @param[out] locus Transformed locus (indices, directions and weights).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_GETLOCUS by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
void transform_locus(const KernelContext &context, const std::size_t entry,
                     const int num_points, const LocusTransform &transform,
                     TransformedLocus &locus) noexcept {
  const InteractionGrid &quad = context.interactions;
  const int iaref = context.settings.reference_direction;
  const std::size_t stride = quad.num_k_entries * quad.num_a_entries;
  const auto n_points = static_cast<std::size_t>(num_points);
  for (std::size_t ui = 0; ui < n_points; ++ui) {
    const std::size_t i = entry + (stride * ui);
    const int k2 = transform.swap_k ? quad.k4_index[i] : quad.k2_index[i];
    const int k4 = transform.swap_k ? quad.k2_index[i] : quad.k4_index[i];
    locus.k2_index[ui] = transform.k_offset + k2;
    locus.k4_index[ui] = transform.k_offset + k4;
    const int a2 =
        transform.swap_direction ? quad.a4_index[i] : quad.a2_index[i];
    const int a4 =
        transform.swap_direction ? quad.a2_index[i] : quad.a4_index[i];
    if (transform.mirror) {
      locus.a2_index[ui] = direction_index(transform.ibeta, iaref, a2);
      locus.a4_index[ui] = direction_index(transform.ibeta, iaref, a4);
    } else {
      locus.a2_index[ui] = a2 + transform.ibeta;
      locus.a4_index[ui] = a4 + transform.ibeta;
    }
    copy_weights(quad, i, ui, transform.swap_weights, locus);
  }
}

// --- get_locus --------------------------------------------------------------
/**
 * @brief Retrieves and transforms the locus of a (k1, k3) combination
 * (Q_GETLOCUS).
 * @details Determines the entry of the reference locus from the differences
 * of the wavenumber and direction indices of k1 and k3 (locus_entry), and
 * transforms it to the actual k1 and k3 with six cases (select_transform and
 * transform_locus): the wavenumber indices are shifted by the wavenumber
 * offset of k1, the directions are rotated and, in the mirrored cases,
 * reflected, with the matching swap of the interpolation weights. The factor
 * Z is scaled with lambda^(13/2), lambda = q_kfac^(ikmin - 1), for geometric
 * grids. The order of operations matches the Fortran original so results are
 * reproducible.
 *
 * Kept as in the WW3 code (please review): in case 3 the weights of k2 are
 * taken from the reference weights of k2 although the wavenumber index of k2
 * is taken from k4, and in case 6 the wavenumber indices of k2 and k4 are
 * swapped but their direction indices are not. If none of the six cases
 * applies, the indices and weights of the transformed locus keep their
 * previous values, as in the Fortran code. This cannot occur, as k1 = k3 is
 * skipped by the caller.
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param ik3 1-based wavenumber index of k3.
 * @param ia3 1-based direction index of k3.
 * @param[out] locus Transformed locus.
 * @return The result of the lookup.
 * @pre The reference loci are valid for the grid.
 * @post If kFound is returned, locus holds the transformed locus.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_GETLOCUS by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
LocusResult get_locus(const KernelContext &context, const int ik1,
                      const int ia1, const int ik3, const int ia3,
                      TransformedLocus &locus) noexcept {
  const InteractionGrid &quad = context.interactions;
  const LocusEntry entry = locus_entry(context, ik1, ia1, ik3, ia3);
  if (entry.a_entry > context.settings.max_direction_entry) {
    return LocusResult::kInvalidAddress;
  }
  const std::size_t index =
      static_cast<std::size_t>(entry.k_entry - 1) +
      (quad.num_k_entries * static_cast<std::size_t>(entry.a_entry - 1));
  const int nloc = quad.num_points[index];
  if (nloc <= 0) {
    return LocusResult::kEmpty;
  }
  locus.num_points = nloc;

  // Scale factors for geometric grids
  int k_offset = 0;
  float lambda = 1.0F;
  if (context.settings.geometric_grid != 0) {
    k_offset = entry.ikmin - 1;
    lambda = std::pow(context.grid.wavenumber_factor,
                      static_cast<float>(entry.ikmin) - 1.0F);
  }
  const float j_lambda = 1.0F / std::sqrt(lambda);
  const float c_lambda = sixth_power(lambda);
  const float zz_lambda = (lambda * c_lambda) / j_lambda;

  LocusTransform transform{};
  if (select_transform(context, ik1, ia1, ik3, ia3, entry, k_offset,
                       transform)) {
    transform_locus(context, index, nloc, transform, locus);
  }
  const std::size_t stride = quad.num_k_entries * quad.num_a_entries;
  const auto n_points = static_cast<std::size_t>(nloc);
  for (std::size_t ui = 0; ui < n_points; ++ui) {
    locus.factor[ui] = zz_lambda * quad.locus_factor[index + (stride * ui)];
  }
  return LocusResult::kFound;
}

// --- T13Result --------------------------------------------------------------
/**
 * @brief Line integral and diagonal terms of one (k1, k3) combination.
 */
struct T13Result {
  float t13 = 0.0F;     ///< Line integral T13.
  float diag_k1 = 0.0F; ///< Contribution to the diagonal term of k1.
  float diag_k3 = 0.0F; ///< Contribution to the diagonal term of k3.
};

// --- density_at -------------------------------------------------------------
/**
 * @brief Bilinear interpolation of N(k,theta) at a locus point.
 * @details Index and direction neighbours as in Q_T13V4: the upper
 * wavenumber index is limited to nk, the directions are periodic, and for
 * sector grids the direction indices are limited to [1, na]. The tail factor
 * T2 = T4 = 1 of Q_T13V4 is omitted.
 * @param context Kernel context.
 * @param k_index 1-based lower wavenumber index.
 * @param a_index 1-based lower direction index (before periodicity).
 * @param weights Weights of (k, a), (k + 1, a), (k, a + 1), (k + 1, a + 1).
 * @return The interpolated action density.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_T13V4 by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
float density_at(const KernelContext &context, const int k_index,
                 const int a_index,
                 const std::array<float, 4> &weights) noexcept {
  const auto nkq = static_cast<int>(context.num_wavenumbers);
  const auto naq = static_cast<int>(context.num_directions);
  const int jk = k_index;
  const int jkp = std::min(jk + 1, nkq);
  // C++ % truncates towards zero, as Fortran MOD
  int ja = ((a_index - 1 + naq) % naq) + 1;
  int jap = ((a_index + naq) % naq) + 1;
  if (context.settings.grid_type < 3) {
    ja = std::min(std::max(ja, 1), naq);
    jap = std::min(std::max(jap, 1), naq);
  }
  const std::span<const float> n = context.density;
  return (((weights[0] * n[context.at(jk, ja)]) +
           (weights[1] * n[context.at(jkp, ja)])) +
          (weights[2] * n[context.at(jk, jap)])) +
         (weights[3] * n[context.at(jkp, jap)]);
}

// --- compute_q_t13v4 --------------------------------------------------------
/**
 * @brief Computes the line integral T13 and the diagonal terms of one
 * (k1, k3) combination (Q_T13V4).
 * @details T13 = sum Z [N1 N3 (N4 - N2) + N2 N4 (N3 - N1)] along the locus,
 * with N2 and N4 interpolated at the locus points. The diagonal terms are
 * the derivatives of the integrand with respect to N1 and N3. Zero if k1 = k3
 * or the locus is empty. The order of operations matches the Fortran
 * original so results are reproducible.
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param ik3 1-based wavenumber index of k3.
 * @param ia3 1-based direction index of k3.
 * @param locus Work space for the transformed locus.
 * @param[out] result Line integral and diagonal terms.
 * @return The result of the locus lookup (kInvalidAddress on error).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_T13V4 by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
LocusResult compute_q_t13v4(const KernelContext &context, const int ik1,
                            const int ia1, const int ik3, const int ia3,
                            TransformedLocus &locus,
                            T13Result &result) noexcept {
  result = T13Result{};
  if (ik1 == ik3 && ia1 == ia3) {
    return LocusResult::kEmpty;
  }
  const LocusResult found = get_locus(context, ik1, ia1, ik3, ia3, locus);
  if (found != LocusResult::kFound) {
    return found;
  }
  const float qn1 = context.density[context.at(ik1, ia1)];
  const float qn3 = context.density[context.at(ik3, ia3)];
  const float qn13p = qn1 * qn3;
  const float qn13d = qn3 - qn1;

  float t13 = 0.0F;
  float diag_k1 = 0.0F;
  float diag_k3 = 0.0F;
  const auto n_points = static_cast<std::size_t>(locus.num_points);
  for (std::size_t i = 0; i < n_points; ++i) {
    const float qn2 = density_at(context, locus.k2_index[i], locus.a2_index[i],
                                 {
                                     locus.k2_weight1[i],
                                     locus.k2_weight2[i],
                                     locus.k2_weight3[i],
                                     locus.k2_weight4[i],
                                 });
    const float qn4 = density_at(context, locus.k4_index[i], locus.a4_index[i],
                                 {
                                     locus.k4_weight1[i],
                                     locus.k4_weight2[i],
                                     locus.k4_weight3[i],
                                     locus.k4_weight4[i],
                                 });
    const float nprod = (qn13p * (qn4 - qn2)) + ((qn2 * qn4) * qn13d);
    const float rterm = locus.factor[i];
    t13 = t13 + (rterm * nprod);
    const float qd1 = (qn3 * (qn4 - qn2)) - (qn2 * qn4);
    const float qd3 = (qn1 * (qn4 - qn2)) + (qn2 * qn4);
    diag_k1 = diag_k1 + (qd1 * rterm);
    diag_k3 = diag_k3 + (qd3 * rterm);
  }
  result.t13 = t13;
  result.diag_k1 = diag_k1;
  result.diag_k3 = diag_k3;
  return LocusResult::kFound;
}

// --- is_computed ------------------------------------------------------------
/**
 * @brief Filter test of a (k1, k3) combination (Q_XNL4V4).
 * @details The combination is computed if filtering is off, or if the
 * directional difference does not exceed QF_DMAX, the wavenumber ratio does
 * not exceed QF_KRAT, and not both densities are below
 * QF_FRAC max(N) (k_1 / k_3)^7.5.
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param ik3 1-based wavenumber index of k3.
 * @param ia3 1-based direction index of k3.
 * @param density_max Maximum of N(k,theta).
 * @return True if the combination is computed.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_XNL4V4 by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
bool is_computed(const KernelContext &context, const int ik1, const int ia1,
                 const int ik3, const int ia3,
                 const float density_max) noexcept {
  const KernelSettings &settings = context.settings;
  const float qn1 = context.density[context.at(ik1, ia1)];
  const float qn3 = context.density[context.at(ik3, ia3)];
  const float a_dif =
      kHalfCircleDeg -
      std::abs(kHalfCircleDeg -
               std::abs(context.angle_at(ia1) - context.angle_at(ia3)));
  const float k_rat = std::max(context.k_at(ik1) / context.k_at(ik3),
                               context.k_at(ik3) / context.k_at(ik1));
  const float qn_min =
      (settings.filter_fraction * density_max) *
      context.grid.filter_power[static_cast<std::size_t>(ik3 - 1)];
  const bool filter_direction = a_dif > settings.filter_max_direction;
  const bool filter_ratio = k_rat > settings.filter_max_k_ratio;
  const bool filter_density = qn1 < qn_min && qn3 < qn_min;
  return (!filter_direction && !filter_ratio && !filter_density) ||
         settings.use_filter == 0;
}

// --- add_contribution -------------------------------------------------------
/**
 * @brief Adds the contribution of one (k1, k3) combination to the transfer
 * and the diagonal term (Q_XNL4V4).
 * @details Scales T13 and the diagonal terms with Q_SCALE (and 2 with the
 * symmetry), adds T13 k3 dtheta dk3 at k1 and, with the symmetry, subtracts
 * T13 k1 dtheta dk1 at k3. With the symmetry, the transfer at k3 is updated
 * at direction JA3, which for symmetric sector grids with IA3 < IAREF is
 * NAQ - IA1 + 1, while the diagonal term is updated at IA3, as in the
 * Fortran code. The order of the four updates matches the Fortran original.
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param ik3 1-based wavenumber index of k3.
 * @param ia3 1-based direction index of k3.
 * @param t13_result Line integral and diagonal terms.
 * @param[in,out] transfer Transfer of N(k,theta).
 * @param[in,out] diagonal Diagonal term.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_XNL4V4 by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
void add_contribution(const KernelContext &context, const int ik1,
                      const int ia1, const int ik3, const int ia3,
                      const T13Result &t13_result,
                      const std::span<float> transfer,
                      const std::span<float> diagonal) noexcept {
  const KernelSettings &settings = context.settings;
  const bool symmetry = settings.use_symmetry == 1;
  float t13 = t13_result.t13 * settings.grid_scale;
  float diag_k1 = t13_result.diag_k1 * settings.grid_scale;
  float diag_k3 = t13_result.diag_k3 * settings.grid_scale;
  if (symmetry) {
    t13 = 2.0F * t13;
    diag_k1 = 2.0F * diag_k1;
    diag_k3 = 2.0F * diag_k3;
  }
  int ja3 = ia3;
  if (settings.grid_type == 1 && ia3 < settings.reference_direction) {
    ja3 = static_cast<int>(context.num_directions) - ia1 + 1;
  }
  const float delta = context.grid.direction_step;
  const std::size_t i13 = context.at(ik1, ia1);
  transfer[i13] = transfer[i13] +
                  (((t13 * context.k_at(ik3)) * delta) * context.dk_at(ik3));
  if (symmetry) {
    const std::size_t i33 = context.at(ik3, ja3);
    transfer[i33] = transfer[i33] -
                    (((t13 * context.k_at(ik1)) * delta) * context.dk_at(ik1));
  }
  diagonal[i13] = diagonal[i13] + (((diag_k1 * context.k_at(ik3)) * delta) *
                                   context.dk_at(ik3));
  if (symmetry) {
    const std::size_t i33 = context.at(ik3, ia3);
    diagonal[i33] = diagonal[i33] - (((diag_k3 * context.k_at(ik1)) * delta) *
                                     context.dk_at(ik1));
  }
}

// --- add_k1_interactions ----------------------------------------------------
/**
 * @brief Adds all interactions of one k1 to the transfer and the diagonal
 * term (Q_XNL4V4).
 * @details Loops over k3 (all directions, and wavenumbers from k1 upward if
 * the symmetry is used, otherwise all), in the order of the Fortran loops.
 * @param context Kernel context.
 * @param ik1 1-based wavenumber index of k1.
 * @param ia1 1-based direction index of k1.
 * @param density_max Maximum of N(k,theta).
 * @param locus Work space for the transformed locus.
 * @param[in,out] transfer Transfer of N(k,theta).
 * @param[in,out] diagonal Diagonal term.
 * @return False on a locus address error.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_XNL4V4 by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
bool add_k1_interactions(const KernelContext &context, const int ik1,
                         const int ia1, const float density_max,
                         TransformedLocus &locus,
                         const std::span<float> transfer,
                         const std::span<float> diagonal) noexcept {
  const auto nkq = static_cast<int>(context.num_wavenumbers);
  const auto naq = static_cast<int>(context.num_directions);
  const int jk3 = context.settings.use_symmetry == 0 ? 1 : ik1;
  T13Result t13_result;
  for (int ik3 = jk3; ik3 <= nkq; ++ik3) {
    for (int ia3 = 1; ia3 <= naq; ++ia3) {
      if (!is_computed(context, ik1, ia1, ik3, ia3, density_max)) {
        continue;
      }
      if (compute_q_t13v4(context, ik1, ia1, ik3, ia3, locus, t13_result) ==
          LocusResult::kInvalidAddress) {
        return false;
      }
      add_contribution(context, ik1, ia1, ik3, ia3, t13_result, transfer,
                       diagonal);
    }
  }
  return true;
}

} // namespace

// --- compute_q_xnl4v4 -------------------------------------------------------
/**
 * @brief Computes the WRT nonlinear transfer and its diagonal term
 * (Q_XNL4V4).
 * @details The steps are:
 * 1. Convert A(sigma,theta) to N(k,theta) = A c_g / k, with the group
 *    velocity of Z_CMPCG for the computational depth.
 * 2. Loop over k1 (wavenumbers, and directions IAQ1 to IAQ2) and k3, skip
 *    filtered combinations (is_computed), and add the contributions of T13
 *    to the transfer and the diagonal term (add_k1_interactions,
 *    add_contribution).
 * 3. Mirror the result for symmetric sector grids, and convert the transfer
 *    back to A(sigma,theta) with the Jacobian k / c_g.
 *
 * All computations are in single precision, and the order of operations and
 * of the updates of the output arrays matches the Fortran original so
 * results are reproducible. The power lambda^6 of Q_GETLOCUS is computed as
 * gfortran does at -O3 (see sixth_power).
 *
 * On a locus address error, the computation stops and the partly computed
 * outputs are returned, as in Q_XNL4V4. The Fortran caller reports the
 * error. Screen output of Q_XNL4V4 inside the loops (IQ_SCREEN >= 1) is not
 * reproduced.
 *
 * Variables of Q_XNL4V4 that are computed but not used are not computed
 * here: the first estimate of QN_MIN and the filter indicators IFIL_TOT.
 * @param[in] action_density Action density spectrum A(sigma,theta) (size
 * nk * na, aspec), see file description.
 * @param[in] radian_frequency Radian frequencies (size nk, sigma).
 * @param[in] num_directions Number of directions na (NAQ).
 * @param[in] depth Water depth used in the computation (Q_DEPTH).
 * @param[in] gravity Acceleration of gravity (Q_GRAV).
 * @param[in] grid Wavenumber and direction grid.
 * @param[in] interactions Precomputed loci.
 * @param[in] settings Settings of the computation.
 * @param[out] transfer Nonlinear transfer of A(sigma,theta) (size nk * na,
 * xnl).
 * @param[out] diagonal Diagonal term (size nk * na, diag).
 * @return KernelStatus::kSuccess, or the error found.
 * @pre grid.wavenumber, grid.wavenumber_step, grid.filter_power and
 * radian_frequency have size nk; grid.direction_deg has size na; the spectra
 * have size nk * na; the QUAD arrays have the sizes in interactions and hold
 * loci computed for this grid.
 * @post On success, transfer and diagonal are fully populated.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from Xnl Q_XNL4V4 by G. Ph. van Vledder.
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
KernelStatus compute_q_xnl4v4(const std::span<const float> action_density,
                              const std::span<const float> radian_frequency,
                              const std::size_t num_directions,
                              const float depth, const float gravity,
                              const WavenumberGrid &grid,
                              const InteractionGrid &interactions,
                              const KernelSettings &settings,
                              const std::span<float> transfer,
                              const std::span<float> diagonal) noexcept {
  const std::size_t nk = grid.wavenumber.size();
  const std::size_t na = num_directions;
  std::ranges::fill(diagonal, 0.0F);

  // === Section: Action density N(k,theta) =================================
  std::vector<float> group_velocity(nk);
  std::vector<float> density(nk * na);
  for (std::size_t ik = 0; ik < nk; ++ik) {
    group_velocity[ik] = compute_z_cmpcg(radian_frequency[ik], depth, gravity);
    for (std::size_t ia = 0; ia < na; ++ia) {
      const std::size_t i = ik + (nk * ia);
      density[i] =
          (action_density[i] / grid.wavenumber[ik]) * group_velocity[ik];
    }
  }
  std::ranges::fill(transfer, 0.0F);
  const float density_max = fortran_maxval(density);

  // === Section: Loop over k1 and k3 =======================================
  const KernelContext context{
      .num_wavenumbers = nk,
      .num_directions = na,
      .grid = grid,
      .interactions = interactions,
      .settings = settings,
      .density = density,
  };
  TransformedLocus locus(interactions.num_locus_points);
  const auto nkq = static_cast<int>(nk);
  for (int ik1 = 1; ik1 <= nkq; ++ik1) {
    for (int ia1 = settings.first_direction; ia1 <= settings.last_direction;
         ++ia1) {
      if (!add_k1_interactions(context, ik1, ia1, density_max, locus, transfer,
                               diagonal)) {
        return KernelStatus::kInvalidLocusAddress;
      }
    }
  }

  // === Section: Mirror for symmetric sector grids =========================
  if (settings.grid_type == 1) {
    const auto naq = static_cast<int>(na);
    for (int ikq = 1; ikq <= nkq; ++ikq) {
      for (int iaq = settings.reference_direction; iaq <= naq; ++iaq) {
        const int jaq = naq + 1 - iaq;
        transfer[context.at(ikq, jaq)] = transfer[context.at(ikq, iaq)];
      }
    }
  }

  // === Section: Transfer of A(sigma,theta) ================================
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const float jacobian = grid.wavenumber[ik] / group_velocity[ik];
    for (std::size_t ia = 0; ia < na; ++ia) {
      const std::size_t i = ik + (nk * ia);
      transfer[i] = transfer[i] * jacobian;
    }
  }
  return KernelStatus::kSuccess;
}

} // namespace ww4_nl2

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- q_xnl4v4_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to ww4_nl2::compute_q_xnl4v4.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_nl2::compute_q_xnl4v4. Bound in mod_xnl4v5.f90 as
 * BIND(C, name="q_xnl4v4_cpp"). Arguments have the names of the Xnl library
 * variables, see NL2_ww3.h.
 * @param[in] aspec Action density spectrum (size nkq * naq).
 * @param[in] sigma Radian frequencies (size nkq).
 * @param[in] nkq Number of wavenumbers.
 * @param[in] naq Number of directions.
 * @param[in] q_depth Water depth used in the computation.
 * @param[in] q_grav Acceleration of gravity.
 * @param[in] q_scale Scale factor from the grid search.
 * @param[in] q_k Wavenumbers (size nkq).
 * @param[in] q_dk Wavenumber steps (size nkq).
 * @param[in] q_kpow Filter powers of the wavenumbers (size nkq).
 * @param[in] q_ad Directions in degrees (size naq).
 * @param[in] q_delta Directional step in radians.
 * @param[in] q_deltad Directional step in degrees.
 * @param[in] q_kfac Geometric wavenumber factor.
 * @param[in] qf_frac Filter fraction of the maximum action density.
 * @param[in] qf_dmax Filter maximum directional difference (degrees).
 * @param[in] qf_krat Filter maximum wavenumber ratio.
 * @param[in] iq_filt Filtering switch.
 * @param[in] iq_sym Symmetry switch.
 * @param[in] iq_grid Grid type.
 * @param[in] iq_geom Geometric scaling switch.
 * @param[in] iaq1 First direction of the k1 loop.
 * @param[in] iaq2 Last direction of the k1 loop.
 * @param[in] iaref Reference direction.
 * @param[in] iamax Largest valid direction entry.
 * @param[in] quad_nloc Points per locus (size mkq * maq).
 * @param[in] quad_ik2 k index of k2 (size mkq * maq * klocus).
 * @param[in] quad_ia2 Direction index of k2 (size mkq * maq * klocus).
 * @param[in] quad_ik4 k index of k4 (size mkq * maq * klocus).
 * @param[in] quad_ia4 Direction index of k4 (size mkq * maq * klocus).
 * @param[in] quad_w1k2 Weight 1 of k2 (size mkq * maq * klocus).
 * @param[in] quad_w2k2 Weight 2 of k2 (size mkq * maq * klocus).
 * @param[in] quad_w3k2 Weight 3 of k2 (size mkq * maq * klocus).
 * @param[in] quad_w4k2 Weight 4 of k2 (size mkq * maq * klocus).
 * @param[in] quad_w1k4 Weight 1 of k4 (size mkq * maq * klocus).
 * @param[in] quad_w2k4 Weight 2 of k4 (size mkq * maq * klocus).
 * @param[in] quad_w3k4 Weight 3 of k4 (size mkq * maq * klocus).
 * @param[in] quad_w4k4 Weight 4 of k4 (size mkq * maq * klocus).
 * @param[in] quad_zz Factor along locus (size mkq * maq * klocus).
 * @param[in] mkq First dimension of the QUAD arrays.
 * @param[in] maq Second dimension of the QUAD arrays.
 * @param[in] klocus Third dimension of the QUAD arrays.
 * @param[out] xnl Nonlinear transfer (size nkq * naq).
 * @param[out] diag Diagonal term (size nkq * naq).
 * @return 0 on success, 1 if a locus direction entry exceeds iamax.
 * @pre All arrays are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-25
 */
int q_xnl4v4_cpp(const float *const aspec, const float *const sigma,
                 const int nkq, const int naq, const float q_depth,
                 const float q_grav, const float q_scale,
                 const float *const q_k, const float *const q_dk,
                 const float *const q_kpow, const float *const q_ad,
                 const float q_delta, const float q_deltad, const float q_kfac,
                 const float qf_frac, const float qf_dmax, const float qf_krat,
                 const int iq_filt, const int iq_sym, const int iq_grid,
                 const int iq_geom, const int iaq1, const int iaq2,
                 const int iaref, const int iamax, const int *const quad_nloc,
                 const int *const quad_ik2, const int *const quad_ia2,
                 const int *const quad_ik4, const int *const quad_ia4,
                 const float *const quad_w1k2, const float *const quad_w2k2,
                 const float *const quad_w3k2, const float *const quad_w4k2,
                 const float *const quad_w1k4, const float *const quad_w2k4,
                 const float *const quad_w3k4, const float *const quad_w4k4,
                 const float *const quad_zz, const int mkq, const int maq,
                 const int klocus, float *const xnl,
                 float *const diag) noexcept {
  const auto nk = static_cast<std::size_t>(nkq);
  const auto na = static_cast<std::size_t>(naq);
  const auto mk = static_cast<std::size_t>(mkq);
  const auto ma = static_cast<std::size_t>(maq);
  const auto nl = static_cast<std::size_t>(klocus);
  const std::size_t n_quad = mk * ma * nl;

  const ww4_nl2::WavenumberGrid grid{
      .wavenumber = std::span<const float>(q_k, nk),
      .wavenumber_step = std::span<const float>(q_dk, nk),
      .filter_power = std::span<const float>(q_kpow, nk),
      .direction_deg = std::span<const float>(q_ad, na),
      .direction_step = q_delta,
      .direction_step_deg = q_deltad,
      .wavenumber_factor = q_kfac,
  };
  const ww4_nl2::InteractionGrid interactions{
      .num_k_entries = mk,
      .num_a_entries = ma,
      .num_locus_points = nl,
      .num_points = std::span<const int>(quad_nloc, mk * ma),
      .k2_index = std::span<const int>(quad_ik2, n_quad),
      .a2_index = std::span<const int>(quad_ia2, n_quad),
      .k4_index = std::span<const int>(quad_ik4, n_quad),
      .a4_index = std::span<const int>(quad_ia4, n_quad),
      .k2_weight1 = std::span<const float>(quad_w1k2, n_quad),
      .k2_weight2 = std::span<const float>(quad_w2k2, n_quad),
      .k2_weight3 = std::span<const float>(quad_w3k2, n_quad),
      .k2_weight4 = std::span<const float>(quad_w4k2, n_quad),
      .k4_weight1 = std::span<const float>(quad_w1k4, n_quad),
      .k4_weight2 = std::span<const float>(quad_w2k4, n_quad),
      .k4_weight3 = std::span<const float>(quad_w3k4, n_quad),
      .k4_weight4 = std::span<const float>(quad_w4k4, n_quad),
      .locus_factor = std::span<const float>(quad_zz, n_quad),
  };
  const ww4_nl2::KernelSettings settings{
      .use_filter = iq_filt,
      .use_symmetry = iq_sym,
      .grid_type = iq_grid,
      .geometric_grid = iq_geom,
      .first_direction = iaq1,
      .last_direction = iaq2,
      .reference_direction = iaref,
      .max_direction_entry = iamax,
      .filter_fraction = qf_frac,
      .filter_max_direction = qf_dmax,
      .filter_max_k_ratio = qf_krat,
      .grid_scale = q_scale,
  };
  const ww4_nl2::KernelStatus status = ww4_nl2::compute_q_xnl4v4(
      std::span<const float>(aspec, nk * na), std::span<const float>(sigma, nk),
      na, q_depth, q_grav, grid, interactions, settings,
      std::span<float>(xnl, nk * na), std::span<float>(diag, nk * na));
  return static_cast<int>(status);
}

} // extern "C"
