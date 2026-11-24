/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file BT4_ww3.cpp
 * @brief Implementation of the single-point SHOWEX movable-bed bottom friction
 * (BT4) source term routines for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the BT4 SHOWEX bottom friction source
 * term (W3SBT4) for one spatial point and of its initialization (TABU_ERF),
 * with C-compatible interfaces called from the WAVEWATCH III (WW3) Fortran
 * module w3sbt4md.F90. Starting point of the WW4 conversion of the BT4 source
 * term. The input and output spectrum (wave action density N(k,theta)) and the
 * parts of the WW3 code that are kept for review are described in BT4_ww3.h.
 * Each routine has its full documentation (details, contracts and heritage) in
 * front of its definition, and a shorter version with its declaration in
 * BT4_ww3.h.
 *
 * References:
 * - Ardhuin, F., T. G. Drake and T. H. C. Herbers, 2002: Observations of
 *   wave-generated vortex ripples on the North Carolina continental shelf.
 *   J. Geophys. Res., 107(C10), 3143.
 * - Ardhuin, F., W. C. O'Reilly, T. H. C. Herbers and P. F. Jessen, 2003:
 *   Swell transformation across the continental shelf. Part I: Attenuation
 *   and directional broadening. J. Phys. Oceanogr., 33, 1921-1939.
 * - Grant, W. D. and O. S. Madsen, 1979: Combined wave and current interaction
 *   with a rough bottom. J. Geophys. Res., 84, 1797-1808.
 * - Tolman, H. L., 1994: Wind-waves and moveable-bed bottom-friction. J. Phys.
 *   Oceanogr., 24, 994-1009.
 * - Tolman, H. L., 1995: Sub-grid modeling of moveable-bed bottom-friction in
 *   wind-wave models. Coastal Eng., 26, 57-75.
 * - Wilson, K. C., 1989: Friction of wave-induced sheet flow. Coastal Eng.,
 *   12, 371-379.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section BT4.
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
 * w3sbt4md.F90 (W3SBT4, INSBT4 and TABU_ERF) by F. Ardhuin and J. Lepesqueur.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */

// Local
#include "BT4_ww3.h"

// System
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace ww4_bt4 {

namespace {

// --- Constants --------------------------------------------------------------
/// Single precision pi, as the WW3 constant PI (3.141592653589793 rounded to
/// single precision).
constexpr float kPi = std::numbers::pi_v<float>;
static_assert(kPi == 3.141592653589793F); // NOLINT(modernize-use-std-numbers)

/// 2 pi (TPI in WW3).
constexpr float kTwoPi = 2.0F * kPi;

/// Relative depth offsets of the three sub-grid depths (XSUB in W3SBT4). The
/// last one is the central depth, whose bulk quantities are used afterwards.
constexpr std::array<float, 3> kDepthOffsets{-0.001F, 0.001F, 0.0F};

/// Bands with k d >= this value are considered deep water (no friction).
constexpr float kDeepWaterLimit = 6.0F;

/// Lower limit of the grain size D_50 (m).
constexpr float kMinGrainSize = 1.0e-5F;

/// Upper limit of the grain size D_50 (m).
constexpr float kMaxGrainSize = 1.0F;

/// Lower limit of the orbital velocity and displacement variances.
constexpr float kMinVariance = 1.0e-7F;

/// Ripple wavelength factor 1.7 = 2 sqrt(2) 0.6 (Ardhuin et al., 2002).
constexpr float kRippleLengthFactor = 1.7F;

/// Small value added to the squared orbital direction vector.
constexpr float kMinDirectionNorm = 1.0e-12F;

/// Lower limit of a_b / k in the friction factor table lookup.
constexpr float kMinRelativeOrbit = 0.3F;

/// Factor in the sheet flow roughness (Wilson, 1989).
constexpr float kSheetFlowFactor = 0.0655F;

/// Exponent in the sheet flow roughness (Wilson, 1989).
constexpr float kSheetFlowExponent = 1.4F;

/// Relative variation of the Shields number below which no sub-grid treatment
/// is used.
constexpr float kMinRelativeShieldsVariation = 0.0001F;

/// Absolute variation of the Shields number below which no sub-grid treatment
/// is used.
constexpr float kMinShieldsVariation = 1.0e-8F;

/// Small value added to the ripple probability in the mean Shields number.
constexpr float kMinProbability = 0.0001F;

/// Ripple probability above which the ripples are considered active.
constexpr float kActiveRippleProbability = 0.5F;

// --- fortran_int ------------------------------------------------------------
/**
 * @brief Conversion of a real value to an integer as Fortran INT on x86.
 * @details Truncates towards zero. Values outside the integer range and NaN
 * give the smallest integer, as the cvttss2si instruction used by gfortran.
 * In C++ such a conversion is undefined, so the range is tested first.
 * @param value The real value.
 * @return The truncated integer value.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
constexpr int fortran_int(const float value) noexcept {
  constexpr float kIntLimit = 2147483648.0F; // 2^31
  if (value > -kIntLimit && value < kIntLimit) {
    return static_cast<int>(value);
  }
  return std::numeric_limits<int>::min();
}

// --- table_value ------------------------------------------------------------
/**
 * @brief Linear interpolation in a table at a fractional index.
 * @details As in W3SBT4: the integer index is limited to the last interval,
 * and the interpolation weight is limited to 1.
 * @param table Table values (size n + 1).
 * @param position Fractional index x >= 0 into the table.
 * @param lower_index Integer index i = min(n - 1, int(x)), limited by the
 * caller where needed.
 * @return table(i) (1 - w) + table(i + 1) w with w = min(1, x - i).
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
float table_value(const std::span<const float> table, const float position,
                  const int lower_index) noexcept {
  const float upper_weight =
      std::min(1.0F, position - static_cast<float>(lower_index));
  const float lower_weight = 1.0F - upper_weight;
  const auto i = static_cast<std::size_t>(lower_index);
  return (table[i] * lower_weight) + (table[i + 1] * upper_weight);
}

// --- friction_factor --------------------------------------------------------
/**
 * @brief Friction factor from the table FWTABLE for a relative orbital
 * displacement a_b / k.
 * @details Table lookup as in W3SBT4 (sections 2.b and 2.c).
 * @param friction_table Friction factor table FWTABLE(0:SIZEFWTABLE).
 * @param relative_orbit Ratio a_b / k of orbital displacement and roughness.
 * @param friction_log_min Lower limit ABMIN of log10(a_b / k) in the table.
 * @param friction_log_step Step DELAB of log10(a_b / k) in the table.
 * @return The friction factor.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
float friction_factor(const std::span<const float> friction_table,
                      const float relative_orbit, const float friction_log_min,
                      const float friction_log_step) noexcept {
  const int num_intervals = static_cast<int>(friction_table.size()) - 1;
  const float position =
      std::max((std::log10(std::max(relative_orbit, kMinRelativeOrbit)) -
                friction_log_min) /
                   friction_log_step,
               1.0F);
  const int lower_index = std::min(num_intervals - 1, fortran_int(position));
  return table_value(friction_table, position, lower_index);
}

} // namespace

// --- compute_tabu_erf -------------------------------------------------------
/**
 * @brief Tabulates the function (1 + erf(x)) / 2 used in the sub-grid
 * treatment (TABU_ERF).
 * @details Values for x = -x_max + i dx, i = 0 ... n, with dx = 2 x_max / n.
 * For x < 0, WW3 evaluates the value as (1 - |erf(x)|) / 2, which is the same
 * as (1 + erf(x)) / 2. The operations match the Fortran original so results
 * are reproducible.
 *
 * Possible bug in the WW3 code, reproduced here: TABU_ERF multiplies by
 * 2**(1/2), which is the integer 1 (1/2 is integer division), and uses
 * erf(x). The table is used as the probability that the Shields number
 * exceeds the ripple threshold, with x in standard deviations, for which the
 * cumulative normal distribution (1 + erf(x / sqrt(2))) / 2 may have been
 * intended. To be reviewed before use in WW4.
 * @param[out] erf_table Tabulated values (size n + 1 for n intervals,
 * ERFTABLE(0:SIZEERFTABLE)).
 * @param[in] max_stdev Range x_max of the table (XERFMAX).
 * @return The spacing dx of the table (DELXERF).
 * @pre erf_table has size n + 1 with n >= 1.
 * @post erf_table is fully populated.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 TABU_ERF by J. Lepesqueur.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
float compute_tabu_erf(const std::span<float> erf_table,
                       const float max_stdev) noexcept {
  const std::size_t num_intervals = erf_table.size() - 1;
  const float step = (2.0F * max_stdev) / static_cast<float>(num_intervals);
  for (std::size_t i = 0; i <= num_intervals; ++i) {
    const float x = -max_stdev + (static_cast<float>(i) * step);
    // As in WW3: erf(x), not erf(x / sqrt(2)) (possible bug, see @details).
    erf_table[i] = x < 0.0F ? (1.0F - std::abs(std::erf(x))) / 2.0F
                            : (1.0F + std::erf(x)) / 2.0F;
  }
  return step;
}

// --- compute_w3sbt4 ---------------------------------------------------------
/**
 * @brief Computes the BT4 SHOWEX bottom friction source term (W3SBT4).
 * @details Movable-bed bottom friction, see file description in BT4_ww3.h.
 * The steps are:
 * 1. Limit D_50 to [1e-5, 1] and set the background roughness
 *    max(A_5, A_6 D_50).
 * 2. For the three sub-grid depths d (1 - 0.001), d (1 + 0.001) and d: bulk
 *    orbital velocity and displacement and their direction from the bands
 *    with k d < 6, the potential ripple wavelength, the skin friction factor
 *    and the normalized Shields number psi / psi_c.
 * 3. Roughness k_N: without sub-grid treatment if the variation of the
 *    Shields number with depth, d(psi / psi_c) / d(d / d) sigma_d, is small,
 *    and with sub-grid treatment (probability of active ripples from the erf
 *    table) otherwise. The bedform parameters are set.
 * 4. Friction factor f_e for k_N, diagonal term
 *    D = -f_e u_b sigma^2 / (2 g sinh^2(k d)), source term S = D N, and the
 *    momentum flux to the bottom boundary layer.
 *
 * All computations are in single precision, and the order of operations
 * matches the Fortran original so results are reproducible. Fortran MAX(a, b)
 * and MIN(a, b) are written as std::max(a, b) and std::min(a, b) with the
 * same argument order. The sub-grid weights WSUB of W3SBT4 are declared but
 * not used, and are not included here.
 *
 * Possible bug in the WW3 code, reproduced here: the momentum flux TAUBBL is
 * summed with ECOS(IS) and ESIN(IS) using the spectral index IS instead of
 * the direction index ITH (cos_table[is] and sin_table[is] here). As only the
 * first nth elements of these WW3 arrays are set and the others are zero,
 * only the first frequency band contributes to TAUBBL. The intended
 * expression is probably the sum over all bands with ECOS(ITH) and ESIN(ITH).
 * To be reviewed before use in WW4.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nspec = nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] wavenumber Wavenumber k of each frequency band (size nk).
 * @param[in] depth Mean water depth d.
 * @param[in] grain_size Median grain size D_50 (D50).
 * @param[in] critical_shields Critical Shields number psi_c (PSIC).
 * @param[out] bottom_stress Momentum flux to the bottom boundary layer
 * (TAUBBL, size 2).
 * @param[out] bedform Bedform roughness and ripple wavelength components
 * (BEDFORM, size 3).
 * @param[out] source_term Bottom friction source term S_bt as 1-D array (size
 * nspec).
 * @param[out] diagonal_term Diagonal term of the source term derivative as 1-D
 * array (size nspec).
 * @param[in] num_directions Number of directions nth.
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] integration_factor Factor DDEN = d_sigma d_theta sigma of each
 * frequency band (size nk).
 * @param[in] cos_table WW3 array ECOS, cosines of the directions in the first
 * nth elements (size nspec), see file description.
 * @param[in] sin_table WW3 array ESIN, sines of the directions in the first nth
 * elements (size nspec).
 * @param[in] parameters Empirical coefficients SBTCX(1:7): A_1, A_2, A_3, A_4,
 * sigma_d, A_5, A_6 (size 7).
 * @param[in] friction_table Friction factor table FWTABLE(0:SIZEFWTABLE).
 * @param[in] friction_log_min Lower limit of log10(a_b / k_N) in the friction
 * factor table (ABMIN).
 * @param[in] friction_log_step Step of log10(a_b / k_N) in the friction factor
 * table (DELAB).
 * @param[in] erf_table Table of (1 + erf(x)) / 2 (ERFTABLE(0:SIZEERFTABLE)).
 * @param[in] erf_max Range x_max of the erf table (XERFMAX).
 * @param[in] erf_step Spacing of the erf table (DELXERF).
 * @param[in] gravity Acceleration of gravity g.
 * @param[in] sediment_gravity Specific gravity s of the sediment (SED_SG).
 * @pre group_velocity, wavenumber, radian_frequency and integration_factor
 * have the same size nk. action_density, source_term, diagonal_term,
 * cos_table and sin_table have size nth * nk with nth = num_directions.
 * bottom_stress has size 2, bedform size 3 and parameters size 7. The tables
 * have at least two values.
 * @post All output arrays are fully populated.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3SBT4 by F. Ardhuin and J.
 * Lepesqueur.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void compute_w3sbt4(const std::span<const float> action_density,
                    const std::span<const float> group_velocity,
                    const std::span<const float> wavenumber, const float depth,
                    const float grain_size, const float critical_shields,
                    const std::span<float> bottom_stress,
                    const std::span<float> bedform,
                    const std::span<float> source_term,
                    const std::span<float> diagonal_term,
                    const std::size_t num_directions,
                    const std::span<const float> radian_frequency,
                    const std::span<const float> integration_factor,
                    const std::span<const float> cos_table,
                    const std::span<const float> sin_table,
                    const std::span<const float> parameters,
                    const std::span<const float> friction_table,
                    const float friction_log_min, const float friction_log_step,
                    const std::span<const float> erf_table, const float erf_max,
                    const float erf_step, const float gravity,
                    const float sediment_gravity) noexcept {
  const std::size_t num_wavenumbers = wavenumber.size();
  const float ripple_factor = parameters[0];    // A_1, SBTCX(1)
  const float ripple_exponent = parameters[1];  // A_2, SBTCX(2)
  const float ripple_threshold = parameters[2]; // A_3, SBTCX(3)
  const float relict_factor = parameters[3];    // A_4, SBTCX(4)
  const float depth_stdev = parameters[4];      // sigma_d, SBTCX(5)
  const float min_roughness = parameters[5];    // A_5, SBTCX(6)
  const float grain_roughness = parameters[6];  // A_6, SBTCX(7)

  // === Section: Grain size and background roughness =======================
  float d50 = std::max(grain_size, kMinGrainSize);
  d50 = std::min(d50, kMaxGrainSize);
  const float background = std::max(min_roughness, grain_roughness * d50);

  // Submerged weight factor (s - 1) g, as in the Fortran expressions.
  const float submerged_gravity = (sediment_gravity - 1.0F) * gravity;

  // === Section: Sub-grid loop =============================================
  std::vector<float> decay(num_wavenumbers); // CBETA
  std::array<float, 3> shields_values{};
  const std::span<float> shields(shields_values);
  const std::span<const float> depth_offsets(kDepthOffsets);
  float velocity_sq = 0.0F; // UORB2, 2 var(u_b)
  float velocity = 0.0F;    // UORB
  float orbit = 0.0F;       // AORB
  float ripple_x = 0.0F;    // LX
  float ripple_y = 0.0F;    // LY

  for (std::size_t isub = 0; isub < depth_offsets.size(); ++isub) {
    const float sub_depth = depth * (1.0F + depth_offsets[isub]);
    float velocity_var = 0.0F;
    float orbit_var = 0.0F;
    float orbit_x = 0.0F;
    float orbit_y = 0.0F;

    for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
      if ((wavenumber[ik] * sub_depth) < kDeepWaterLimit) {
        float energy = 0.0F;
        float energy_x = 0.0F;
        float energy_y = 0.0F;
        for (std::size_t ith = 0; ith < num_directions; ++ith) {
          const float value = action_density[ith + (num_directions * ik)];
          energy = energy + value;
          energy_x = energy_x + (value * cos_table[ith]);
          energy_y = energy_y + (value * sin_table[ith]);
        }
        const float sigma = radian_frequency[ik];
        const float sinh_kd = std::sinh(wavenumber[ik] * sub_depth);
        decay[ik] = (0.5F * (sigma * sigma)) / (gravity * (sinh_kd * sinh_kd));
        const float factor =
            (((integration_factor[ik] / group_velocity[ik]) * 2.0F) *
             decay[ik]) *
            gravity;
        const float band_var = energy * factor;
        velocity_var = velocity_var + band_var;
        orbit_var = orbit_var + (band_var / (sigma * sigma));
        orbit_x = orbit_x + (energy_x * factor);
        orbit_y = orbit_y + (energy_y * factor);
      } else {
        decay[ik] = 0.0F;
      }
    }

    // RMS orbital velocity and displacement
    velocity_sq = 2.0F * velocity_var;
    velocity = std::sqrt(std::max(kMinVariance, velocity_sq));
    orbit = std::sqrt(std::max(kMinVariance, 2.0F * orbit_var));

    // Potential ripple wavelength (Ardhuin et al., 2002)
    const float norm = std::sqrt(((orbit_x * orbit_x) + (orbit_y * orbit_y)) +
                                 kMinDirectionNorm);
    ripple_x = ((orbit * kRippleLengthFactor) * orbit_x) / norm;
    ripple_y = ((orbit * kRippleLengthFactor) * orbit_y) / norm;

    // Skin friction factor and normalized Shields number
    const float skin_friction = friction_factor(
        friction_table, orbit / d50, friction_log_min, friction_log_step);
    const float shields_number =
        (skin_friction * velocity_sq) / ((2.0F * submerged_gravity) * d50);
    shields[isub] = shields_number / critical_shields;
  }

  // === Section: Roughness =================================================
  const float shields_slope =
      ((shields[1] - shields[0]) / (depth_offsets[1] - depth_offsets[0])) *
      depth_stdev;
  // Sheet flow roughness factor (u_b^2 / ((s - 1) g a_b))^1.4 (Wilson, 1989)
  const auto sheet_flow = [&] {
    return std::pow(velocity_sq / (submerged_gravity * orbit),
                    kSheetFlowExponent);
  };

  float roughness = 0.0F; // KSUBN
  if ((std::abs(shields_slope) < (kMinRelativeShieldsVariation * shields[2])) ||
      (std::abs(shields_slope) < kMinShieldsVariation)) {
    // No sub-grid treatment
    if (shields[2] > ripple_threshold) {
      const float ripple_roughness =
          (orbit * ripple_factor) * std::pow(shields[2], ripple_exponent);
      const float sheet_roughness = (orbit * kSheetFlowFactor) * sheet_flow();
      roughness = ripple_roughness + sheet_roughness;
      bedform[1] = ripple_x;
      bedform[2] = ripple_y;
    } else {
      roughness = std::max(background, orbit * relict_factor);
      bedform[1] = -ripple_x;
      bedform[2] = -ripple_y;
    }
  } else {
    // Sub-grid treatment (Tolman, 1995)
    const float shields_offset =
        (ripple_threshold - shields[2]) / shields_slope;
    const int num_erf_intervals = static_cast<int>(erf_table.size()) - 1;
    const float erf_position =
        std::max((shields_offset + erf_max) / erf_step, 0.0F);
    const int erf_lower =
        std::max(std::min(num_erf_intervals - 1, fortran_int(erf_position)), 0);
    const float ripple_probability = std::max(
        std::min(table_value(erf_table, erf_position, erf_lower), 1.0F), 0.0F);
    const float relict_probability = 1.0F - ripple_probability;
    // Mean normalized Shields number with ripples (Tolman, 1995)
    const float shields_ripples = std::max(
        shields[2] + (((std::exp(-(0.5F * (shields_offset * shields_offset))) /
                        std::sqrt(kTwoPi)) *
                       shields_slope) /
                      (ripple_probability + kMinProbability)),
        ripple_threshold);
    roughness =
        (relict_probability * std::max(background, orbit * relict_factor)) +
        ((ripple_probability * orbit) *
         ((ripple_factor * std::pow(shields_ripples, ripple_exponent)) +
          (kSheetFlowFactor * sheet_flow())));
    if (ripple_probability > kActiveRippleProbability) {
      bedform[1] = ripple_x;
      bedform[2] = ripple_y;
    } else {
      bedform[1] = -ripple_x;
      bedform[2] = -ripple_y;
    }
  }
  bedform[0] = roughness;

  // === Section: Source term and momentum flux =============================
  const float friction = friction_factor(friction_table, orbit / roughness,
                                         friction_log_min, friction_log_step);
  bottom_stress[0] = 0.0F;
  bottom_stress[1] = 0.0F;
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    // Jacobian to band energy times coefficient to momentum
    const float momentum_factor =
        ((integration_factor[ik] / group_velocity[ik]) * gravity) /
        (radian_frequency[ik] / wavenumber[ik]);
    const float diagonal_value = -((friction * velocity) * decay[ik]);
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      const std::size_t is = ith + (num_directions * ik);
      diagonal_term[is] = diagonal_value;
      const float momentum =
          (momentum_factor * diagonal_value) * action_density[is];
      // As in WW3: index is instead of ith, so only the first band
      // contributes (possible bug, see @details).
      bottom_stress[0] = bottom_stress[0] - (momentum * cos_table[is]);
      bottom_stress[1] = bottom_stress[1] - (momentum * sin_table[is]);
      source_term[is] = diagonal_value * action_density[is];
    }
  }
}

} // namespace ww4_bt4

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- tabu_erf_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to ww4_bt4::compute_tabu_erf.
 * @details Bound in w3sbt4md.F90 as BIND(C, name="tabu_erf_cpp").
 * @param[out] erftable Tabulated values ERFTABLE(0:sizeerftable) (size
 * sizeerftable + 1).
 * @param[in] sizeerftable Number of intervals of the table (SIZEERFTABLE).
 * @param[in] xerfmax Range of the table (XERFMAX).
 * @return The spacing of the table (DELXERF).
 * @pre sizeerftable >= 1 and erftable has size sizeerftable + 1.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
float tabu_erf_cpp(float *const erftable, const int sizeerftable,
                   const float xerfmax) noexcept {
  return ww4_bt4::compute_tabu_erf(
      std::span<float>(erftable, static_cast<std::size_t>(sizeerftable) + 1),
      xerfmax);
}

// --- w3sbt4_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_bt4::compute_w3sbt4.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_bt4::compute_w3sbt4. Bound in w3sbt4md.F90 as
 * BIND(C, name="w3sbt4_cpp").
 * @param[in] A Action density spectrum as 1-D array (size nth * nk).
 * @param[in] CG Group velocities (size nk).
 * @param[in] WN Wavenumbers (size nk).
 * @param[in] depth Water depth.
 * @param[in] d50 Median grain size.
 * @param[in] psic Critical Shields number.
 * @param[out] TAUBBL Momentum flux to the bottom boundary layer (size 2).
 * @param[out] BEDFORM Bedform parameters (size 3).
 * @param[out] S Bottom friction source term (size nth * nk).
 * @param[out] D Diagonal term of the derivative (size nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Integration factor DDEN (size nk).
 * @param[in] ecos Array ECOS (size nth * nk).
 * @param[in] esin Array ESIN (size nth * nk).
 * @param[in] sbtcx Coefficients SBTCX(1:7) (size 7).
 * @param[in] fwtable Friction factor table (size sizefwtable + 1).
 * @param[in] sizefwtable Number of intervals of the friction factor table.
 * @param[in] abmin Lower limit of the friction factor table (ABMIN).
 * @param[in] delab Step of the friction factor table (DELAB).
 * @param[in] erftable Erf table (size sizeerftable + 1).
 * @param[in] sizeerftable Number of intervals of the erf table.
 * @param[in] xerfmax Range of the erf table (XERFMAX).
 * @param[in] delxerf Spacing of the erf table (DELXERF).
 * @param[in] grav Acceleration of gravity.
 * @param[in] sed_sg Specific gravity of the sediment (SED_SG).
 * @pre nth > 0, nk > 0, sizefwtable >= 1, sizeerftable >= 1, and all arrays
 * are allocated with the sizes given above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void w3sbt4_cpp(const float *const A, const float *const CG,
                const float *const WN, const float depth, const float d50,
                const float psic, float *const TAUBBL, float *const BEDFORM,
                float *const S, float *const D, const int nth, const int nk,
                const float *const sig, const float *const dden,
                const float *const ecos, const float *const esin,
                const float *const sbtcx, const float *const fwtable,
                const int sizefwtable, const float abmin, const float delab,
                const float *const erftable, const int sizeerftable,
                const float xerfmax, const float delxerf, const float grav,
                const float sed_sg) noexcept {
  constexpr std::size_t kNumParameters = 7;
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  const std::size_t nspec = nthz * nkz;
  ww4_bt4::compute_w3sbt4(
      std::span<const float>(A, nspec), std::span<const float>(CG, nkz),
      std::span<const float>(WN, nkz), depth, d50, psic,
      std::span<float>(TAUBBL, 2), std::span<float>(BEDFORM, 3),
      std::span<float>(S, nspec), std::span<float>(D, nspec), nthz,
      std::span<const float>(sig, nkz), std::span<const float>(dden, nkz),
      std::span<const float>(ecos, nspec), std::span<const float>(esin, nspec),
      std::span<const float>(sbtcx, kNumParameters),
      std::span<const float>(fwtable,
                             static_cast<std::size_t>(sizefwtable) + 1),
      abmin, delab,
      std::span<const float>(erftable,
                             static_cast<std::size_t>(sizeerftable) + 1),
      xerfmax, delxerf, grav, sed_sg);
}

} // extern "C"
