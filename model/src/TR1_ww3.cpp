/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file TR1_ww3.cpp
 * @brief Implementation of the single-point Lumped Triad Approximation (TR1)
 * triad interaction source term routine for WAVEWATCH IV (WW4).
 * @details Pure C++20 implementation of the TR1 triad nonlinear interaction
 * source term (W3STR1) for one spatial point, with a C-compatible interface
 * called from the WAVEWATCH III (WW3) Fortran module w3str1md.F90. Starting
 * point of the WW4 conversion of the TR1 source term. The input and output
 * spectrum (wave action density N(k,theta)) is described in TR1_ww3.h. Each
 * routine has its full documentation (details, contracts and heritage) in
 * front of its definition, and a shorter version with its declaration in
 * TR1_ww3.h.
 *
 * References:
 * - Eldeberky, Y., 1996: Nonlinear transformation of wave spectra in the
 *   nearshore zone. Ph.D. thesis, Delft University of Technology, Delft, The
 *   Netherlands, 203 pp.
 * - Madsen, P. A. and O. R. Sorensen, 1993: Bound waves and triad interactions
 *   in shallow water. Ocean Eng., 20(4), 359-388.
 * - WAVEWATCH III Development Group, WAVEWATCH III user manual, section TR1.
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
 * w3str1md.F90 (W3STR1) by A. J. van der Westhuysen, with contributions by A.
 * Roland, based on the SWAN code by M. Zijlema.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */

// Local
#include "TR1_ww3.h"

// System
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace ww4_tr1 {

namespace {

// --- Constants --------------------------------------------------------------
/// Single precision pi, as the WW3 constant PI (3.141592653589793 rounded to
/// single precision).
constexpr float kPi = std::numbers::pi_v<float>;
static_assert(kPi == 3.141592653589793F); // NOLINT(modernize-use-std-numbers)

/// 2 pi (TPI in WW3).
constexpr float kTwoPi = 2.0F * kPi;

/// pi / 2, factor in the biphase.
constexpr float kHalfPi = 0.5F * kPi;

/// Smallest positive normal float (TINY(1.) in WW3). Nothing is computed if
/// the action density is below this value for all spectral components.
constexpr float kMinSpectrum = std::numeric_limits<float>::min();

/// Action density below which the source and diagonal terms are set to zero
/// (JACEPS in W3STR1).
constexpr float kMinAction = 1.0e-12F;

/// Proportionality constant alpha (PTRIAD(1) in W3STR1).
constexpr float kAlpha = 1.0F;

/// Interactions are computed for frequencies below this multiple of the mean
/// frequency sigma_m01 (PTRIAD(2) in W3STR1).
constexpr float kMaxFrequencyFactor = 10.0F;

/// Constant 0.2 in the biphase beta = (pi / 2) (tanh(0.2 / Ur) - 1)
/// (PTRIAD(4) in W3STR1).
constexpr float kBiphaseConstant = 0.2F;

/// Interactions are only computed for Ursell numbers of at least this value
/// (PTRIAD(5) in W3STR1).
constexpr float kMinUrsell = 0.01F;

/// Ratio H_s / sqrt(E) = 4.
constexpr float kSignificantHeightFactor = 4.0F;

/// Factor 2 / 15 in the denominator of the interaction coefficient J.
constexpr float kTwoFifteenths = 2.0F / 15.0F;

/// Factor 2 / 5 in the denominator of the interaction coefficient J.
constexpr float kTwoFifths = 2.0F / 5.0F;

// --- integer_power ----------------------------------------------------------
/**
 * @brief Integer power of a value for a run-time exponent.
 * @details Uses the square-and-multiply algorithm of libgcc (__powisf2) that
 * gfortran uses for x**n with variable n. The algorithm defines the rounding,
 * so it must not be changed.
 * @param value The value x.
 * @param exponent The integer exponent n, which may be negative.
 * @return The power x**n.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
float integer_power(float value, const int exponent) noexcept {
  unsigned int remaining = exponent < 0
                               ? 0U - static_cast<unsigned int>(exponent)
                               : static_cast<unsigned int>(exponent);
  float result = ((remaining % 2U) != 0U) ? value : 1.0F;
  while ((remaining >>= 1U) != 0U) {
    value = value * value;
    if ((remaining % 2U) != 0U) {
      result = result * value;
    }
  }
  return exponent < 0 ? 1.0F / result : result;
}

// --- positive_part ----------------------------------------------------------
/**
 * @brief Positive part max(0, x) as evaluated by gfortran for MAX(0., x).
 * @details Returns x if x > 0 and +0 otherwise, also for x = -0 and NaN. This
 * is the result of the maxss instruction that gfortran generates at -O3, the
 * WW3 release optimization. At -O0, gfortran returns -0 for x = -0.
 * @param value The value x.
 * @return The positive part of x.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
constexpr float positive_part(const float value) noexcept {
  return value > 0.0F ? value : 0.0F;
}

} // namespace

// --- compute_w3str1 ---------------------------------------------------------
/**
 * @brief Computes the TR1 Lumped Triad Approximation source term (W3STR1).
 * @details Self-sum triad interactions of Eldeberky (1996) per direction, see
 * file description in TR1_ww3.h. The steps are:
 * 1. Total energy and first moment of the spectrum, including the parametric
 *    tail (FTE, FTF), giving H_s, sigma_m01 and the Ursell number Ur.
 * 2. Indices and weights of the components at 2 sigma and sigma / 2 on the
 *    logarithmic frequency grid, from the frequency ratio XIS of two bands in
 *    the middle of the grid.
 * 3. Highest band for the interactions: the last band below 10 sigma_m01, and
 *    at least the band at 2 sigma_1.
 * 4. If Ur >= 0.01: per direction, the energy density spectrum E and the
 *    contributions S+ for the bands up to the highest band, then
 *    S = S+(sigma) - 2 S+(2 sigma) converted to action density for all bands
 *    except the highest. Otherwise S and D are set to zero.
 *
 * All computations are in single precision, and the order of operations
 * matches the Fortran original so results are reproducible. x**n with a
 * run-time exponent n uses the gfortran (libgcc) algorithm, and max(0, x)
 * returns +0 for x = -0 as gfortran at -O3.
 *
 * As in W3STR1, nothing is written to source_term and diagonal_term if the
 * spectrum is zero (all values below TINY(1.)), and if Ur >= 0.01, nothing is
 * written for the highest wavenumber band.
 *
 * Variables of W3STR1 that are computed but not used are not computed here:
 * the band energies ED, the directional sum EF, the interpolated energy EDM
 * and the maximum wave height HMAX.
 * @param[in] action_density Action density spectrum N(k,theta) as 1-D array
 * (size nspec = nth * nk), see file description.
 * @param[in] group_velocity Group velocity c_g of each frequency band (size
 * nk).
 * @param[in] wavenumber Wavenumber k of each frequency band (size nk).
 * @param[in] depth Mean water depth d.
 * @param[in,out] source_term Triad source term S in action density form as
 * 1-D array (size nspec).
 * @param[in,out] diagonal_term Diagonal term D of the source term derivative
 * as 1-D array (size nspec).
 * @param[in] num_directions Number of directions nth.
 * @param[in] radian_frequency Radian frequency sigma of each frequency band
 * (size nk).
 * @param[in] integration_factor Factor DDEN = d_sigma d_theta sigma of each
 * frequency band for the spectral integrals (size nk).
 * @param[in] tail_energy_factor Tail factor for the energy integral (FTE).
 * @param[in] tail_frequency_factor Tail factor for the first moment (FTF).
 * @param[in] gravity Acceleration of gravity g.
 * @pre group_velocity, wavenumber, radian_frequency and integration_factor
 * have the same size nk >= 4. action_density, source_term and diagonal_term
 * have size nth * nk with nth = num_directions. The frequency grid is
 * logarithmic with a ratio such that the band at 2 sigma_1 is within the grid.
 * @post If the spectrum is not zero, source_term and diagonal_term are set for
 * all bands except the highest, and for all bands if Ur < 0.01.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @author Code Heritage: Converted from WW3 W3STR1 by A. J. van der
 * Westhuysen, with contributions by A. Roland, based on the SWAN code by M.
 * Zijlema.
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void compute_w3str1(const std::span<const float> action_density,
                    const std::span<const float> group_velocity,
                    const std::span<const float> wavenumber, const float depth,
                    const std::span<float> source_term,
                    const std::span<float> diagonal_term,
                    const std::size_t num_directions,
                    const std::span<const float> radian_frequency,
                    const std::span<const float> integration_factor,
                    const float tail_energy_factor,
                    const float tail_frequency_factor,
                    const float gravity) noexcept {
  const std::size_t num_wavenumbers = wavenumber.size();
  const auto num_wavenumbers_int = static_cast<int>(num_wavenumbers);

  // === Section: Zero spectrum, nothing to compute =========================
  if (std::ranges::all_of(action_density, [](const float value) {
        return value < kMinSpectrum;
      })) {
    return;
  }

  // === Section: Energy, mean frequency and Ursell number ==================
  std::vector<float> band_energy(num_wavenumbers);
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    float sum = 0.0F;
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      sum += action_density[ith + (num_directions * ik)];
    }
    band_energy[ik] = (sum * integration_factor[ik]) / group_velocity[ik];
  }
  float total_energy = 0.0F;
  float first_moment = 0.0F;
  for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
    total_energy = total_energy + band_energy[ik];
    first_moment = first_moment + (band_energy[ik] * radian_frequency[ik]);
  }
  const float tail_band = band_energy[num_wavenumbers - 1] /
                          integration_factor[num_wavenumbers - 1];
  total_energy = total_energy + (tail_band * tail_energy_factor);
  first_moment = first_moment + (tail_band * tail_frequency_factor);
  const float mean_frequency = first_moment / total_energy;

  const float significant_height =
      kSignificantHeightFactor * std::sqrt(positive_part(total_energy));
  // 2 sqrt(2); sqrt2_v is sqrt(2) rounded to single precision, as SQRT(2.).
  const float two_sqrt_two = 2.0F * std::numbers::sqrt2_v<float>;
  const float ursell =
      (gravity * significant_height) /
      ((two_sqrt_two * (mean_frequency * mean_frequency)) * (depth * depth));

  const float depth_sq = depth * depth;
  const float depth_cu = (depth * depth) * depth;
  const float gravity_depth = gravity * depth;

  // === Section: Indices and weights of the interacting components =========
  const int band2 = num_wavenumbers_int / 2; // I2, 1-based
  const int band1 = band2 - 1;               // I1, 1-based
  const float ratio = radian_frequency[static_cast<std::size_t>(band2 - 1)] /
                      radian_frequency[static_cast<std::size_t>(band1 - 1)];
  const float log_ratio = std::log(ratio);

  // ln2_v is log(2) rounded to single precision, as LOG(2.).
  const auto plus = static_cast<int>(std::numbers::ln2_v<float> / log_ratio);
  const int plus1 = plus + 1;
  const float plus_weight =
      (2.0F - integer_power(ratio, plus)) /
      (integer_power(ratio, plus1) - integer_power(ratio, plus));
  const float plus1_weight = 1.0F - plus_weight;

  const auto minus = static_cast<int>(std::log(0.5F) / log_ratio);
  const int minus1 = minus - 1;
  const float minus_weight =
      (integer_power(ratio, minus) - 0.5F) /
      (integer_power(ratio, minus) - integer_power(ratio, minus1));
  const float minus1_weight = 1.0F - minus_weight;

  // === Section: Highest band for the interactions =========================
  int max_band = 1; // ISMAX, 1-based
  for (int ik = 1; ik <= num_wavenumbers_int; ++ik) {
    if (radian_frequency[static_cast<std::size_t>(ik - 1)] <
        (kMaxFrequencyFactor * mean_frequency)) {
      max_band = ik;
    }
  }
  max_band = std::max(max_band, plus1);

  // === Section: Triad interactions ========================================
  // Written as !(Ur >= 0.01) so that Ur = NaN also gives zero, as in W3STR1.
  if (!(ursell >= kMinUrsell)) {
    std::ranges::fill(diagonal_term, 0.0F);
    std::ranges::fill(source_term, 0.0F);
    return;
  }

  const float biphase = kHalfPi * (std::tanh(kBiphaseConstant / ursell) - 1.0F);
  const float sin_biphase = std::abs(std::sin(biphase));

  // Contributions S+ per direction and band, zero beyond the highest band.
  // The extra plus1 bands hold the zeros read for 2 sigma near the end of the
  // grid (SA is dimensioned (NTH,100) and set to zero in W3STR1).
  const std::size_t num_bands_plus =
      num_wavenumbers + static_cast<std::size_t>(plus1);
  std::vector<float> contribution(num_directions * num_bands_plus, 0.0F);
  std::vector<float> energy(num_wavenumbers);

  for (std::size_t ith = 0; ith < num_directions; ++ith) {
    for (std::size_t ik = 0; ik < num_wavenumbers; ++ik) {
      energy[ik] = ((action_density[ith + (num_directions * ik)] * kTwoPi) *
                    radian_frequency[ik]) /
                   group_velocity[ik];
    }
    for (int ik = 1; ik <= max_band; ++ik) {
      const auto ikz = static_cast<std::size_t>(ik - 1);
      const float energy0 = energy[ikz];
      const float frequency0 = radian_frequency[ikz];
      const float wavenumber0 = wavenumber[ikz];
      const float phase_velocity0 = frequency0 / wavenumber0;
      float energy_m = 0.0F;
      float frequency_m = 0.0F;
      float wavenumber_m = 0.0F;
      float phase_velocity_m = 0.0F;
      if (ik > -minus1) {
        const auto im1 = static_cast<std::size_t>(ik + minus1 - 1);
        const auto im = static_cast<std::size_t>(ik + minus - 1);
        energy_m = (minus_weight * energy[im1]) + (minus1_weight * energy[im]);
        frequency_m = (minus_weight * radian_frequency[im1]) +
                      (minus1_weight * radian_frequency[im]);
        wavenumber_m =
            (minus_weight * wavenumber[im1]) + (minus1_weight * wavenumber[im]);
        phase_velocity_m = frequency_m / wavenumber_m;
      }
      const float aux1 =
          (wavenumber_m * wavenumber_m) *
          (gravity_depth + (2.0F * (phase_velocity_m * phase_velocity_m)));
      const float aux2 =
          (wavenumber0 * depth) *
          ((gravity_depth + (((kTwoFifteenths * gravity) * depth_cu) *
                             (wavenumber0 * wavenumber0))) -
           ((kTwoFifths * (frequency0 * frequency0)) * depth_sq));
      const float interaction = aux1 / aux2;
      const float factor = (((kAlpha * phase_velocity0) * group_velocity[ikz]) *
                            (interaction * interaction)) *
                           sin_biphase;
      contribution[ith + (num_directions * ikz)] = positive_part(
          factor * ((energy_m * energy_m) - ((2.0F * energy_m) * energy0)));
    }
  }

  // === Section: Source term in action density form ========================
  const auto plusz = static_cast<std::size_t>(plus);
  const auto plus1z = static_cast<std::size_t>(plus1);
  for (std::size_t ik = 0; ik + 1 < num_wavenumbers; ++ik) {
    const float jacobian =
        (radian_frequency[ik] * kTwoPi) / group_velocity[ik]; // SIGPICG
    for (std::size_t ith = 0; ith < num_directions; ++ith) {
      const std::size_t is = ith + (num_directions * ik);
      const float total =
          contribution[is] -
          (2.0F * ((plus_weight *
                    contribution[ith + (num_directions * (ik + plus1z))]) +
                   (plus1_weight *
                    contribution[ith + (num_directions * (ik + plusz))])));
      if (action_density[is] > kMinAction) {
        diagonal_term[is] = total / (action_density[is] * jacobian);
        source_term[is] = total / jacobian;
      } else {
        diagonal_term[is] = 0.0F;
        source_term[is] = 0.0F;
      }
    }
  }
}

} // namespace ww4_tr1

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- w3str1_cpp -------------------------------------------------------------
/**
 * @brief Fortran interface to ww4_tr1::compute_w3str1.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * ww4_tr1::compute_w3str1. Bound in w3str1md.F90 as
 * BIND(C, name="w3str1_cpp"). S and D are only partly set (see
 * ww4_tr1::compute_w3str1), so they are INTENT(INOUT) in the Fortran
 * interface.
 * @param[in] A Action density spectrum as 1-D array (size nth * nk).
 * @param[in] CG Group velocities (size nk).
 * @param[in] WN Wavenumbers (size nk).
 * @param[in] depth Mean water depth.
 * @param[in,out] S Triad source term as 1-D array (size nth * nk).
 * @param[in,out] D Diagonal term of the derivative as 1-D array (size
 * nth * nk).
 * @param[in] nth Number of directions.
 * @param[in] nk Number of wavenumbers.
 * @param[in] sig Radian frequencies (size nk).
 * @param[in] dden Integration factor DDEN (size nk).
 * @param[in] fte Tail factor for the energy integral (FTE).
 * @param[in] ftf Tail factor for the first moment (FTF).
 * @param[in] grav Acceleration of gravity.
 * @pre nth > 0, nk >= 4 and all arrays are allocated with the sizes given
 * above.
 * @author Main Author(s): Aldgisl (AI Persona), Hendrik L. Tolman
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-24
 * @date Last update: 2026-11-24
 */
void w3str1_cpp(const float *const A, const float *const CG,
                const float *const WN, const float depth, float *const S,
                float *const D, const int nth, const int nk,
                const float *const sig, const float *const dden,
                const float fte, const float ftf, const float grav) noexcept {
  const auto nkz = static_cast<std::size_t>(nk);
  const auto nthz = static_cast<std::size_t>(nth);
  ww4_tr1::compute_w3str1(
      std::span<const float>(A, nthz * nkz), std::span<const float>(CG, nkz),
      std::span<const float>(WN, nkz), depth, std::span<float>(S, nthz * nkz),
      std::span<float>(D, nthz * nkz), nthz, std::span<const float>(sig, nkz),
      std::span<const float>(dden, nkz), fte, ftf, grav);
}

} // extern "C"
