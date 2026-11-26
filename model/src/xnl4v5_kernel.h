/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file xnl4v5_kernel.h
 * @brief Declarations of the single-point interaction kernel of the Xnl
 * library for the WRT (exact) nonlinear four-wave interactions.
 * @details Declares the C++20 kernel of the Webb-Resio-Tracy (WRT) method of
 * the Xnl library by G. Ph. van Vledder (mod_xnl4v5.f90), used by WAVEWATCH
 * III (WW3) for the NL2 source term, and its C-compatible interface called
 * from Q_XNL4V4 in mod_xnl4v5.f90. The Xnl library is kept as a separate
 * library (namespace xnl4v5, CMake target xnl4v5_lib), independent of the
 * WW3 and WW4 source term code. Implementations are in xnl4v5_kernel.cpp.
 * The xnl4v5 routine uses descriptive variable names. The extern "C" wrapper
 * keeps the Xnl Fortran names (aspec, q_k, quad_zz, ...) for easy
 * cross-referencing with mod_xnl4v5.f90.
 *
 * This is the first stage of the translation of the Xnl library. The kernel
 * covers the computations of Q_XNL4V4, Q_T13V4 and Q_GETLOCUS (and the
 * service routines Z_CMPCG and Z_WNUMB): the nonlinear transfer and its
 * diagonal term for a given spectrum, using an interaction grid that is
 * computed (or read from a BQF file) by the Fortran code. Grid generation,
 * grid selection, depth scaling, file I/O, logging and error handling remain
 * in Fortran for now.
 *
 * Method: the rate of change of action density N_1 at wavenumber k_1 is
 * written as dN_1/dt = 2 int T(k_1, k_3) dk_3, where the function T13 is a
 * line integral along the locus of wavenumbers k_2 (and k_4 = k_1 + k_2 -
 * k_3) that satisfy the resonance conditions for the given k_1 and k_3:
 * T(k_1, k_3) = sum over the locus of Z [N_1 N_3 (N_4 - N_2) + N_2 N_4 (N_3 -
 * N_1)]. The precomputed factor Z combines the coupling coefficient, the
 * step along the locus, the symmetry factor and the Jacobian. N_2 and N_4 are
 * bilinear interpolations of the spectrum at the locus points. Loci are only
 * stored for a reference k_1 (the first wavenumber and reference direction),
 * and are transformed to other k_1 and k_3 using the geometric scaling of the
 * wavenumber grid (factor lambda^(13/2) on Z for deep water) and rotation and
 * mirroring of the directions. Interactions are filtered on the directional
 * difference and wavenumber ratio of k_1 and k_3 and on the action density,
 * and the symmetry of the integrand can be used to halve the work. The
 * diagonal term for implicit integration is the derivative of the integrand
 * with respect to N_1.
 *
 * Input and output spectrum: the input is the action density spectrum
 * A(sigma,theta) as a function of radian frequency and direction, as passed
 * to xnl_main by W3SNL2, stored in Fortran (column-major) order
 * aspec(ik, ia), so that the 1-D index of wavenumber ik and direction ia
 * (both 0-based) is ik + nk * ia. It is converted to N(k,theta) internally
 * and the transfer is converted back, so that xnl is the rate of change of
 * A(sigma,theta), stored in the same order. Note that this is the transposed
 * order of the WW3 spectra.
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
 * Vledder: mod_xnl4v5.f90 (Q_XNL4V4, Q_T13V4, Q_GETLOCUS) and serv_xnl4v5.f90
 * (Z_CMPCG, Z_WNUMB), as included in WAVEWATCH III (WW3) for the NL2 source
 * term (w3snl2md.F90 by Hendrik L. Tolman and G. Ph. van Vledder).
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */

#ifndef XNL4V5_KERNEL_H
#define XNL4V5_KERNEL_H

// System
#include <cstddef>
#include <cstdint>
#include <span>

namespace xnl4v5 {

// --- WavenumberGrid ---------------------------------------------------------
/**
 * @brief Wavenumber and direction grid of the interaction computation.
 * @details Non-owning views of the grid variables set by Q_INIT in the Xnl
 * library.
 */
struct WavenumberGrid {
  std::span<const float> wavenumber;      ///< Wavenumbers k (Q_K, size nk).
  std::span<const float> wavenumber_step; ///< Wavenumber steps (Q_DK, size nk).
  std::span<const float> filter_power;    ///< (k_1 / k)^7.5 (Q_KPOW, size nk).
  std::span<const float> direction_deg;   ///< Directions in degrees (Q_AD).
  float direction_step;     ///< Directional step in radians (Q_DELTA).
  float direction_step_deg; ///< Directional step in degrees (Q_DELTAD).
  float wavenumber_factor;  ///< Geometric wavenumber factor (Q_KFAC).
};

// --- InteractionGrid --------------------------------------------------------
/**
 * @brief Precomputed loci of the reference interactions.
 * @details Non-owning views of the QUAD_* arrays of the Xnl library, with
 * dimensions (mk, ma, nl) in Fortran order: index (km, am, il) (all 0-based)
 * is km + mk * (am + ma * il). quad_nloc has dimensions (mk, ma).
 */
struct InteractionGrid {
  std::size_t num_k_entries;           ///< First dimension mk (MKQ).
  std::size_t num_a_entries;           ///< Second dimension ma (MAQ).
  std::size_t num_locus_points;        ///< Third dimension nl (KLOCUS).
  std::span<const int> num_points;     ///< Points per locus (QUAD_NLOC).
  std::span<const int> k2_index;       ///< Lower k index of k2 (QUAD_IK2).
  std::span<const int> a2_index;       ///< Lower direction of k2 (QUAD_IA2).
  std::span<const int> k4_index;       ///< Lower k index of k4 (QUAD_IK4).
  std::span<const int> a4_index;       ///< Lower direction of k4 (QUAD_IA4).
  std::span<const float> k2_weight1;   ///< Weight 1 of k2 (QUAD_W1K2).
  std::span<const float> k2_weight2;   ///< Weight 2 of k2 (QUAD_W2K2).
  std::span<const float> k2_weight3;   ///< Weight 3 of k2 (QUAD_W3K2).
  std::span<const float> k2_weight4;   ///< Weight 4 of k2 (QUAD_W4K2).
  std::span<const float> k4_weight1;   ///< Weight 1 of k4 (QUAD_W1K4).
  std::span<const float> k4_weight2;   ///< Weight 2 of k4 (QUAD_W2K4).
  std::span<const float> k4_weight3;   ///< Weight 3 of k4 (QUAD_W3K4).
  std::span<const float> k4_weight4;   ///< Weight 4 of k4 (QUAD_W4K4).
  std::span<const float> locus_factor; ///< Factor Z along locus (QUAD_ZZ).
};

// --- KernelSettings ---------------------------------------------------------
/**
 * @brief Settings of the interaction computation.
 * @details Switches, filter values and loop limits of the Xnl library. All
 * direction indices are 1-based, as in the Fortran code.
 */
struct KernelSettings {
  int use_filter;      ///< Filtering on (1) or off (0) (IQ_FILT).
  int use_symmetry;    ///< Symmetry reduction on (1) or off (0) (IQ_SYM).
  int grid_type;       ///< 1 symmetric sector, 2 sector, 3 circle (IQ_GRID).
  int geometric_grid;  ///< Geometric scaling of loci used (1) (IQ_GEOM).
  int first_direction; ///< First direction of the k1 loop (IAQ1).
  int last_direction;  ///< Last direction of the k1 loop (IAQ2).
  int reference_direction;    ///< Reference direction of k1 (IAREF).
  int max_direction_entry;    ///< Largest valid direction entry (IAMAX).
  float filter_fraction;      ///< Fraction of the maximum density (QF_FRAC).
  float filter_max_direction; ///< Largest |theta1 - theta3| (QF_DMAX).
  float filter_max_k_ratio;   ///< Largest k ratio of k1 and k3 (QF_KRAT).
  float grid_scale;           ///< Scale factor from the grid search (Q_SCALE).
};

// --- KernelStatus -----------------------------------------------------------
/**
 * @brief Result status of the interaction kernel.
 */
enum class KernelStatus : std::uint8_t {
  kSuccess = 0,             ///< Computation completed.
  kInvalidLocusAddress = 1, ///< Locus direction entry above IAMAX (MEMORY).
};

// --- compute_q_xnl4v4 -------------------------------------------------------
/**
 * @brief Computes the WRT nonlinear transfer and its diagonal term
 * (Q_XNL4V4).
 * @details See file description.
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
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */
[[nodiscard]] KernelStatus compute_q_xnl4v4(
    std::span<const float> action_density,
    std::span<const float> radian_frequency, std::size_t num_directions,
    float depth, float gravity, const WavenumberGrid &grid,
    const InteractionGrid &interactions, const KernelSettings &settings,
    std::span<float> transfer, std::span<float> diagonal) noexcept;

} // namespace xnl4v5

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

// --- q_xnl4v4_cpp -----------------------------------------------------------
/**
 * @brief Fortran interface to xnl4v5::compute_q_xnl4v4.
 * @details Wraps the raw Fortran arrays in non-owning spans and forwards to
 * xnl4v5::compute_q_xnl4v4. Bound in mod_xnl4v5.f90 as
 * BIND(C, name="q_xnl4v4_cpp"). Arguments have the names of the Xnl library
 * variables.
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
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */
int q_xnl4v4_cpp(const float *aspec, const float *sigma, int nkq, int naq,
                 float q_depth, float q_grav, float q_scale, const float *q_k,
                 const float *q_dk, const float *q_kpow, const float *q_ad,
                 float q_delta, float q_deltad, float q_kfac, float qf_frac,
                 float qf_dmax, float qf_krat, int iq_filt, int iq_sym,
                 int iq_grid, int iq_geom, int iaq1, int iaq2, int iaref,
                 int iamax, const int *quad_nloc, const int *quad_ik2,
                 const int *quad_ia2, const int *quad_ik4, const int *quad_ia4,
                 const float *quad_w1k2, const float *quad_w2k2,
                 const float *quad_w3k2, const float *quad_w4k2,
                 const float *quad_w1k4, const float *quad_w2k4,
                 const float *quad_w3k4, const float *quad_w4k4,
                 const float *quad_zz, int mkq, int maq, int klocus, float *xnl,
                 float *diag) noexcept;

} // extern "C"

#endif // XNL4V5_KERNEL_H
