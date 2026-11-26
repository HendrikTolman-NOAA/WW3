/**
 *       +--------------------------------------------------------+
 *       | WAVEWATCH IV, open source, code management by NOAA/NWS |
 *       +--------------------------------------------------------+
 *
 * @file xnl4v5_physics.h
 * @brief Declarations of the physics functions of the Xnl library for the WRT
 * (exact) nonlinear four-wave interactions.
 * @details Declares the C++20 translations of the physics functions of the
 * Xnl library by G. Ph. van Vledder (mod_xnl4v5.f90 and serv_xnl4v5.f90),
 * used by WAVEWATCH III (WW3) for the NL2 source term, and their C-compatible
 * interfaces called from the corresponding Fortran functions in
 * mod_xnl4v5.f90. The Xnl library is kept as a separate library (namespace
 * xnl4v5, CMake target xnl4v5_lib). Implementations are in
 * xnl4v5_physics.cpp. The xnl4v5 routines use descriptive variable names. The
 * extern "C" wrappers keep the Xnl Fortran names for easy cross-referencing
 * with mod_xnl4v5.f90.
 *
 * This is the second stage of the translation of the Xnl library. The
 * functions are used in the computation of the interaction grid (loci of
 * resonant wavenumbers), which itself remains in Fortran for now:
 * - Dispersion relation (X_DISPER) and wavenumber (Z_WNUMB) and group
 *   velocity (Z_CMPCG) from the radian frequency, also used by the kernel.
 * - Locus functions, whose zeros define the locus of k2 for given k1 and k3:
 *   along the symmetry axis (X_LOCUS1), perpendicular to it (X_LOCUS2), in
 *   the wavenumber plane (X_FLOCUS), and the direction of the locus in polar
 *   coordinates (X_COSK).
 * - Jacobian of the delta function in the frequency (X_JACOBIAN), the
 *   magnitude of the difference of the group velocities of k2 and k4.
 * - Coupling coefficient (X_CPLE): deep water coefficient of Webb (1978)
 *   (XC_WEBB) or finite depth coefficient of Herterich and Hasselmann (1980)
 *   (XC_HH).
 *
 * The locus functions use the parameters of the current (k1, k3)
 * combination, which are module variables in the Fortran code: the
 * difference of the radian frequencies q = sigma_1 - sigma_3, the
 * difference vector p = k1 - k3, and the midpoint of the locus on the
 * symmetry axis. They are passed in a LocusParameters structure.
 *
 * Kept as in the WW3 code (please review):
 * - The deep water locus functions (X_FLOCUS, X_LOCUS1, X_LOCUS2 with
 *   IQ_DISP = 1) use sqrt(g) of module m_constants (g = 9.81), while the
 *   dispersion relation X_DISPER uses the gravity of the Xnl library (9.806
 *   in WW3). XC_HH uses g = 9.81 internally.
 * - X_CPLE returns an undefined value for IQ_CPLE = 3 and 4 in the Fortran
 *   code (no coefficient for these values, while Q_CHKCONFIG accepts 3). The
 *   C++ code returns zero.
 * - X_COSK returns an undefined value for IQ_DISP other than 1 and 2 in the
 *   Fortran code. This cannot occur (Q_CHKCONFIG requires 1 or 2). The C++
 *   code uses zero before the limitation to [-1, 1].
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

#ifndef XNL4V5_PHYSICS_H
#define XNL4V5_PHYSICS_H

namespace xnl4v5 {

// --- LocusParameters --------------------------------------------------------
/**
 * @brief Parameters of the locus of one (k1, k3) combination.
 * @details Module variables of the Xnl library set for each (k1, k3)
 * combination in the grid computation, and the depth and dispersion
 * settings.
 */
struct LocusParameters {
  float frequency_difference; ///< q = sigma_1 - sigma_3 (Q).
  float difference_x;         ///< x-component of p = k1 - k3 (PX).
  float difference_y;         ///< y-component of p = k1 - k3 (PY).
  float difference_magnitude; ///< Magnitude of p (PMAG).
  float midpoint_x;           ///< x-component of the locus midpoint (KMIDX).
  float midpoint_y;           ///< y-component of the locus midpoint (KMIDY).
  float depth;                ///< Water depth (Q_DEPTH).
  float gravity;              ///< Acceleration of gravity (Q_GRAV).
  float sqrt_gravity;         ///< sqrt(g) of module m_constants (SQRTG).
  int dispersion;             ///< 1 deep water, 2 finite depth (IQ_DISP).
};

// --- compute_z_wnumb --------------------------------------------------------
/**
 * @brief Wavenumber from radian frequency and depth (Z_WNUMB).
 * @param radian_frequency Radian frequency sigma.
 * @param depth Water depth d.
 * @param gravity Acceleration of gravity g.
 * @return The wavenumber, or -10 for zero depth or frequency.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_z_wnumb(float radian_frequency, float depth,
                                    float gravity) noexcept;

// --- compute_z_cmpcg --------------------------------------------------------
/**
 * @brief Group velocity from radian frequency and depth (Z_CMPCG).
 * @param radian_frequency Radian frequency sigma.
 * @param depth Water depth d.
 * @param gravity Acceleration of gravity g.
 * @return The group velocity, or -10 for zero depth or frequency.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-25
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_z_cmpcg(float radian_frequency, float depth,
                                    float gravity) noexcept;

// --- compute_x_disper -------------------------------------------------------
/**
 * @brief Radian frequency from wavenumber and depth (X_DISPER).
 * @param wavenumber Wavenumber k.
 * @param depth Water depth d.
 * @param dispersion 1 deep water, 2 finite depth (IQ_DISP).
 * @param gravity Acceleration of gravity g (Q_GRAV).
 * @return The radian frequency, or -1 for an invalid dispersion type.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_disper(float wavenumber, float depth,
                                     int dispersion, float gravity) noexcept;

// --- compute_x_cosk ---------------------------------------------------------
/**
 * @brief Cosine of the polar angle of the locus for a wavenumber (X_COSK).
 * @param wavenumber Wavenumber k along the symmetry axis of the locus.
 * @param locus Parameters of the locus.
 * @return The cosine, limited to [-1, 1].
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_cosk(float wavenumber,
                                   const LocusParameters &locus) noexcept;

// --- compute_x_flocus -------------------------------------------------------
/**
 * @brief Locus function in the wavenumber plane (X_FLOCUS).
 * @param kx x-component of the wavenumber k2.
 * @param ky y-component of the wavenumber k2.
 * @param locus Parameters of the locus.
 * @return q + sigma_2 - sigma_4, zero on the locus.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_flocus(float kx, float ky,
                                     const LocusParameters &locus) noexcept;

// --- compute_x_locus1 -------------------------------------------------------
/**
 * @brief Locus function along the symmetry axis (X_LOCUS1).
 * @param wavenumber Magnitude of the wavenumber k2 on the symmetry axis.
 * @param locus Parameters of the locus.
 * @return q + sigma_2 - sigma_4, zero on the locus.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_locus1(float wavenumber,
                                     const LocusParameters &locus) noexcept;

// --- compute_x_locus2 -------------------------------------------------------
/**
 * @brief Locus function perpendicular to the symmetry axis (X_LOCUS2).
 * @param lambda Distance from the midpoint, in units of |p|.
 * @param locus Parameters of the locus.
 * @return q + sigma_2 - sigma_4, zero on the locus.
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_locus2(float lambda,
                                     const LocusParameters &locus) noexcept;

// --- compute_x_jacobian -----------------------------------------------------
/**
 * @brief Jacobian term |c_g2 - c_g4| of the locus integral (X_JACOBIAN).
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
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_jacobian(float k2x, float k2y, float k4x,
                                       float k4y, float depth,
                                       float gravity) noexcept;

// --- compute_xc_webb --------------------------------------------------------
/**
 * @brief Deep water coupling coefficient of Webb (1978) (XC_WEBB).
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
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_xc_webb(float k1x, float k1y, float k2x, float k2y,
                                    float k3x, float k3y, float k4x, float k4y,
                                    float gravity) noexcept;

// --- compute_xc_hh ----------------------------------------------------------
/**
 * @brief Finite depth coupling coefficient of Herterich and Hasselmann (1980)
 * (XC_HH).
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
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_xc_hh(float w1x, float w1y, float w2x, float w2y,
                                  float w3x, float w3y, float depth) noexcept;

// --- compute_x_cple ---------------------------------------------------------
/**
 * @brief Coupling coefficient for a quadruplet (X_CPLE).
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
 * @return The coupling coefficient (zero for other values of coupling).
 * @author Main Author(s): Gerbrant Ph. van Vledder (Xnl library), Aldgisl
 * (AI Persona) and Hendrik L. Tolman (C++ translation)
 * @author Contributors: Kiro (Agentic AI/SI)
 * @date Initial: 2026-11-26
 * @date Last update: 2026-11-26
 */
[[nodiscard]] float compute_x_cple(float k1x, float k1y, float k2x, float k2y,
                                   float k3x, float k3y, float k4x, float k4y,
                                   int coupling, float depth,
                                   float gravity) noexcept;

} // namespace xnl4v5

// --- C-compatible ABI interfaces for Fortran interoperability ---------------
extern "C" {

/**
 * @brief Fortran interface to xnl4v5::compute_x_disper (X_DISPER).
 * @param k Wavenumber.
 * @param d Water depth.
 * @param iq_disp Dispersion type.
 * @param q_grav Acceleration of gravity.
 * @return The radian frequency.
 */
float x_disper_cpp(float k, float d, int iq_disp, float q_grav) noexcept;

/**
 * @brief Fortran interface to xnl4v5::compute_x_cosk (X_COSK).
 * @param k Wavenumber along the symmetry axis.
 * @param q Difference of radian frequencies.
 * @param pmag Magnitude of k1 - k3.
 * @param iq_disp Dispersion type.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The cosine of the polar angle.
 */
float x_cosk_cpp(float k, float q, float pmag, int iq_disp, float q_depth,
                 float q_grav) noexcept;

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
 */
float x_flocus_cpp(float kxx, float kyy, float q, float px, float py,
                   float sqrtg, int iq_disp, float q_depth,
                   float q_grav) noexcept;

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
 */
float x_locus1_cpp(float k2, float q, float pmag, float sqrtg, int iq_disp,
                   float q_depth, float q_grav) noexcept;

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
 */
float x_locus2_cpp(float lambda, float q, float px, float py, float kmidx,
                   float kmidy, float sqrtg, int iq_disp, float q_depth,
                   float q_grav) noexcept;

/**
 * @brief Fortran interface to xnl4v5::compute_x_jacobian (X_JACOBIAN).
 * @param x2 x-component of k2.
 * @param y2 y-component of k2.
 * @param x4 x-component of k4.
 * @param y4 y-component of k4.
 * @param q_depth Water depth.
 * @param q_grav Acceleration of gravity.
 * @return The Jacobian term.
 */
float x_jacobian_cpp(float x2, float y2, float x4, float y4, float q_depth,
                     float q_grav) noexcept;

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
 */
float x_cple_cpp(float k1x, float k1y, float k2x, float k2y, float k3x,
                 float k3y, float k4x, float k4y, int iq_cple, float depth,
                 float grav) noexcept;

} // extern "C"

#endif // XNL4V5_PHYSICS_H
