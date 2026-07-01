"""
stark.py -- the optical chain -> the trap. CLOSED FORM (no master equation, no solver):
this is the arithmetic of the chapter-01 README trap/Stark section, runnable so you can
plug in your own power/waist.

Light shift of a level of polarizability alpha in intensity I:   U = - alpha * I / (2 eps0 c)  (the optical-dipole potential, Grimm et al. 2000).
alpha > 0 is pulled DOWN (trapped); alpha < 0 is pushed UP (anti-trapped).
Two 1 W beams counter-propagate; at a lattice antinode the field doubles, so I = 4 x the
single-beam peak 2P/(pi w0^2).
"""
import numpy as np
import config as c

eps0 = 8.8541878e-12; clight = 2.99792458e8; h = 6.62607015e-34; kB = 1.380649e-23
a_au = 1.64877727e-41          # atomic unit of polarizability (C^2 m^2 / J)

I1     = 2 * c.P_beam / (np.pi * c.w0**2)   # single-beam peak intensity
I_anti = 4 * I1                             # standing-wave antinode (field doubles)

def shift(alpha_au):
    """Signed light shift in 2pi*MHz at the antinode (negative = trapped, positive = up)."""
    return -alpha_au * a_au * I_anti / (2 * eps0 * clight) / h / 1e6


# --- per-(F',m') light shift of the 5P_3/2 manifold (scalar common + tensor split) ----------
# The scalar shift (~+38) is common to every (F',m'). The TENSOR part splits by m' and is
# F'-dependent. In the |F',m'> basis (Arora-Sandars / Wigner-6j) the tensor of the stretched
# |F'=3,m'=+-3>=|m_J=3/2> state equals the bare J' tensor, i.e. shift(alpha2) ~ -18.6 MHz.
# Geometric factor (normalized so F_geom(3,3)=1):
#   F_geom(F',m') = [ (2F'+1){F' 2 F'; 3/2 3/2 3/2} / (7 {3 2 3; 3/2 3/2 3/2}) ]
#                   * (3 m'^2 - F'(F'+1)) / (F'(2F'-1))
# The two Wigner-6j (recomputed exactly with sympy in 01_three_level/src/stark_validate.py, which checks these):
#   {2 2 2; 3/2 3/2 3/2} = 0          -> F'=2 tensor NULL: the EIT target is pure scalar (any geometry)
#   {1 2 1; 3/2 3/2 3/2} = -0.16330 ; {3 2 3; 3/2 3/2 3/2} = +0.13093
#   -> F'-prefactor ratios (relative to F'=3): F'=0: 0, F'=1: -0.5345, F'=2: 0, F'=3: +1.
tensor_F_prefactor = {0: 0.0, 1: -0.534522, 2: 0.0, 3: 1.0}    # (2F'+1){F' 2 F';3/2..}/(7{3 2 3;3/2..})


def tensor_m_factor(Fp, mp):
    return (3 * mp**2 - Fp * (Fp + 1)) / (Fp * (2 * Fp - 1)) if Fp >= 1 else 0.0


def stark_tensor(Fp, mp, theta_deg=90.0):
    """Differential TENSOR 1064 light shift of |5P3/2 F',m'> relative to |F'2,0> (2pi*MHz).
       theta = angle of the trap LINEAR polarization to the quantization axis (B || fibre axis).
       The axial lattice is transverse-polarized about the axial B, so theta=90 deg (default);
       theta=0 (pol || B) reproduces the polarizability-authority +19.4 (F'=3 stretched) etc."""
    ang = (3 * np.cos(np.deg2rad(theta_deg))**2 - 1) / 2.0     # +1 at 0 deg, -1/2 at 90 deg
    return shift(c.alpha2_5P32) * tensor_F_prefactor[Fp] * tensor_m_factor(Fp, mp) * ang


def stark_level(Fp, mp, theta_deg=90.0):
    """Total 1064 shift of |5P3/2 F',m'> = common scalar + tensor(F',m') (2pi*MHz)."""
    return shift(c.alpha0_5P32) + stark_tensor(Fp, mp, theta_deg)


def stark_vector(alpha1_au, F, m, ellipticity):
    """Rank-1 (VECTOR) 1064 light shift (2pi*MHz): a fictitious B ~ ellipticity * m_F, ZERO for a linear
       trap (ellipticity=0). alpha1_au is the vector polarizability; the linear-in-m form
       U_vec = shift(alpha1) * ellipticity * (m/F) is illustrative (the exact prefactor needs a reference)."""
    if F == 0 or ellipticity == 0.0:
        return 0.0
    return shift(alpha1_au) * ellipticity * (m / F)


if __name__ == "__main__":
    U0  = -shift(c.alpha0_5S)                       # trap depth (the ground shift is downward)
    sca =  shift(c.alpha0_5P32)                     # 5P_3/2 scalar (alpha<0 -> up)
    lo  =  shift(c.alpha0_5P32 + c.alpha2_5P32)     # F'=3 stretched |m_J|=3/2 : alpha0+alpha2
    hi  =  shift(c.alpha0_5P32 - c.alpha2_5P32)     # most anti-trapped         : alpha0-alpha2

    print("1064 nm lattice from 1 W + 1 W counter-propagating, w0 = %.0f um" % (c.w0*1e6))
    print(f"  antinode intensity            I = {I_anti:.2e} W/m^2")
    print(f"  GROUND 5S (a0={c.alpha0_5S:+.0f})       -> trap depth U0 = {U0:.1f} MHz = {U0*1e6*h/kB*1e6:.0f} uK")
    print(f"  EXCITED 5P_3/2 scalar (a0={c.alpha0_5P32:+.0f}) -> {sca:+.1f} MHz  (anti-trapped)")
    print( "  tensor (a2=%+.0f) splits 5P_3/2:" % c.alpha2_5P32)
    print(f"     |F'=2, all m'>  : tensor NULL (6j{{2 2 2; 3/2 3/2 3/2}}=0) -> pure scalar {sca:+.1f} MHz  <-- EIT target")
    print(f"     |F'=3, m'=+-3>  : a0+a2 -> {lo:+.1f} MHz  (stretched, least anti-trapped)")
    print(f"     |F'=3, |m_J|=1/2>: a0-a2 -> {hi:+.1f} MHz  (most anti-trapped)")
    print(f"  => the whole 5P_3/2 manifold is expelled ({lo:+.0f} to {hi:+.0f} MHz); the EIT")
    print(f"     target |F'2,0> sits at the pure-scalar {sca:+.1f} MHz, geometry-independent.")

    print("\nper-(F',m') 1064 shift  (scalar + tensor),  tensor REL to |F'2,0> in [ ]:")
    print("  theta=0 (pol || B)  validates vs the polarizability authority; theta=90 is the real trap.")
    for th in (0.0, 90.0):
        print(f"  theta = {th:.0f} deg:")
        for Fp in (1, 2, 3):
            row = "    F'=%d: " % Fp + "  ".join(
                f"m'={mp:+d}:{stark_level(Fp, mp, th):+5.1f}[{stark_tensor(Fp, mp, th):+5.1f}]"
                for mp in range(-Fp, Fp + 1))
            print(row)
    print("  CHECK theta=0: F'=3 m'=+-3 total -> +19.4 (doc); F'=2 tensor == 0 (6j null).")
    print("  The REPUMP target is |F'=1,m'=-1>: tensor REL |F'2,0> = "
          f"{stark_tensor(1, -1, 90.0):+.1f} MHz at the real theta=90 trap.")
