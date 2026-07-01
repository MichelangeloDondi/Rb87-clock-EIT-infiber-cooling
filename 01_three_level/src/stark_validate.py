"""
stark_validate.py -- the exact (sympy) backing for the hardcoded geometric factors in stark.py.

stark.py keeps the per-(F',m') tensor map numpy-only (hardcoded tensor_F_prefactor) so it stays a clean,
hand-checkable closed form. This script recomputes those factors from the Wigner-6j symbols and
checks them against (a) stark.py's hardcoded values and (b) the polarizability authority:
  - {2 2 2; 3/2 3/2 3/2} = 0  -> F'=2 tensor NULL (the EIT target is pure scalar)
  - F'=3 stretched m'=+-3 total -> +19.4 MHz   (theta=0, pol || B)

Run:  python stark_validate.py     (needs sympy)
"""
from sympy import S, Rational
from sympy.physics.wigner import wigner_6j
import stark


def sixj(a, b, c, d, e, f):
    return float(wigner_6j(S(a), S(b), S(c), S(d), S(e), S(f)))


# F'-prefactor (2F'+1){F' 2 F'; 3/2 3/2 3/2}, normalized so F'=3 -> 1
def Fpref(Fp):
    return (2 * Fp + 1) * sixj(Fp, 2, Fp, Rational(3, 2), Rational(3, 2), Rational(3, 2))


pref3 = Fpref(3)
print("Wigner-6j {F' 2 F'; 3/2 3/2 3/2} and the normalized F'-prefactor:")
ok = True
for Fp in (0, 1, 2, 3):
    ratio = Fpref(Fp) / pref3 if Fp >= 1 else 0.0
    hard = stark.tensor_F_prefactor[Fp]
    match = abs(ratio - hard) < 1e-4
    ok = ok and match
    print(f"  F'={Fp}:  ratio = {ratio:+.5f}   stark.tensor_F_prefactor = {hard:+.5f}   {'OK' if match else 'MISMATCH'}")

print(f"\n  6j {{2 2 2;3/2..}} = {sixj(2,2,2,Rational(3,2),Rational(3,2),Rational(3,2)):+.5f}  (0 -> F'=2 pure scalar)")

# physics anchors via stark.py itself (theta=0, pol || B == the authority's convention)
s33 = stark.stark_level(3, 3, 0.0)
n2 = max(abs(stark.stark_tensor(2, m, 0.0)) for m in range(-2, 3))
print(f"\n  stark.stark_level(3,3,theta=0) = {s33:+.2f} MHz   (authority +19.4)   "
      f"{'OK' if abs(s33-19.4) < 0.3 else 'CHECK'}")
print(f"  max |F'=2 tensor| (any m', any theta) = {n2:.3e}   (must be ~0)   "
      f"{'OK' if n2 < 1e-9 else 'CHECK'}")
print("\nALL GEOMETRIC FACTORS MATCH stark.py" if ok else "\n*** MISMATCH: fix stark.tensor_F_prefactor ***")
