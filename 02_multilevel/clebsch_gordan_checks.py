"""
clebsch_gordan_checks.py -- prove the Clebsch-Gordan / line-strength conventions in cooling_multilevel.py
are correct, by reconstructing KNOWN 87Rb D2 facts from the raw CG + Wigner-6j.

clebsch(F,m,q,F',m') = <F m; 1 q | F' m'>  (ground F,m  +  photon 1,q  ->  excited F',m'); m'=m+q.
line_strength(F,F')        = (2F'+1) {1/2 3/2 1; F' F 3/2}^2  = the relative line strength of the F->F' line.

Run:  python clebsch_gordan_checks.py     (needs sympy)
"""
import numpy as np
from sympy import S
from sympy.physics.wigner import wigner_6j
import cooling_multilevel as m

clebsch, line_strength, decay_branching = m.clebsch, m.line_strength, m.decay_branching
ok_all = True


def check(name, got, exp, tol=1e-9):
    global ok_all
    good = abs(got - exp) < tol
    ok_all = ok_all and good
    print(f"   {name:42s} = {got:+.5f}   expect {exp:+.5f}   {'OK' if good else '*** MISMATCH ***'}")


print("1. clebsch spot-checks  <F m; 1 q | F' m'>:")
check("stretched <2,2;1,1|3,3>",          clebsch(2, 2, 1, 3, 3),   1.0)
check("pi-forbidden <2,0;1,0|2,0>",       clebsch(2, 0, 0, 2, 0),   0.0)
check("dF=2 forbidden <1,1;1,-1|3,0>",    clebsch(1, 1, -1, 3, 0),  0.0)
check("<1,0;1,1|2,1> = 1/sqrt2",          clebsch(1, 0, 1, 2, 1),   np.sqrt(0.5))
check("<2,2;1,-1|2,1> = 1/sqrt3",         clebsch(2, 2, -1, 2, 1),  np.sqrt(1.0/3))

print("\n2. Decay m-sum rule  sum_mg clebsch(Fg,mg,mp-mg,Fp,mp)^2  = 1 (allowed) or 0 (dF=2 forbidden):")
for (Fp, mp, Fg, exp) in [(2, 0, 1, 1.0), (2, 1, 2, 1.0), (1, -1, 1, 1.0),
                          (3, 0, 2, 1.0), (3, 0, 1, 0.0), (0, 0, 2, 0.0)]:
    s = sum(clebsch(Fg, mg, mp - mg, Fp, mp)**2 for mg in (mp - 1, mp, mp + 1) if abs(mg) <= Fg)
    check(f"Fp={Fp} mp={mp:+d} -> Fg={Fg}", s, exp)
print("   (=1 makes decay_branching[Fp][Fg]*clebsch^2/tot reproduce decay_branching exactly as the F-level branching)")

print("\n3. F'->F decay branching from raw CG/6j  vs the hardcoded decay_branching table:")
w6 = lambda Fg, Fp: float(wigner_6j(S(1)/2, S(Fg), S(3)/2, S(Fp), S(3)/2, 1))
for Fp in (0, 1, 2, 3):
    raw = {Fg: (2 * Fg + 1) * w6(Fg, Fp)**2 for Fg in (1, 2)}
    tot = sum(raw.values()) or 1.0
    for Fg in (1, 2):
        check(f"F'={Fp} -> F={Fg}", raw[Fg] / tot, decay_branching[Fp][Fg])

print("\n4. Line strengths:")
s2 = [line_strength(2, f) for f in (1, 2, 3)]
print(f"   F=2 -> F'=1:2:3  =  {s2[0]/s2[0]:.0f} : {s2[1]/s2[0]:.0f} : {s2[2]/s2[0]:.0f}   (iconic D2 ratio 1 : 5 : 14)")
check("   F=2 ratio is 1:5:14", s2[1] / s2[0] + 100 * s2[2] / s2[0], 5 + 1400)
s1 = [line_strength(1, f) for f in (0, 1, 2)]
print(f"   F=1 -> F'=0:1:2  =  {s1[0]/s1[0]*2:.0f} : {s1[1]/s1[0]*2:.0f} : {s1[2]/s1[0]*2:.0f}   (D2 ratio 2 : 5 : 5)")
# cross-validate F=1 absorption strengths against the (known) decay branching:
Sym = lambda F, Fp: (2 * F + 1) * line_strength(F, Fp)           # symmetric |<F||d||F'>|^2 (relative)
check("S(1,2) == S(2,2)  [F'=2 decays 1:1]",    Sym(1, 2), Sym(2, 2))
check("S(1,1) == 5*S(2,1) [F'=1 decays 5:1]",   Sym(1, 1), 5 * Sym(2, 1))
check("sum-rule sum_F' line_strength(1,F') = 1/2",       sum(line_strength(1, f) for f in (0, 1, 2)), 0.5)
check("sum-rule sum_F' line_strength(2,F') = 1/2",       sum(line_strength(2, f) for f in (1, 2, 3)), 0.5)

print("\nALL CG / LINE-STRENGTH CHECKS PASS" if ok_all else "\n*** SOME CHECK FAILED ***")
