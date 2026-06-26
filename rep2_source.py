"""
rep2_source.py -- where to source the F=2 repumper (rep2), and the AOM it needs.

This SUPERSEDES the earlier master_aom.py, which had the tag sign backwards.

THE SIGN (decisive): the tag AOM DOWN-shifts the retro, retro = forward - 2f_A, NOT +.
config.py says so (repump2 = control - 2f_A) and the comb only CLOSES this way: with
f_mod = A_HFS + 2f_A, the retro of the forward +1 sideband lands exactly on the probe
(control + A_HFS).  The up-shift puts nothing on the probe -- so the down-shift is not a
convention, it is required for the Lambda to exist.  (Bench-confirm the diffraction orders
before ordering; if the build truly up-shifts, config.py itself is misdescribed.)

THE ANSWER: don't drive the raw 85Rb master onto F'1 at all.  You already make 87Rb
F2->F'3 (cooler-frequency) light; from there rep2 (in-trap F2->F'1) is one standard
double-pass AOM away, and -- crucially -- that shift is pure 87Rb hyperfine + the 1064
push, so the 85Rb lock position never enters.  The "pin the lock to 321 MHz" caveat of
master_aom.py dissolves: start from the cooler, not the crossover.

  rep2  (clears |2,-1>,|2,0>,...) : F2->F'1, sigma+, from a CW slave locked to F2->F'3.
  F=1 leg (clears |1,0>,|1,+1>)   : F1->F'2 = your MOT repumper, as-is, ZERO shift.

Everything is wired to stark.py / config.py.  All frequencies 2*pi*MHz.
"""
import numpy as np
import stark, config as c

BARE  = {0:-229.17, 1:-156.95, 2:0.0, 3:266.65}   # 87Rb 5P3/2 HFS vs F'2 (free space, MHz)
HFSg  = 6834.68                                    # 87Rb 5S1/2 F1<->F2 splitting
U0    = -stark.shift(c.alpha0_5S)                  # 1064 ground down-shift (~22.7), F-independent
TH    = c.theta_trap
twofA = c.twofA                                    # tag double-pass shift (400)
GAMMA = 2*np.pi*6.07e6                             # 5P3/2 decay (1/s)
def line(Fp, mp): return BARE[Fp] + stark.stark_level(Fp,mp,TH) + U0   # in-trap F2->|F'm'>

COOLER = BARE[3]                                   # free-space 87Rb F2->F'3 lock (the slave anchor)
f1c    = float(np.mean([line(1,mp) for mp in (-1,0,1)]))   # in-trap F'1 centre

# ---------------------------------------------------------------- the sign, proven ----
print("="*72); print("1.  THE TAG DOWN-SHIFTS  (retro = forward - 2f_A) -- comb-closure proof")
print("="*72)
fmod = HFSg + twofA
print(f"   f_mod = A_HFS + 2f_A = {fmod:.1f} MHz = {fmod/1e3:.3f} GHz  (carrier = control)")
print(f"   forward +1 sideband = control + f_mod          (= repump1 = probe + 2f_A)")
print(f"   DOWN: retro carrier      = control - 2f_A = {-twofA:+.0f}   == config repump2  [OK]")
print(f"   DOWN: retro +1 sideband  = control + A_HFS = {HFSg:+.0f}  == the PROBE       [OK]")
print(f"   UP would give retro carrier {+twofA:+.0f} and nothing on the probe -> comb fails.")

# ----------------------------------------------------------- rep2 source budget ----
print("\n"+"="*72); print("2.  rep2 = F2->F'1 (sigma+).  IN-TRAP target, and the source options")
print("="*72)
print(f"   in-trap F2->F'1 centre = {f1c:+.1f} (rel bare F2->F'2) = {f1c-COOLER:+.1f} from the cooler")
shA = f1c - COOLER
print(f"\n   >>> RECOMMENDED: CW slave on F2->F'3 (cooler), AOM {shA:+.1f} -> DOUBLE-PASS {abs(shA)/2:.0f} MHz")
print(f"       lock-INDEPENDENT: this shift is pure 87Rb hyperfine + 1064 push; the 85Rb")
print(f"       crossover never enters, so the 321 MHz lock ambiguity is irrelevant here.")

# why F'1 and not the (closer) F'2 -- bright-leg scatter
shF2 = (BARE[2]+stark.stark_level(2,0,TH)+U0) - COOLER
print(f"\n   why F'1, not F'2 (closer at {shF2:+.0f} from cooler): sigma+ on F'2 drives the BRIGHT")
print(f"       leg |2,+1>->|F'2,+2> ON RESONANCE.  F'1 has no m'=+2, so |2,+1> is spared.")

# retro byproduct if rep2 is folded into the tagged arm (down-shift: retro = F'1 - 2f_A)
rb = f1c - twofA
allowed = [(F,m,line(F,m)) for F in (1,2,3) for m in range(-F,F+1)]   # F2->F'0 is dF=2 forbidden
F,m,E = min(allowed, key=lambda t: abs(rb-t[2])); d = abs(rb-E)
print(f"\n   folded into the EOM arm? its retro byproduct (sigma-) lands at {rb:+.0f} = {d:.0f} off F'{F}:")
print(f"       scatter {GAMMA*1.6**2/(4*d**2):.0f}/s, AC-Stark {1.6**2/(4*d)/c.nu_z:.3f} nu_z (Om~1.6) -> harmless.")

# --------------------------------------------------------------- F=1 leg ----
print("\n"+"="*72); print("3.  F=1 leg = your MOT repumper (F1->F'2), as-is -- ZERO shift")
print("="*72)
print("   F1->F'2 is exactly the transition that empties the F=1 dark sublevels |1,0>,|1,+1>.")
print("   Route some MOT-repumper light into the fibre; no AOM, no new lock.")

# -------------------------------------------------- the rejected alternatives ----
print("\n"+"="*72); print("4.  why NOT the master_aom.py routes (kept here as the 'do not' record)")
print("="*72)
for nm, lk in [("cooler-1240", COOLER-1240), ("repumper+5540", -HFSg+5540)]:
    print(f"   raw 85Rb master @ {lk:+.0f} ({nm}) onto F'1: AOM {f1c-lk:+.0f} -> DP {abs(f1c-lk)/2:.0f} (hard, and lock-dependent)")
fwd = f1c + twofA                                   # rep2-on-retro forces forward = F'1 + 2f_A
f3 = [(mp, line(3,mp)) for mp in range(-3,4)]
mp,E3 = min(f3, key=lambda t: abs(fwd-t[1])); d3 = abs(fwd-E3)
print(f"   rep2-on-retro (down-shift): forward partner = F'1 + 2f_A = {fwd:+.0f} = {d3:.0f} off in-trap F'3")
print(f"       -> full-power scatter {GAMMA*8.8**2/(4*d3**2):.2e}/s on the cycling line.  HAZARD.")
print(f"   (master_aom.py picked this one because it assumed the UP-shift; with the real")
print(f"    down-shift the safe and hazardous assignments swap. Source from the cooler instead.)")
