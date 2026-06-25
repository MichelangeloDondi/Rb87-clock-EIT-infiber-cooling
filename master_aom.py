"""
master_aom.py -- required AOM for the rep2-on-retro master assignment.

In the single-EOM tagged-retro chain the master laser is retro-reflected through
the shared tag AOM (double-passed: +2f_A shift on the round trip).  For rep2-on-retro
we want the retro to land on the in-trap F'=1 manifold (rep2: F=2->F'1, sigma+):

    retro = forward + 2f_A = F'1_centre  -->  forward = F'1_centre - 2f_A

This script derives the required master AOM shift from the 85Rb lock frequency
and checks the forward byproduct (must be benign).

WARNING: two stated descriptions of the 85Rb lock disagree by ~320 MHz (see below).
         Measure the lock frequency before ordering.

Conventions: all frequencies in 2*pi*MHz; bare 87Rb F=2->F'=2 = 0; blue = positive.
"""
import numpy as np
import stark
import config as c

# ── 87Rb reference transitions used to specify the 85Rb F3-F4 crossover lock ──
BARE  = {0: -229.17, 1: -156.95, 2: 0.0, 3: 266.65}  # 5P3/2 HFS vs F'=2 (free-space, MHz)
HFSg  = 6834.68            # 87Rb 5S1/2 F=1<->F=2 splitting (MHz)
U0    = -stark.shift(c.alpha0_5S)   # ground 5S down-shift (positive; same for F=1 and F=2)
TH    = c.theta_trap

# "cooler"   = F=2->F'=3 cycling/imaging transition  (+266.65 MHz vs F2->F'2)
# "repumper" = F=1->F'=2 repumper transition         (-6834.68 MHz vs same reference)
COOLER = BARE[3]
REPUMP = -HFSg + BARE[2]

# ── 85Rb F3-F4 crossover: two stated offsets that disagree ─────────────────────
lock_A = COOLER - 1240.0    # "1240 MHz below the cooler"    -> -973 MHz
lock_B = REPUMP + 5540.0    # "5540 MHz above the repumper"  -> -1295 MHz

print("85Rb F3-F4 crossover lock (vs bare 87Rb F2->F'2 = 0):")
print(f"  'cooler - 1240'    ->  {lock_A:+.1f} MHz")
print(f"  'repumper + 5540'  ->  {lock_B:+.1f} MHz")
print(f"  *** these disagree by {abs(lock_A - lock_B):.0f} MHz -- measure before ordering ***\n")

# ── in-trap F'=1 transition frequencies (F=2->F'1 sigma+, relevant for rep2) ──
def line(Fp, mp):
    """In-trap laser freq for F=2->|F'=Fp,m'=mp> (vs bare F2->F'2 = 0, in 2pi*MHz)."""
    return BARE[Fp] + stark.stark_level(Fp, mp, TH) + U0

fp1 = {mp: line(1, mp) for mp in (-1, 0, 1)}
f1c = float(np.mean(list(fp1.values())))           # centre of F'=1 sigma+ accessible lines

print("In-trap F'=1 lines for F=2->F'1 sigma+ (laser freq vs bare F2->F'2 = 0):")
for mp, f in fp1.items():
    print(f"  m'={mp:+d}  {f:+.1f} MHz")
print(f"  centre  {f1c:+.1f} MHz  <- retro target\n")

# ── rep2-on-retro: forward target and required master AOM ─────────────────────
twofA = c.twofA                    # 400 MHz  (config.py: double-passed 200 MHz tag AOM)
fwd   = f1c - twofA               # forward must sit here; retro = fwd + 2f_A = F'1 centre

print(f"rep2-on-retro  (forward at {fwd:+.1f} MHz,  retro = forward + {twofA:.0f} = {f1c:+.1f} MHz):")
print(f"  {'interpretation':<20}  {'lock':>9}  {'AOM':>8}  {'DP cell':>9}")
for name, lock in [("cooler-1240", lock_A), ("repumper+5540", lock_B)]:
    aom = fwd - lock
    dp  = aom / 2
    ok  = "  <- standard DP" if 80 < dp < 350 else "  <- non-standard"
    print(f"  {name:<20}  {lock:>+9.0f}  {aom:>+8.0f}  {dp:>7.0f} MHz{ok}")
print()

# ── forward byproduct: scatter rate and AC-Stark on the cooled state ───────────
# The forward (sigma-) beam at fwd MHz passes through the F'=1,2,3 field.
# Find the nearest dipole-allowed transition (F'=0 omitted: F=2->F'0 is delta-F=2 forbidden).
nearby = [(Fp, mp, line(Fp, mp)) for Fp in (1, 2, 3) for mp in range(-Fp, Fp + 1)]
Fp_near, _, E_near = min(nearby, key=lambda t: abs(fwd - t[2]))
D   = abs(fwd - E_near)                            # detuning from nearest allowed line (MHz)
Om  = 3.0                                          # typical repumper Rabi (2pi*MHz)
# Gamma in rad/s; Om/D in MHz -- the (2pi*1e6)^2 factors cancel in Omega^2/Delta^2
GAMMA = 2 * np.pi * 6.07e6
scat  = GAMMA * Om**2 / (4 * D**2)                # scatter rate (1/s)
# AC-Stark in kHz: (Omega/2)^2 / Delta with both in MHz gives MHz; *1e3 -> kHz
ac_kHz = Om**2 / (4 * D) * 1e3

print(f"Forward byproduct at {fwd:+.1f} MHz (sigma-):")
print(f"  nearest allowed line: F'={Fp_near},  {D:.0f} MHz off  ({D/c.Gamma:.0f} linewidths)")
print(f"  scatter rate:  {scat:.0f} /s  (< 1 per motional period at nu_z = {c.nu_z*1e3:.0f} kHz)")
print(f"  AC-Stark:      {ac_kHz:.1f} kHz = {ac_kHz / (c.nu_z * 1e3):.3f} nu_z  (servo-suppressed)\n")

# ── 2f_A sweep: raising the tag AOM shrinks the master AOM ────────────────────
# Useful if a smaller DP cell is available; each row also shifts f_mod = 6.835 + 2f_A GHz.
print(f"2f_A sweep  (lock A = {lock_A:+.0f} MHz assumed; use lock B if that is correct):")
print(f"  {'2f_A':>5}  {'forward':>8}  {'AOM':>7}  {'DP cell':>8}  f_mod")
for tag in (300, 350, 400, 450, 500, 550, 600):
    fwd_s = f1c - tag
    aom_s = fwd_s - lock_A
    dp_s  = aom_s / 2
    ok    = "ok " if 80 < dp_s < 350 else "nonstd"
    fmod  = (HFSg + tag) / 1e3          # EOM drive in GHz
    print(f"  {tag:>5}  {fwd_s:>+8.1f}  {aom_s:>+7.0f}  {dp_s:>7.0f}  {ok}  f_mod={fmod:.3f} GHz")

print("\nNote: changing twofA in config.py shifts ALL tone positions (rep1, rep2, f_mod).")
