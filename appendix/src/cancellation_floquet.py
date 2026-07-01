"""
cancellation_floquet.py -- can the F'1 leak (chapter 03) be cancelled?  The Floquet verdict.

Chapter 03 showed that *if* the dark state's residual coupling onto |F'1,0> could be removed, the floor
would recover from ~0.055 toward the ~0.03 no-leak ideal -- but it tested this by SCALING the probe's
coherent F'1 edge, a proxy that quietly assumed a knob you do not physically have. This module asks the
question honestly, with a real cancellation tone at its real frequency, and finds the answer is NO: a
single co-propagating tone cannot cancel the leak.

THE SETUP (a minimal internal model, motion dropped -- the leak is an internal-coherence effect):
  0 = |2,+1>   (control leg of the dark state)
  1 = |1,-1>   (probe leg)
  2 = |F'2,0>  (the intended vertex)
  3 = |F'1,0>  (the leak vertex, 157 MHz below)
Control sigma- and probe sigma+ drive (0,1)->2 with the dark-state ratio, and ALSO (0,1)->3 with the
atomic ratios R_c=-0.346, R_p=+1.732. The dark state |D2> = (Op|0> - Oc|1>)/N is dark on 2 but leaks on 3.

THE TWO CANCELLERS (both add a |1>->|3> coupling of amplitude canceller_amp):
  * STATIC (in-frame): a hypothetical tone at the probe's OWN frame frequency, so its |3> coupling is
    time-independent. At canceller_amp = Op*(R_c - R_p) it makes |D2> dark on 3 as well -> a PERFECT dark state ->
    all scattering stops. The principle works. But it is unphysical: a real tone at the probe frequency is
    sigma+ on F=1 and drives |3> AND |2> together at the fixed atomic ratio -- it only rescales Op, it
    cannot move the F'1/F'2 ratio, and it is that ratio (R_p != R_c) that makes the leak.
  * RESONANT (the real tone): an actual laser tuned to |1,-1>->|F'1,0> (the F1->F'1 line). In the
    control/probe frame it BEATS at Delta + 157, so its |3> coupling oscillates and cannot statically null
    the (static) leak. It just adds its own off-resonant scatter, making the floor WORSE. This is the only
    physically realizable single tone.

Run `python cancellation_floquet.py` (needs numpy, qutip).
"""

import numpy as np, qutip as qt

Delta, split, Gamma = 45.0, 157.0, 6.07          # detuning, F'2-F'1 spacing, linewidth (2pi MHz)
Om_c, Om_p = 4.0, 2.0                       # representative dark-state Rabis
R_c, R_p = -0.346, 1.732                    # the atomic F'1/F'2 coupling ratios (from chapter 03)

_k = lambda i: qt.basis(4, i)
_op = lambda i, j: _k(i) * _k(j).dag()
_H0 = (Delta*_op(2, 2) + (Delta+split)*_op(3, 3)
       + Om_c*_op(2, 0) + Om_c*R_c*_op(3, 0)        # control: ->F'2 and ->F'1
       + Om_p*_op(2, 1) + Om_p*R_p*_op(3, 1))       # probe:   ->F'2 and ->F'1
_H0 = _H0 + _H0.dag()
_cops = [np.sqrt(Gamma/2)*_op(0, 2), np.sqrt(Gamma/2)*_op(1, 2),
         np.sqrt(Gamma/2)*_op(0, 3), np.sqrt(Gamma/2)*_op(1, 3)]
_Pe1, _Pe2 = _op(3, 3), _op(2, 2)
beat_freq = Delta + split                                       # the beat frequency of the real tone


def scatter(canceller_amp=0.0, resonant=True):
    """Total steady-state excited population (= total scatter / Gamma) with a canceller of amplitude canceller_amp
       on |1>->|3>.  resonant=True is the real F1->F'1 tone (beats); resonant=False is the in-frame proxy."""
    if canceller_amp == 0.0:
        rho = qt.steadystate(_H0, _cops)
        return float(qt.expect(_Pe1, rho) + qt.expect(_Pe2, rho))
    if not resonant:
        H = _H0 + canceller_amp*_op(3, 1) + canceller_amp*_op(1, 3)
        rho = qt.steadystate(H, _cops)
        return float(qt.expect(_Pe1, rho) + qt.expect(_Pe2, rho))
    H = [_H0, [canceller_amp*_op(3, 1), lambda t: np.exp(-1j*beat_freq*t)],
              [canceller_amp*_op(1, 3), lambda t: np.exp(1j*beat_freq*t)]]
    tl = np.linspace(0, 14.0, 7000)   # to t=14 (~85x the 1/Gamma relaxation time -> steady state); 7000 pts resolve the Delta+157 beat
    r = qt.mesolve(H, qt.steadystate(_H0, _cops), tl, c_ops=_cops, e_ops=[_Pe1, _Pe2])
    n = int((2*np.pi/beat_freq) / (tl[1]-tl[0])) * 4        # average over the last ~4 beat periods
    return float(np.mean(r.expect[0][-n:]) + np.mean(r.expect[1][-n:]))


if __name__ == "__main__":
    print(__doc__.split("Run ")[0])
    g_null = Om_p * (R_c - R_p)
    base = scatter(0.0)
    print("=" * 73)
    print("no canceller:                          total scatter = %.3e" % base)
    print("static-null amplitude canceller_amp = Op(R_c-R_p) = %.2f\n" % g_null)
    print("  canceller_amp            STATIC (in-frame proxy)      RESONANT (the real tone)")
    print("  " + "-" * 68)
    for f in (0.5, 1.0, 1.5):
        canceller_amp = f * g_null
        print("  %+6.2f         total = %.3e             total = %.3e"
              % (canceller_amp, scatter(canceller_amp, False), scatter(canceller_amp, True)))
    print("  " + "-" * 68)
    print("VERDICT:")
    print("  static proxy at canceller_amp=%.2f  ->  scatter ~ 0  (perfect dark state; the principle)." % g_null)
    print("  but that tone is unphysical (it only rescales Op).")
    print("  the real resonant tone -> scatter RISES.  The Delta+157 beat defeats cancellation.")
    print("  => the F'1 leak is irreducible with a single co-propagating tone.")
