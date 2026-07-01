"""
config.py -- chapter 03 reuses chapter 02's operating point.

The F'1-leak study runs on exactly the same 87Rb parameters as chapter 02 (same trap,
detuning, Rabis, magnetic field). Rather than copy the numbers here -- where they could
drift out of sync with chapter 02 -- this file LOADS them from 02_multilevel/src/config.py.
So every knob lives in one place; read it there.

To try a chapter-03-only value, set it in the OVERRIDES block at the bottom. It then takes
effect everywhere chapter 03 reuses the chapter-02 solver (cooling_multilevel.solve).
"""
import os

_ch02_config = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "..", "02_multilevel", "src", "config.py")
exec(open(_ch02_config).read())          # pull in chapter 02's numbers verbatim

# --- chapter-03 overrides (none by default) -------------------------------------------
# e.g.  Delta = 30.0     # a chapter-03-only single-photon detuning, to explore
