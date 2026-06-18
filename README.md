# squinterferometer

**Build a snap-together squeezed light interferometer or other quantum-optics experiment**

Build an optical setup the way you'd sketch it on a whiteboard. Start with some input lasers or vacuum ports, then chain components one after another, adding squeezing, optical amplification, phase shifts, beam splitters, or loss.

The toolkit calculates output fields, photon numbers, homodyne measurements, noise, and phase sensitivity using exact symbolic expressions. Behind the scenes each component is represented as a matrix of creation and annihilation operators. The inputs use symbolic sympy variables and the results are given as sympy algebraic formulas. Substitute numbers for each variable to get a numerical answer.

Example interferometers
- Mach–Zehnder interferometer
- Squeezed light interferometer
- boosted squeezed light interferometer, aka SU(1,1) interferometer
- Mach–Zehnder interferometer with loss

Example answers
- how many photons land at each output port
- what a homodyne detector would read at a given phase
- how precisely can you measure a phase, i.e. the whole reason to use squeezing

## Example

```python
from squinterferometer import Squint
from sympy import symbols, pi

α = symbols('α', real=True)   # laser amplitude
ϕ = symbols('ϕ', real=True)   # the phase you're measuring

c = (Squint().input('a', α).input('b', 0)
             .bs('a', 'b').phase('a', ϕ).bs('a', 'b'))

c.dphisqr(c.N('a'), ϕ, ϕ0=pi/2)   # exact phase sensitivity Δϕ²
```

## Requirements

SymPy
