# squinter

**Build and model a squeezed light interferometer**

A toolkit for modeling conventional and squeezed-light interferometers and other Gaussian quantum-optics experiments. Build an optical setup the way you'd sketch it on a whiteboard. Start with input lasers or vacuum ports, then snap together a chain of components—squeezers, parametric amplifiers, phase shifters, beam splitters, and loss channels. Computations include the resulting output fields, photon numbers, homodyne measurements, noise, and phase sensitivity.

Behind the scenes each component is represented as a transfer matrix acting on the creation and annihilation operators in the Heisenberg picture. The inputs use symbolic sympy variables and the results are given as sympy algebraic formulas. Substitute numbers for each variable to get a numerical answer.

Example interferometers
- Mach–Zehnder interferometer
- squeezed light interferometer
- boosted squeezed light interferometer, aka SU(1,1) interferometer
- Mach–Zehnder interferometer with loss

Example answers
- how many photons land at each output port
- what a homodyne detector would read at a given phase
- how precisely can you measure the phase, i.e. the whole reason to use squeezing

## Example

```python
from squinter import Squinter
from sympy import symbols, pi

α = symbols('α', real=True)   # laser amplitude
ϕ = symbols('ϕ', real=True)   # the phase you're measuring

c = Squinter().input('a', α).input('b', 0).beamsplitter('a', 'b').phase('a', ϕ).beamsplitter('a', 'b')

c.dphisqr(c.N('a'), ϕ, ϕ0=pi/2)   # exact phase sensitivity Δϕ²
```

## Requirements

SymPy
