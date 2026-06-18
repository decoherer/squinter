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




#

Included examples cover a Mach–Zehnder interferometer, a two-amplifier interferometer, squeezed-vacuum injection, and interferometers with internal loss.


Then ask the questions interferometers are actually built to answer:
how many photons land at each output port
what a homodyne detector would read at a given phase
how precisely you can measure a phase — the sensitivity Δϕ², which is the whole reason you put squeezing into an interferometer in the first place

New light paths and loss channels can be added anywhere in a setup. 

Squinterferometer is intended for researchers and students who understand photon interferometry and squeezing but want to model experiments without implementing the underlying operator calculations from scratch.




# squinterferometer

squinterferometer is a Python toolkit for building and analyzing photon interferometers with coherent light, vacuum, squeezing, optical amplification, phase shifts, beam splitters, and loss. Experiments are assembled one optical element at a time, making it easy to explore both familiar interferometers and squeezing-enhanced designs.

The toolkit calculates output fields, photon numbers, homodyne measurements, averages, noise, and phase sensitivity using exact symbolic expressions. New light paths and loss channels can be added anywhere in a setup, and built-in checks help catch transformations that do not obey the required quantum rules.

Included examples cover a Mach–Zehnder interferometer, a two-amplifier interferometer, squeezed-vacuum injection, and interferometers with internal loss. squinterferometer is intended for researchers and students who understand photon interferometry and squeezing but want to model experiments without implementing the underlying operator calculations from scratch.

One-line repository description

Symbolic Python tools for photon interferometers with squeezing, amplification, loss, and phase-sensitivity calculations.






# squinterferometer

**Snap together a quantum-optics experiment and get the exact answer — no numerical simulator, no approximations.**

`squinterferometer` lets you build an optical setup the way you'd sketch it on a whiteboard: start with some input beams, then chain components one after another — beam splitters, phase shifters, squeezers, parametric amplifiers, and loss. Behind the scenes each component is a matrix and the whole experiment is their product, so the result comes back as a real algebraic formula instead of a number out of a numerical simulation.

It all runs on one simple idea: everything is a component, not a special input. Your inputs are just laser light or vacuum. Squeezing isn't a fancy state you have to prepare — it's a squeezer you drop into the line. Loss isn't a separate model — it's a beam splitter that lets some light leak away where you don't look. You can add new beams partway through, and a one-line check confirms that whatever you've wired up is physically allowed.

Then ask the questions interferometers are actually built to answer:

- how many photons land at each output port
- what a homodyne detector would read at a given phase
- how precisely you can measure a phase — the sensitivity `Δϕ²`, which is the whole reason you put squeezing into an interferometer in the first place

## Worked examples, a few lines each

- an ordinary Mach–Zehnder interferometer
- the squeezed-light "boosted" interferometer (two parametric amplifiers back to back) that beats the shot-noise limit
- the classic squeezed-vacuum trick: feed squeezed vacuum into the dark port
- a leaky interferometer that shows exactly how much loss eats your advantage

## Example

```python
from squinterferometer import Squint
from sympy import symbols, pi

α = symbols('α', real=True)   # laser amplitude
ϕ = symbols('ϕ', real=True)   # the phase you're measuring

c = (Squint().input('a', α).input('b', 0)
             .bs('a', 'b').phase('a', ϕ).bs('a', 'b'))

c.dphisqr(c.N('a'), ϕ, ϕ0=pi/2)   # exact phase sensitivity Δϕ²
assert c.check().is_zero_matrix    # is this circuit physical? yes
```

## Requirements

Just SymPy. There's no numerical backend — every answer is an exact formula you can read, differentiate, or simplify by hand.


