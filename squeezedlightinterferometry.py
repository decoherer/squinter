
from sympy import sqrt,exp,cos,sin,atan,atan2,sinh,cosh,tanh
from sympy import pi,oo,zoo,I,E,re,im,conjugate,sign,arg,Abs,Min,Max
from sympy import pprint,S,Eq,diff,symbols,Matrix,limit,sympify,expand,factor,simplify,trigsimp,powsimp,exptrigsimp,expand_trig
from sympy import Rational as R
from sympy.physics.quantum.boson import BosonOp, BosonFockKet, BosonFockBra
from sympy.physics.quantum.operatorordering import normal_ordered_form
from sympy.physics.quantum import qapply, Dagger

def nof(expr):
    return normal_ordered_form(expr.expand(),independent=True) # wrong result if not expanded! e.g. (a+b)*(at+bt) is unchanged
def mziexample(α=None,β=None,ϕ0=None): # α,β = input amplitudes, ϕ0 = phase shift
    α = α if α is not None else symbols('α',real=True)
    β = β if β is not None else symbols('β',real=True)
    ϕ = symbols('ϕ',real=True)
    αc,βc = α.conjugate(),β.conjugate()
    a,b = BosonOp('a'),BosonOp('b') # a,b are two input channels, a=top path, b=bottom path
    at,bt = a.adjoint(),b.adjoint()
    BS = sqrt(R(1,2))*Matrix([[1,-1],[1,1]]) # 50/50 beam splitter
    PH = Matrix([[exp(I*ϕ),0],[0,1]]) # phase shifters
    v = Matrix([a,b]) # input state
    A,B = BS*PH*BS*v # output state
    def ev(expr,subs={a:α,at:αc,b:β,bt:βc}): # expectation value
        return nof(expr).subs(subs)
    # evaluate phase noise at output port A
    N = ev(Dagger(A)*A) # print('N',N)
    NN = ev(Dagger(A)*A*Dagger(A)*A) # print('NN',NN)
    ΔNΔN = NN - N**2 # print('ΔN²',ΔNΔN)
    dNdϕ = diff(N,ϕ) # print('dN/dϕ',dNdϕ)
    ΔϕΔϕ = (ΔNΔN/dNdϕ**2).rewrite(exp,cos,sin).expand().simplify() # print('Δϕ²',ΔϕΔϕ)
    ΔϕΔϕ = ΔϕΔϕ if ϕ0 is None else ΔϕΔϕ.subs(ϕ,ϕ0)
    print('Δϕ² = ΔN²/(dN/dϕ)²')
    pprint(ΔϕΔϕ)
    return ΔϕΔϕ
def boostedexample(α=None,β=None,r=None,ϕ0=None): # coherent light boosted SU(1,1) interferometer
    α = α if α is not None else symbols('α',real=True)
    β = β if β is not None else symbols('β',real=True)
    r = r if r is not None else symbols('r',positive=True)
    ϕ = symbols('ϕ',real=True)
    αc,βc = α.conjugate(),β.conjugate()
    a,b = BosonOp('a'),BosonOp('b') # a,b are two input channels, a=signal input, b=idler input
    at,bt = a.adjoint(),b.adjoint()
    µ,𝜈 = cosh(r),sinh(r)
    v = Matrix([a,at,b,bt]) # input state
    OPA1 = Matrix([[µ,0,0,+𝜈],[0,µ,+𝜈,0],[0,+𝜈,µ,0],[+𝜈,0,0,µ]])
    PH = Matrix([[exp(I*ϕ),0,0,0],[0,exp(-I*ϕ),0,0],[0,0,1,0],[0,0,0,1]])
    OPA2 = Matrix([[µ,0,0,-𝜈],[0,µ,-𝜈,0],[0,-𝜈,µ,0],[-𝜈,0,0,µ]])
    A,At,B,Bt = OPA2*PH*OPA1*v
    def ev(expr,subs={a:α,at:αc,b:β,bt:βc}): # expectation value
        return nof(expr).subs(subs)
    N = ev(At*A+Bt*B) # print('N',N.subs(ϕ,ϕ0))
    dNdϕ = diff(N,ϕ) # print('dN/dϕ',dNdϕ.subs(ϕ,ϕ0))
    NN = ev((At*A+Bt*B)**2) # print('NN',NN.subs(ϕ,ϕ0))
    ΔNΔN = NN - N**2 # print('ΔN²',ΔNΔN.subs(ϕ,ϕ0))
    ΔϕΔϕ = (ΔNΔN/dNdϕ**2).rewrite(exp,cos,sin).expand()#.simplify() # 
    ΔϕΔϕ = ΔϕΔϕ if ϕ0 is None else ΔϕΔϕ.subs(ϕ,ϕ0)
    print('Δϕ² = <ΔN²>/(d<N>/dϕ)² =',ΔϕΔϕ.subs(α,0).subs(β,0).limit(ϕ,0)) # equals 1/sinh²(2r) ref:Plick15,pg5
    return ΔϕΔϕ
    # SU(1,1) boosted interferometer notes
    #  - phase of the SU(1,1) interferometer depends on sum of the phases φs + φi in each path
    #    the signal and idler can share the same modulation path, halves Vπ
    #  - can also phase modulate the pump alone. or all three with pump dominating but probably not efficiently. net phase is ϕs+ϕi-ϕp
    #  - R = 1-η + η*exp(±2√(P ηSHG))
    #    ηSHG from R fit in many papers was usually much less than measured ηSHG (Kashiwazaki20 excepted)
    #    (due to multimode factor of 2?)
    #  - Internal loss has same effect as in squeezing, but external loss ~no little effect
    #  - but boosted has exponentially better phase sensitivity than squeezing alone

if __name__ == '__main__':
    mziexample()
    mziexample(α=symbols('α'),β=0,ϕ0=pi/2)
    boostedexample()

