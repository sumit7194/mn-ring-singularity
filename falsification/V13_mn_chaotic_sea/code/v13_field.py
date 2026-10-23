"""V13 reduced 2-DOF field for MN timelike geodesics, regularised at the ergosurface A = 0 (amendment A2).
H = 1/2 [g^xx p_x^2 + g^yy p_y^2 + W(x, y)] = -1/2,  W = g^{ab} p_a p_b over (t, phi) with p_t = -E, p_phi = L.
Derivation (the author's own, from the V9 transcription v9_nve.mn_functions):
  W = f (L - w E)^2 / rho^2 - E^2 / f,   rho^2 = k^2 (x^2-1)(1-y^2),   f = e^{2psi} A/B,   w = 2k e^{-2psi} C/A - w0.
  With L' = L + w0 E and the exact identity (x^2-1)(1-y^2) B^2 - 4 C^2 = A Q (verified by polynomial division):
  W = [ e^{2psi} A L'^2 - 4 k C E L' - e^{-2psi} k^2 E^2 Q ] / (rho^2 B)                  (no 1/A)
  g^xx = e^{2psi} (x^2-1) (1-al^2)^2 / (k^2 B e^G),  g^yy = e^{2psi} (1-y^2) (1-al^2)^2 / (k^2 B e^G)   (A cancels)
  g_phiphi = [ e^{-2psi} k^2 Q (x^2-1)(1-y^2)... ] -- computed as -f w^2 + rho^2/f in regularised form below.
where e^{2gamma} = (x^2-1)/(x^2-y^2) e^G A/((x^2-1)(1-al^2)^2), as in V9.
Gate A2 (prereg): regularised H equals the unregularised H to <= 1e-12 relative at 50 random points with |A| > 0.1."""
import sys, os, math
import numpy as np, sympy as sp, numba
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../V9_mn_obstruction_check/code'))
from v9_nve import mn_functions, x, y
a_s, b_s = sp.symbols('a_s b_s')

def _polys():
    Am = (x**2 - 1)*(1 + a_s*b_s)**2 - (1 - y**2)*(b_s - a_s)**2
    Bm = (x + 1 + (x - 1)*a_s*b_s)**2 + ((1 + y)*a_s + (1 - y)*b_s)**2
    Cm = (x**2 - 1)*(1 + a_s*b_s)*(b_s - a_s - y*(a_s + b_s)) + (1 - y**2)*(b_s - a_s)*(1 + a_s*b_s + x*(1 - a_s*b_s))
    num = sp.expand((x**2 - 1)*(1 - y**2)*Bm**2 - 4*Cm**2)
    Q, r = sp.div(sp.Poly(num, x, y, a_s, b_s), sp.Poly(Am, x, y, a_s, b_s))
    assert r.is_zero, 'identity failed'
    return Am, Bm, Cm, Q.as_expr()

def build(M, a_, beta):
    """Symbolic pieces (regularised and unregularised) for one MN point; exact parameters."""
    M, a_, beta = map(sp.nsimplify, (M, a_, beta))
    chi = a_/M; al = sp.nsimplify((-1 + sp.sqrt(1 - chi**2))/chi); k = sp.sqrt(M**2 - a_**2)
    R = x*sp.sqrt(1 + (y**2 - 1)/x**2); u = x*y/R
    P = [sp.Integer(1), u, (3*u**2 - 1)/2, (5*u**3 - 3*u)/2]
    psi = beta*P[2]/R**3
    A_ = -al*sp.exp(-2*beta*(-1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))))
    B_ = al*sp.exp(2*beta*(1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))))
    G = 3*beta**2*(P[3]**2 - P[2]**2)/R**6 + 2*beta*sum((x - y + (-1)**l*(x + y))*P[l]/R**(l + 1) - 2 for l in range(3)) + 8*beta
    Am, Bm, Cm, Q = _polys()
    sub = {a_s: A_, b_s: B_}
    Ax, Bx, Cx, Qx = [e.subs(sub) for e in (Am, Bm, Cm, Q)]
    w0 = 4*k*al/(1 - al**2)
    rho2 = k**2*(x**2 - 1)*(1 - y**2)
    F = mn_functions(M, a_, beta)
    return dict(k=k, al=al, psi=psi, G=G, A=Ax, B=Bx, C=Cx, Q=Qx, w0=w0, rho2=rho2, F=F)

def W_reg(S, E, L):
    Lp = L + S['w0']*E
    return (sp.exp(2*S['psi'])*S['A']*Lp**2 - 4*S['k']*S['C']*E*Lp - sp.exp(-2*S['psi'])*S['k']**2*E**2*S['Q'])/(S['rho2']*S['B'])

def ginv_reg(S):
    pref = sp.exp(2*S['psi'])*(1 - S['al']**2)**2/(S['k']**2*S['B']*sp.exp(S['G']))
    return pref*(x**2 - 1), pref*(1 - y**2)

def gpp_reg(S):
    """g_phiphi = -f w^2 + rho^2/f, regularised: f w^2 = e^{2psi}A/B (2k e^{-2psi}C/A - w0)^2
    = [4k^2 e^{-2psi} C^2/A - 4 k w0 C + w0^2 e^{2psi} A]/B, and rho^2/f = rho^2 B e^{-2psi}/A.
    Using rho^2 B^2 - 4k^2 C^2 = k^2 A Q:  g_phiphi = [e^{-2psi} k^2 Q + 4 k w0 C - w0^2 e^{2psi} A]/B."""
    return (sp.exp(-2*S['psi'])*S['k']**2*S['Q'] + 4*S['k']*S['w0']*S['C'] - S['w0']**2*sp.exp(2*S['psi'])*S['A'])/S['B']

def unreg(S, E, L):
    f, om, e2g, k = S['F']['f'], S['F']['om'], S['F']['e2g'], S['F']['k']
    rho2 = S['rho2']
    W = f*(L - om*E)**2/rho2 - E**2/f
    gxx = f*(x**2 - 1)/(k**2*e2g*(x**2 - y**2)); gyy = f*(1 - y**2)/(k**2*e2g*(x**2 - y**2))
    gpp = -f*om**2 + rho2/f
    return W, gxx, gyy, gpp

def gate_A2(name, M, a_, beta, E, L, n=50, seed=13):
    """Regularised vs unregularised at random (x, y) with |A| > 0.1, mpmath 40 digits."""
    import mpmath as mp
    mp.mp.dps = 40
    S = build(M, a_, beta); E, L = sp.nsimplify(E), sp.nsimplify(L)
    reg = [W_reg(S, E, L), *ginv_reg(S), gpp_reg(S)]; ur = unreg(S, E, L)
    fr = [sp.lambdify((x, y), e, 'mpmath') for e in reg]; fu = [sp.lambdify((x, y), e, 'mpmath') for e in ur]
    fA = sp.lambdify((x, y), S['A'], 'mpmath')
    rng = np.random.default_rng(seed); worst = [0.0]*4; used = 0
    while used < n:
        xv = mp.mpf(1.05 + 10*rng.random()); yv = mp.mpf(2*rng.random() - 1)*mp.mpf('0.98')
        if abs(fA(xv, yv)) <= 0.1: continue
        used += 1
        for i in range(4):
            r, u = fr[i](xv, yv), fu[i](xv, yv)
            worst[i] = max(worst[i], float(abs(r - u)/max(abs(u), mp.mpf('1e-30'))))
    return dict(point=name, E=str(E), L=str(L), n=n, worst_rel=dict(W=worst[0], gxx=worst[1], gyy=worst[2], gpp=worst[3]),
                PASS=all(w <= 1e-12 for w in worst))

if __name__ == '__main__':
    import json
    R_ = sp.Rational
    pts = {'C0': (1, R_(9, 10), R_(19, 20)/(R_(19, 100))**R_(3, 2), R_(19, 20), 3),
           'F5R': (89, 80, R_(13394411, 1186380), R_(19, 20), 3*89),
           'p1': (5, 3, R_(1, 5), R_(37, 40), R_(29, 10)*5),
           'p2': (13, 5, R_(-1, 3), R_(9, 10), R_(13, 10)*13)}
    out = []
    for nm, (M, a_, beta, E, Lz) in pts.items():
        res = gate_A2(nm, M, a_, beta, E, Lz); out.append(res); print(json.dumps(res), flush=True)
    json.dump(out, open(os.path.join(HERE, '../results/v13_gate_A2.json'), 'w'), indent=1)
    print('ALL PASS' if all(o['PASS'] for o in out) else 'GATE A2 FAILED')
