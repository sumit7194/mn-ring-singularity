"""V15 Q4: Komar mass of the MN equatorial ring (pre-registered; see PREREGISTRATION.md Q4)."""
import sys, os
import sympy as sp, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../V13_mn_chaotic_sea/code'))
sys.path.insert(0, os.path.join(HERE, '../../V9_mn_obstruction_check/code'))
from v13_field import build
from v9_nve import x, y
mp.mp.dps = 30

def integrand_funcs(M, a, beta):
    S = build(M, a, beta)
    k, al, psi, G, A, B, C, Q, w0 = (S[n] for n in ('k', 'al', 'psi', 'G', 'A', 'B', 'C', 'Q', 'w0'))
    gtt = -sp.exp(2*psi)*A/B
    gtp = 2*k*C/B - w0*sp.exp(2*psi)*A/B                       # f*omega, regular at A = 0
    gpp = (sp.exp(-2*psi)*k**2*Q + 4*k*w0*C - w0**2*sp.exp(2*psi)*A)/B
    Fx = 2*(gpp*sp.diff(gtt, x) - gtp*sp.diff(gtp, x))/(k*(1 - y**2))      # (*F)_{y phi}
    Fy = -2*(gpp*sp.diff(gtt, y) - gtp*sp.diff(gtp, y))/(k*(x**2 - 1))     # (*F)_{x phi}
    L = lambda e: sp.lambdify((x, y), e, 'mpmath')
    return L(Fy), L(Fx), (gtt, gtp, gpp)

def komar(fxphi, fyphi, curve, dcurve, s0, s1):
    """Komar mass: -(1/16pi) * 2pi * closed line integral of (*F)_{x phi} dx + (*F)_{y phi} dy along curve(s)."""
    def g(s):
        X, Y = curve(s); dX, dY = dcurve(s)
        return fxphi(X, Y)*dX + fyphi(X, Y)*dY
    return -mp.quad(g, [s0, (s0 + s1)/2, s1])/8      # standard Komar: -(1/8pi) oint *nabla xi = -(1/16pi) oint *d xi (see Q4 note)

def circle(xc, yc, r):
    return (lambda s: (xc + r*mp.cos(s), yc + r*mp.sin(s)), lambda s: (-r*mp.sin(s), r*mp.cos(s)))

if __name__ == '__main__':
    chi, q = sp.Rational(9, 10), sp.Rational(95, 100)
    kap = sp.sqrt(1 - chi**2); beta = q/kap**3
    out = {}
    for name, b in (('Kerr (beta=0)', 0), ('MN q=0.95', beta)):
        fxp, fyp, _ = integrand_funcs(1, chi, b)
        X = mp.mpf(200)
        asym = komar(fxp, fyp, lambda s: (X, s), lambda s: (0, 1), -1, 1)
        print(f"C1 {name}: asymptotic Komar mass = {mp.nstr(asym, 15)}")
        out[name] = (fxp, fyp, asym)
    fxp, fyp, asym = out['MN q=0.95']
    reg = komar(fxp, fyp, *circle(mp.mpf(3), mp.mpf('0.3'), mp.mpf('0.05')), 0, 2*mp.pi)
    print("C2 regular-point circle (x=3, y=0.3, r=0.05):", mp.nstr(reg, 6))
    # ring location: B = A = 0 on y = 0 (V15 Q1); refine from the email value
    _, _, (gtt, gtp, gpp) = integrand_funcs(1, chi, beta)
    S = build(1, chi, beta); Bf = sp.lambdify(x, S['B'].subs(y, 0), 'mpmath')
    xs = mp.findroot(lambda t: Bf(t), mp.mpf('1.92983'))
    print("ring x_s =", mp.nstr(xs, 20))
    for r in ('1e-2', '3e-3', '1e-3'):
        m = komar(fxp, fyp, *circle(xs, mp.mpf(0), mp.mpf(r)), 0, 2*mp.pi)
        print(f"ring circle r={r}: Komar mass = {mp.nstr(m, 12)}")
