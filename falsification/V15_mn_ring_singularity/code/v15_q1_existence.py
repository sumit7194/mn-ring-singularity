#!/usr/bin/env python
"""V15 Q1: existence map of the Manko-Novikov B = 0 set (equatorial ring and any off-equator zeros).

Own implementation (formulas read from the GLM 2008 / Brink 2008 transcription in V9 v9_nve.mn_functions and
V13 v13_field.py; no code imported from there).

  a = -alpha exp(-2 beta Sa),  Sa = -1 + sum_{l=0..2} (x-y) P_l(u) / R^{l+1}
  b =  alpha exp( 2 beta Sb),  Sb =  1 + sum_{l=0..2} (-1)^{3-l} (x+y) P_l(u) / R^{l+1}
  R = sqrt(x^2+y^2-1), u = x y / R,  alpha = (-1+kappa)/chi, kappa = sqrt(1-chi^2) = k/M, beta = q/kappa^3 (M = 1)
  A = (x^2-1)(1+ab)^2 - (1-y^2)(b-a)^2,   B = [x+1+(x-1)ab]^2 + [(1+y)a+(1-y)b]^2.

Parts: 1 symbolic reduction (sympy), 2 PASS checks (mpmath, 40 digits, two routes), 3 (chi, q) map,
4 off-equator analysis (exact parametrisation of the zero set + numeric scans + direct minimisation of B).
Single process; thread pools pinned to 1 (machine is shared; limit is 2 cores)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[v] = '1'
import sys, json, time, platform, random
import numpy as np, sympy as sp, mpmath as mp, scipy
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
os.makedirs(RES, exist_ok=True)
mp.mp.dps = 40
T_START = time.time()
OUT = dict(task='V15 Q1 existence map of the MN B = 0 set', date='2026-10-17',
           versions=dict(python=sys.version.split()[0], sympy=sp.__version__, mpmath=mp.__version__,
                         numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()),
           mp_dps=mp.mp.dps)


# ============================================================ Part 1: symbolic reduction
def part1():
    r = {}
    x, y, be, al, bsym = sp.symbols('x y beta alpha b', real=True)
    R = sp.sqrt(x**2 + y**2 - 1); u = x*y/R
    P = [sp.Integer(1), u, (3*u**2 - 1)/2]
    Sa = -1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))
    Sb = 1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))
    a = -al*sp.exp(-2*be*Sa); b = al*sp.exp(2*be*Sb)
    Am = lambda a, b, y: (x**2 - 1)*(1 + a*b)**2 - (1 - y**2)*(b - a)**2
    Bm = lambda a, b, y: (x + 1 + (x - 1)*a*b)**2 + ((1 + y)*a + (1 - y)*b)**2
    # (1) a(x,0) = -b(x,0): the two exponents coincide at y = 0
    Sa0, Sb0 = sp.simplify(Sa.subs(y, 0)), sp.simplify(Sb.subs(y, 0))
    r['Sa(x,0)+Sb(x,0)... (-Sa0 - Sb0 should be 0 => exponents -2beta*Sa0 == 2beta*Sb0)'] = str(sp.simplify(Sa0 + Sb0))
    r['a(x,0)+b(x,0)'] = str(sp.simplify(a.subs(y, 0) + b.subs(y, 0)))
    R0 = sp.sqrt(x**2 - 1)
    g = 1 - x/R0 + x/(2*R0**3)
    r['Sb(x,0) - (1 - x/R + x/(2R^3))'] = str(sp.simplify(Sb0 - g))
    # (2) B, A on y = 0 with a = -b  (b symbol)
    B0 = sp.factor(sp.expand(Bm(-bsym, bsym, 0)))
    A0 = sp.factor(sp.expand(Am(-bsym, bsym, 0)))
    r['B(x,0) with a=-b'] = str(B0)
    r['A(x,0) with a=-b'] = str(A0)
    bsol = sp.sqrt((x + 1)/(x - 1))
    r['A(x,0) at b^2=(x+1)/(x-1)'] = str(sp.simplify(A0.subs(bsym, bsol)))
    r['B(x,0) at b^2=(x+1)/(x-1)'] = str(sp.simplify(B0.subs(bsym, bsol)))
    # (3) x = coth s, C = cosh s = x/R:  2g = (C-1)^2 (C+2)
    C = sp.symbols('C', positive=True)
    xs_, Rs_ = C/sp.sqrt(C**2 - 1), 1/sp.sqrt(C**2 - 1)
    G = sp.factor(sp.simplify(2*g.subs(x, xs_).subs(sp.sqrt(xs_**2 - 1), Rs_)))
    G2 = sp.factor(sp.simplify(2 - 2*C + C*(C**2 - 1)))
    r['2g in terms of C=x/R'] = str(G2)
    # numeric cross-check of 2g(coth s) = (cosh s-1)^2(cosh s+2) at 20 points, 50 digits
    with mp.workdps(50):
        errs = []
        for s in [mp.mpf(i)/7 for i in range(1, 21)]:
            xx = mp.coth(s); RR = mp.sqrt(xx**2 - 1)
            errs.append(abs(2*(1 - xx/RR + xx/(2*RR**3)) - (mp.cosh(s) - 1)**2*(mp.cosh(s) + 2)))
        r['max |2g(coth s) - (cosh s-1)^2(cosh s+2)| (20 pts, 50 dig)'] = mp.nstr(max(errs), 5)
    # (4) -ln|alpha| = artanh(kappa)
    kap = sp.symbols('kappa', positive=True)
    absal = (1 - kap)/sp.sqrt(1 - kap**2)
    r['-ln|alpha| - artanh(kappa)'] = str(sp.simplify(sp.expand_log(-sp.log(absal) - sp.atanh(kap).rewrite(sp.log), force=True)))
    # (5) monotonicity lemma: d ln G/ds = 2 coth(s/2) + sinh s/(cosh s+2)
    s = sp.symbols('s', positive=True)
    Gs = (sp.cosh(s) - 1)**2*(sp.cosh(s) + 2)
    dlnG = sp.diff(Gs, s)/Gs
    r['dlnG/ds - [2coth(s/2) + sinh s/(cosh s+2)] (numeric, 10 pts)'] = float(max(abs(sp.N((dlnG - 2*sp.coth(s/2) - sp.sinh(s)/(sp.cosh(s) + 2)).subs(s, sv), 40)) for sv in [sp.Rational(i, 3) for i in range(1, 11)]))
    with mp.workdps(50):
        r['max |-ln|alpha| - artanh kappa| (20 chi, 50 dig)'] = mp.nstr(max(abs(-mp.log((1 - mp.sqrt(1 - c**2))/c) - mp.atanh(mp.sqrt(1 - c**2))) for c in [mp.mpf(i)/21 for i in range(1, 21)]), 5)
    # (6) R = 1/2 boundary: s* = asinh 2, G(s*) = 2(1+sqrt5)
    r['G(asinh 2)'] = str(sp.nsimplify(sp.simplify(Gs.subs(s, sp.asinh(2)))))
    # (7) off-equator: S = Sb - Sa, D = Sb + Sa, stable forms
    S = Sb - Sa; D = Sb + Sa
    Dst = y*(3*(x**2 - 1)*(1 - y**2)/R**2 - 2*y**2)/R**3
    Sst = (1 - y**2)*(x*(x**2 - 1 - 2*y**2)/R**5 - 2/(R*(R + x)))
    r['D (factored)'] = str(sp.factor(sp.simplify(D)))
    r['D - stable form'] = str(sp.simplify(D - Dst))
    r['S - stable form'] = str(sp.simplify(sp.radsimp(S - Sst)))
    with mp.workdps(50):
        fS = sp.lambdify((x, y), S - Sst, 'mpmath'); fD = sp.lambdify((x, y), D - Dst, 'mpmath')
        rng = random.Random(7); m = 0
        for _ in range(30):
            xx = 1 + mp.mpf(10)**mp.mpf(rng.uniform(-6, 2)); yy = mp.mpf(rng.uniform(-0.999, 0.999))
            m = max(m, abs(fS(xx, yy)), abs(fD(xx, yy)))
        r['max |S,D - stable forms| (30 random pts, 50 dig)'] = mp.nstr(m, 5)
    X = sp.symbols('X', positive=True)
    r['S as x->inf (leading)'] = str(sp.factor(sp.series(S.subs(x, 1/X), X, 0, 5).removeO()))
    r['D as x->inf (leading)'] = str(sp.factor(sp.series(D.subs(x, 1/X), X, 0, 4).removeO()))
    # (8) B = 0 (a,b real, a>0>b) implies A = 0 everywhere: substitute the two conditions
    p = sp.sqrt((x + 1)*(1 - y)/((x - 1)*(1 + y))); q_ = sp.sqrt((x + 1)*(1 + y)/((x - 1)*(1 - y)))
    with mp.workdps(50):
        fA = sp.lambdify((x, y), Am(p, -q_, y), 'mpmath'); fB = sp.lambdify((x, y), Bm(p, -q_, y), 'mpmath')
        rng = random.Random(11); m = 0
        for _ in range(30):
            xx = 1 + mp.mpf(10)**mp.mpf(rng.uniform(-6, 2)); yy = mp.mpf(rng.uniform(-0.999, 0.999))
            m = max(m, abs(fA(xx, yy)), abs(fB(xx, yy)))
        r['max |A|,|B| at a=p, b=-r solving both B-squares (30 random pts)'] = mp.nstr(m, 5)
    r['A at B-zero, by hand'] = ('ab = -(x+1)/(x-1) => (x^2-1)(1+ab)^2 = 4(x+1)/(x-1); (b-a)^2 = (p+r)^2 = '
                                 '4(x+1)/((x-1)(1-y^2)) => (1-y^2)(b-a)^2 = 4(x+1)/(x-1) => A = 0.')
    # (9) y -> 0 limit of the off-equator branch: beta -> R^3/3 ; T0(s) = tanh(s/2)(1+sech^2(s/2)/2)/3 - s
    with mp.workdps(50):
        m = 0
        for sv in [mp.mpf(i)/5 for i in range(1, 31)]:
            lhs = (mp.cosh(sv) - 1)**2*(mp.cosh(sv) + 2)/mp.sinh(sv)**3
            rhs = mp.tanh(sv/2)*(1 + mp.sech(sv/2)**2/2)
            m = max(m, abs(lhs - rhs))
        r['max |G/sinh^3 - tanh(s/2)(1+sech^2(s/2)/2)| (30 pts)'] = mp.nstr(m, 5)
    r['derivation'] = [
        'y=0: P1=0, P2=-1/2, Sa=-1+x/R-x/(2R^3) = -Sb, so a=-b identically (alpha<0 => a>0>b).',
        'B(x,0) = (x+1-(x-1)b^2)^2 ; A(x,0) = -(x+1-(x-1)b^2)((x+1)b^2-(x-1)) ; b^2 = alpha^2 exp(4 beta g).',
        'B=0 on y=0  <=>  alpha^2 exp{4 beta [1 - x/R + x/(2R^3)]} = (x+1)/(x-1)  (=> A=0 too).',
        'With s = arcoth x (x/R = cosh s, R = 1/sinh s) and -ln|alpha| = artanh(kappa):',
        '   beta (cosh s - 1)^2 (cosh s + 2) = s + artanh(kappa),  beta = q/kappa^3, kappa = sqrt(1-chi^2).',
        'Equivalently 4 beta sinh^4(s/2)(cosh s+2) = s + artanh(kappa); in x: beta (x/R-1)^2 (x/R+2) = arcoth x + artanh kappa.',
        'Uniqueness: beta = H(s) := (s+T)/G(s), T = artanh kappa > 0. H decreases strictly: (s+T) G\'/G >= s G\'/G '
        '= 2 s coth(s/2) + s sinh s/(cosh s+2) > 4 > 1. H -> inf (s->0, x->inf), H -> 0 (s->inf, x->1+).',
        '=> exactly one root x_s > 1 for every q > 0 and every chi in (0,1); none for q <= 0 (G > 0, s+T > 0).',
        'R = 1/2 <=> x = sqrt5/2 <=> s = asinh 2, G = 2(1+sqrt5): root has R < 1/2 iff q < q*(chi) = '
        'kappa^3 [asinh 2 + artanh kappa] / (2(1+sqrt5)).',
        'Off equator: with a = p > 0, b = -r < 0, B = 0 <=> pr = (x+1)/(x-1) and (1+y)p = (1-y)r <=> '
        'E1: beta S = arcoth x + artanh kappa,  E2: beta D = artanh y, with S = Sb - Sa, D = Sb + Sa.',
        'Hence every B=0 point (x,y), y != 0, fixes beta = artanh(y)/D and T = beta S - arcoth x uniquely; it is a zero for '
        'a physical chi iff T > 0 (chi = sech T, q = beta tanh^3 T). B = 0 implies A = 0 identically (shown above).']
    return r


# ============================================================ helpers (mpmath)
def params(chi, q):
    chi, q = mp.mpf(chi), mp.mpf(q)
    kap = mp.sqrt(1 - chi**2)
    return dict(chi=chi, q=q, kap=kap, al=(-1 + kap)/chi, beta=q/kap**3, T=mp.atanh(kap))


def full_ab(x, y, Pm):
    """a, b straight from the definitions (independent of the reduction)."""
    R = mp.sqrt((x - 1)*(x + 1) + y*y); u = x*y/R; P = [1, u, (3*u*u - 1)/2]
    Sa = -1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))
    Sb = 1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))
    return -Pm['al']*mp.exp(-2*Pm['beta']*Sa), Pm['al']*mp.exp(2*Pm['beta']*Sb)


def full_AB(x, y, Pm):
    a, b = full_ab(x, y, Pm)
    s1 = x + 1 + (x - 1)*a*b; s2 = (1 + y)*a + (1 - y)*b
    A = (x*x - 1)*(1 + a*b)**2 - (1 - y*y)*(b - a)**2
    scale = (x + 1)**2 + ((x - 1)*a*b)**2 + ((1 + y)*a)**2 + ((1 - y)*b)**2
    return A, s1**2 + s2**2, s1, s2, scale


def F_red(s, Pm):
    """Reduced equatorial function in s = arcoth x (root <=> B(x,0) = 0)."""
    return 4*Pm['beta']*mp.sinh(s/2)**4*(mp.cosh(s) + 2) - s - Pm['T']


def F_red_np(xm1, beta, T):
    s = 0.5*np.log1p(2/xm1)
    return 4*beta*np.sinh(s/2)**4*(np.cosh(s) + 2) - s - T


def equatorial_roots(Pm, lo=-8, hi=3, n=4000):
    """All roots x>1: sign changes of F_red on a log grid in x-1 in [1e lo, 1e hi], refined at 40 digits.
    Outside the grid no root is possible (H(s) strictly decreasing: if F has sign + at the lower end it stays +
    below it; F -> -T < 0 as x -> inf).  Returns roots and the end-point signs."""
    e = [mp.mpf(lo) + (hi - lo)*mp.mpf(i)/(n - 1) for i in range(n)]
    xm1 = [mp.mpf(10)**ei for ei in e]
    ss = [mp.acoth(1 + v) for v in xm1]
    Fv = [F_red(s, Pm) for s in ss]
    roots = []
    for i in range(n - 1):
        if Fv[i] == 0 or Fv[i]*Fv[i + 1] < 0:
            sr = mp.findroot(lambda s: F_red(s, Pm), (ss[i + 1], ss[i]), solver='anderson', tol=mp.mpf(10)**(-35))
            roots.append(mp.coth(sr))
    return roots, (int(mp.sign(Fv[0])), int(mp.sign(Fv[-1])))


# ============================================================ Part 2: PASS checks
def part2():
    r = {}
    M, aa, qq = sp.Integer(89), sp.Integer(80), sp.Rational(19, 20)
    k = sp.sqrt(M**2 - aa**2)
    betaF = sp.nsimplify(qq*M**3/k**3)
    r['F5-R k'] = str(k); r['F5-R beta = qM^3/k^3'] = str(betaF)
    r['F5-R beta == 13394411/1186380'] = bool(betaF == sp.Rational(13394411, 1186380))
    pts = {'C0': (mp.mpf('0.9'), mp.mpf('0.95'), mp.mpf('1.92982825852'), mp.mpf('1e-9')),
           'F5-R': (mp.mpf(80)/89, mp.mpf(19)/20, mp.mpf('1.92192258792'), mp.mpf('1e-9')),
           'chi0.6_q0.1024': (mp.mpf('0.6'), mp.mpf('0.1024'), mp.mpf('1.0759622'), mp.mpf('1e-6')),
           'Kerr_q0_chi0.9 (negative control)': (mp.mpf('0.9'), mp.mpf(0), None, None),
           'Kerr_q0_chi0.6 (negative control)': (mp.mpf('0.6'), mp.mpf(0), None, None)}
    for name, (chi, q, xref, tol) in pts.items():
        Pm = params(chi, q)
        roots, ends = equatorial_roots(Pm)
        d = dict(chi=mp.nstr(chi, 15), q=mp.nstr(q, 15), beta=mp.nstr(Pm['beta'], 25), n_roots_x_gt_1=len(roots),
                 F_sign_at_x_minus_1_eq_1em8=ends[0], F_sign_at_x_minus_1_eq_1e3=ends[1],
                 roots=[mp.nstr(xr, 30) for xr in roots])
        if name == 'F5-R':
            d['beta_mp_minus_exact'] = mp.nstr(Pm['beta'] - mp.mpf(13394411)/1186380, 5)
        for xr in roots:
            # route 2: unreduced first square of B at y=0, straight from the a,b definitions
            x2 = mp.findroot(lambda X: full_AB(X, mp.mpf(0), Pm)[2]/(X + 1), xr*(1 + mp.mpf('1e-6')), tol=mp.mpf(10)**(-35))
            A, B, s1, s2, sc = full_AB(xr, mp.mpf(0), Pm)
            d.setdefault('route2_unreduced_root', []).append(mp.nstr(x2, 30))
            d.setdefault('route1_minus_route2', []).append(mp.nstr(xr - x2, 5))
            d.setdefault('B/scale at root', []).append(mp.nstr(B/sc, 5))
            d.setdefault('A/scale^(1/2) at root', []).append(mp.nstr(A/mp.sqrt(sc), 5))
            d.setdefault('R at root', []).append(mp.nstr(mp.sqrt(xr**2 - 1), 15))
            By = full_AB(xr, mp.mpf('1e-3'), Pm)
            d.setdefault('B/scale at (x_s, y=1e-3) (nonzero => isolated in y)', []).append(mp.nstr(By[1]/By[4], 5))
        if xref is not None:
            dev = [abs(xr - xref) for xr in roots]
            d['reference'] = mp.nstr(xref, 15); d['tolerance'] = mp.nstr(tol, 3)
            d['|root - reference|'] = [mp.nstr(v, 5) for v in dev]
            d['PASS'] = bool(len(roots) >= 1 and min(dev) <= tol)
            d['other roots'] = len(roots) - 1
        else:
            d['PASS'] = bool(len(roots) == 0 and ends == (-1, -1))
            # explicit: (x+1)/(x-1) > 1 > alpha^2 for all x > 1
            d['alpha^2'] = mp.nstr(Pm['al']**2, 15)
        qs = Pm['kap']**3*(mp.asinh(2) + Pm['T'])/(2*(1 + mp.sqrt(5)))
        d['q*(chi) R=1/2 threshold'] = mp.nstr(qs, 15)
        r[name] = d
    r['Q1_PASS'] = bool(all(r[k_]['PASS'] for k_ in pts))
    return r


# ============================================================ off-equator (numpy, stable forms)
def SD_np(xm1, w, a_override=None):
    """S, D at x = 1+xm1, y = tanh w, stable forms.  xm1 may underflow; pass a = x^2-1 explicitly if so."""
    x = 1 + xm1
    a = xm1*(2 + xm1) if a_override is None else a_override
    y = np.tanh(w); om = 1/np.cosh(w)**2
    R2 = a + y*y; R = np.sqrt(R2)
    D = y*(3*a*om/R2 - 2*y*y)/R**3
    S = om*(x*(a - 2*y*y)/R**5 - 2/(R*(R + x)))
    return S, D


def trace_offeq(chi, nx=1500, nw=20000, log10_lo=-20000.0, log10_hi=3.0):
    """For fixed chi trace the curve {E1 = E2 = 0 eliminated: Phi(x,y) = artanh(y) S/D - arcoth x - T = 0}, y > 0.
    Parametrised by l = log10(x-1) (x-1 may be far below float range: then a = x^2-1 is set from l, arcoth from l).
    Returns arrays (l, w, beta, q) for D<0 roots, and diagnostics."""
    kap = np.sqrt(1 - chi**2); T = np.arctanh(kap)
    W = np.linspace(1e-3, 18, nw)
    ls = np.concatenate([np.linspace(log10_lo, -300, 300), np.linspace(-300, log10_hi, nx)[1:]])
    out = []; nroots = []; pos_roots = 0
    for l in ls:
        ln_e = l*np.log(10)
        xm1 = np.exp(ln_e)                       # may underflow to 0 (then a=0 to double precision; harmless)
        a = 2*xm1 + xm1*xm1
        ac = 0.5*(np.log(2 + xm1) - ln_e)        # arcoth x
        S, D = SD_np(xm1, W, a_override=a)
        Phi = W*S/D - ac - T
        idx = np.where((np.sign(Phi[:-1]) != np.sign(Phi[1:])) & (np.sign(D[:-1]) == np.sign(D[1:])))[0]
        k = 0
        for i in idx:
            f = lambda w: (lambda S, D: w*S/D - ac - T)(*SD_np(xm1, np.array(w), a_override=a))
            w0 = brentq(f, W[i], W[i + 1], xtol=1e-15, rtol=1e-15)
            S0, D0 = SD_np(xm1, np.array(w0), a_override=a)
            beta = w0/D0
            if D0 < 0:
                out.append((l, w0, beta, beta*kap**3)); k += 1
            else:
                pos_roots += 1
        nroots.append(k)
    arr = np.array(out)
    return arr, np.array(nroots), ls, pos_roots


def part4_global_beta_pos():
    """beta > 0 branch: off-equator zero for a physical chi needs T(x,y) = beta S - arcoth x > 0 with
    beta = artanh(y)/D > 0 (D > 0, which needs y^2 < 3/5).  Scan T over the D > 0 region at 30 digits."""
    r = {}
    with mp.workdps(30):
        best = (-mp.inf, None)
        n = 0
        ls = [mp.mpf(-40) + mp.mpf(46)*i/219 for i in range(220)]            # log10(x-1) in [-40, 6]
        ws = [mp.mpf(10)**(mp.mpf(-12) + 12*mp.mpf(i)/59) for i in range(60)] + \
             [mp.mpf('1.04')*i/180 for i in range(1, 181)]                   # w = artanh y, small-y log part + linear
        grid_T = []
        for l in ls:
            xm1 = mp.mpf(10)**l; x = 1 + xm1
            for w in ws:
                y = mp.tanh(w)
                R = mp.sqrt(xm1*(2 + xm1) + y*y); u = x*y/R; P = [1, u, (3*u*u - 1)/2]
                Sa = -1 + sum((x - y)*P[i]/R**(i + 1) for i in range(3))
                Sb = 1 + sum((-1)**(3 - i)*(x + y)*P[i]/R**(i + 1) for i in range(3))
                S, D = Sb - Sa, Sb + Sa
                if D <= 0:
                    continue
                n += 1
                T = w*S/D - mp.acoth(x)
                if T > best[0]:
                    best = (T, (l, w, y, w/D))
        r['points_scanned_with_D>0'] = n
        r['max_T_over_D>0_grid'] = mp.nstr(best[0], 10)
        r['argmax'] = dict(log10_xm1=mp.nstr(best[1][0], 6), y=mp.nstr(best[1][2], 6), beta=mp.nstr(best[1][3], 6))
        # local refinement from the argmax (maximise T in (l, w))
        def negT(v):
            l, w = v
            xm1 = mp.mpf(10)**l; x = 1 + xm1; y = mp.tanh(w)
            R = mp.sqrt(xm1*(2 + xm1) + y*y); u = x*y/R; P = [1, u, (3*u*u - 1)/2]
            Sa = -1 + sum((x - y)*P[i]/R**(i + 1) for i in range(3))
            Sb = 1 + sum((-1)**(3 - i)*(x + y)*P[i]/R**(i + 1) for i in range(3))
            S, D = Sb - Sa, Sb + Sa
            return float(-(w*S/D - mp.acoth(x))) if D > 0 else 1e9
        from scipy.optimize import minimize
        mres = minimize(negT, [float(best[1][0]), float(best[1][1])], method='Nelder-Mead',
                        options=dict(xatol=1e-10, fatol=1e-14, maxiter=4000))
        r['refined_max_T'] = -mres.fun; r['refined_at'] = dict(log10_xm1=mres.x[0], w=mres.x[1], y=float(np.tanh(mres.x[1])))
    r['edges (analytic)'] = [
        'y->0: beta -> R^3/3 and T -> tanh(s/2)(1+sech^2(s/2)/2)/3 - s <= s/4 - s < 0.',
        'x->inf (needs y^2<3/5): S ~ (3/4)(1-y^2)(1-5y^2)/x^4, D ~ y(3-5y^2)/x^3, so x T -> c(y)-1 with '
        'c(y) = (3/4) artanh(y)(1-y^2)(1-5y^2)/(y(3-5y^2)) < 1 (c<=0 for y^2>=1/5; c <= (3/4)(1.078)/3 otherwise).',
        'D->0+ (y^2<3/5, x^2-1 = 2y^4/(3-5y^2)): S < 0 there (checked at 30 dig, 7 pts), beta -> +inf, T -> -inf.',
        'R->0 (x->1, y->0) with D>0: beta S stays bounded or -> -inf while -arcoth x -> -inf.',
        'x->1 at fixed y>0: D -> -2 < 0 (not in the D>0 region).']
    with mp.workdps(30):
        chk = []
        for yv in ['0.05', '0.2', '0.4', '0.6', '0.7', '0.75', '0.77']:
            y = mp.mpf(yv); x = mp.sqrt(1 + 2*y**4/(3 - 5*y**2))
            R = mp.sqrt(x*x + y*y - 1); u = x*y/R; P = [1, u, (3*u*u - 1)/2]
            Sa = -1 + sum((x - y)*P[i]/R**(i + 1) for i in range(3))
            Sb = 1 + sum((-1)**(3 - i)*(x + y)*P[i]/R**(i + 1) for i in range(3))
            chk.append((yv, mp.nstr(Sb - Sa, 6)))
        r['S on the D=0 curve'] = chk
    r['verdict'] = 'no off-equator zero for any q > 0, chi in (0,1)' if (best[0] < 0 and -mres.fun < 0) else 'CHECK: T>0 found'
    return r


# ============================================================ Part 3 + 4: map
def part3_4_map():
    chis = np.concatenate([[0.001, 0.005], np.round(np.linspace(0.01, 0.99, 99), 4), [0.995, 0.999]])
    qs = np.round(np.linspace(-5, 5, 401), 6)
    xm1g = np.logspace(-14, 6, 6000)                    # wider than the required 1e-8..1e3
    nE = np.zeros((len(chis), len(qs)), int)            # equatorial roots
    cls = np.zeros((len(chis), len(qs)), int)           # 0 none, 1 eq R<0.5, 2 eq R>=0.5, 3 off-eq pair
    nO = np.zeros((len(chis), len(qs)), int)            # off-equator roots (y>0; each has a mirror at -y)
    xs_eq = np.full((len(chis), len(qs)), np.nan)
    off_loc = {}
    mism_R = 0; mism_cnt = 0; trace_diag = []
    for i, chi in enumerate(chis):
        kap = np.sqrt(1 - chi**2); T = np.arctanh(kap)
        qstar = kap**3*(np.arcsinh(2) + T)/(2*(1 + np.sqrt(5)))
        # ---- equatorial
        for j, q in enumerate(qs):
            beta = q/kap**3
            Fv = F_red_np(xm1g, beta, T)
            sc = np.where(np.sign(Fv[:-1]) != np.sign(Fv[1:]))[0]
            nE[i, j] = len(sc)
            if q > 0 and Fv[0] <= 0: mism_cnt += 1          # root would lie below 1e-14 (not expected on this grid)
            if len(sc) == 1:
                k = sc[0]
                xr = 1 + brentq(lambda v: F_red_np(v, beta, T), xm1g[k], xm1g[k + 1], xtol=1e-300, rtol=1e-15)
                xs_eq[i, j] = xr
                Rr = np.sqrt(xr*xr - 1)
                cls[i, j] = 1 if Rr < 0.5 else 2
                if (Rr < 0.5) != (q < qstar): mism_R += 1
            expected = 1 if q > 0 else 0
            if len(sc) != expected: mism_cnt += 1
        # ---- off-equator (beta < 0 branch) curve trace
        arr, nroots, ls, pos_roots = trace_offeq(chi)
        one = nroots == 1
        lq = arr[:, 0]; qq = arr[:, 3]
        dq = np.diff(qq)
        mono = bool(np.all(dq < 0))
        trace_diag.append(dict(chi=float(chi), n_x_samples=int(len(ls)), roots_per_x_hist={int(k): int(v) for k, v in zip(*np.unique(nroots, return_counts=True))},
                               x_samples_without_root_log10_xm1_min=(float(ls[nroots == 0].min()) if np.any(nroots == 0) else None),
                               q_first=float(qq[0]), q_last=float(qq[-1]), log10_xm1_last=float(lq[-1]),
                               q_monotone_decreasing_in_x=mono, D_positive_roots=int(pos_roots)))
        for j, q in enumerate(qs):
            if q >= 0: continue
            if q > qq[0]:
                # beyond the deepest sample (x-1 < 1e-20000): asymptotic tail beta ~ -y/2 -> 0-, monotone
                nO[i, j] = 1; off_loc[(i, j)] = ('tail', None, None)
                continue
            c = np.where((qq[:-1] - q)*(qq[1:] - q) <= 0)[0]
            nO[i, j] = len(c)
            if len(c):
                k = c[0]; t = (q - qq[k])/(qq[k + 1] - qq[k])
                off_loc[(i, j)] = ('interp', lq[k] + t*(lq[k + 1] - lq[k]), arr[k, 1] + t*(arr[k + 1, 1] - arr[k, 1]))
            if nO[i, j] >= 1 and nE[i, j] == 0: cls[i, j] = 3
        print(f'  chi={chi:.4f} done  eq mismatches so far {mism_cnt}, R-class mismatches {mism_R}, '
              f'trace mono={mono} q_first={qq[0]:.3g} q_last={qq[-1]:.3g}', flush=True)
    return chis, qs, nE, nO, cls, xs_eq, off_loc, dict(eq_count_mismatches_vs_analytic=mism_cnt,
                                                       R_class_mismatches_vs_closed_form=mism_R), trace_diag


def verify_offeq_mp(chis, qs, off_loc, nsample=40):
    """Refine a sample of off-equator roots from the map on the FULL (unreduced) B with mpmath (60 digits)."""
    rng = random.Random(5)
    keys = [k for k, v in off_loc.items() if v[0] == 'interp' and v[1] > -25]
    rng.shuffle(keys); out = []
    with mp.workdps(60):
        for (i, j) in keys[:nsample]:
            Pm = params(mp.mpf(str(chis[i])), mp.mpf(str(qs[j])))
            _, l0, w0 = off_loc[(i, j)]
            def eqs(l, w):
                x = 1 + mp.mpf(10)**l; y = mp.tanh(w)
                a, b = full_ab(x, y, Pm)
                return [mp.log(-a*b) - mp.log((x + 1)/(x - 1)), mp.log(-b/a) - mp.log((1 + y)/(1 - y))]
            try:
                sol = mp.findroot(eqs, (mp.mpf(l0), mp.mpf(w0)), tol=mp.mpf(10)**(-50), maxsteps=200)
                x = 1 + mp.mpf(10)**sol[0]; y = mp.tanh(sol[1])
                A, B, s1, s2, sc = full_AB(x, y, Pm)
                out.append(dict(chi=float(chis[i]), q=float(qs[j]), x_minus_1=mp.nstr(x - 1, 12), y=mp.nstr(y, 12),
                                R=mp.nstr(mp.sqrt(x*x + y*y - 1), 8), B_over_scale=mp.nstr(B/sc, 4),
                                A_over_sqrt_scale=mp.nstr(A/mp.sqrt(sc), 4),
                                start_vs_refined_dlog10x=mp.nstr(sol[0] - l0, 3)))
            except Exception as e:
                out.append(dict(chi=float(chis[i]), q=float(qs[j]), error=str(e)[:80]))
    return out


def direct_B_minimisation():
    """Independent of the S/D reduction: minimise B/scale over x>1, 0<=y<1 at several parameter points (mpmath
    grid + local refinement of every grid local minimum)."""
    pts = {'C0': ('0.9', '0.95'), 'F5-R': (str(mp.mpf(80)/89), '0.95'), 'chi0.6_q0.1024': ('0.6', '0.1024'),
           'chi0.99_q3': ('0.99', '3'), 'chi0.3_q0.02': ('0.3', '0.02'),
           'chi0.9_q-0.95': ('0.9', '-0.95'), 'chi0.6_q-1': ('0.6', '-1'), 'chi0.3_q-0.5': ('0.3', '-0.5')}
    res = {}
    with mp.workdps(30):
        L = [mp.mpf(-8) + 11*mp.mpf(i)/109 for i in range(110)]            # log10(x-1) in [-8, 3]
        Wg = [mp.mpf(0)] + [mp.mpf(10)**(mp.mpf(-4) + 3.3*mp.mpf(i)/19) for i in range(20)] + \
             [mp.mpf('0.2') + mp.mpf('5.8')*i/79 for i in range(80)]     # w = artanh y in [0, 6]
        for name, (c, q) in pts.items():
            Pm = params(c, q)
            G = [[None]*len(Wg) for _ in L]
            for a_, l in enumerate(L):
                x = 1 + mp.mpf(10)**l
                for b_, w in enumerate(Wg):
                    _, B, _, _, sc = full_AB(x, mp.tanh(w), Pm)
                    G[a_][b_] = B/sc
            mins = []
            for a_ in range(len(L)):
                for b_ in range(len(Wg)):
                    v = G[a_][b_]
                    nb = [G[a2][b2] for a2 in (a_ - 1, a_, a_ + 1) for b2 in (b_ - 1, b_, b_ + 1)
                          if 0 <= a2 < len(L) and 0 <= b2 < len(Wg) and (a2, b2) != (a_, b_)]
                    if v < mp.mpf('0.99') and all(v < u for u in nb):     # strict minima, off the B/scale ~ 1 plateau
                        mins.append((v, a_, b_))
            mins.sort()
            refined = []
            min_offeq_grid = min(G[a_][b_] for a_ in range(len(L)) for b_ in range(len(Wg)) if Wg[b_] >= mp.mpf('0.01'))
            for v, a_, b_ in mins[:8]:
                # local minimisation of log(B/scale) in (l, w), w >= 0 via w = |w|
                def f(vv):
                    l, w = vv
                    l = min(max(l, -8.0), 3.0); w = min(abs(w), 6.0)        # stay inside the scanned box
                    x = 1 + mp.mpf(10)**mp.mpf(l); y = mp.tanh(mp.mpf(w))
                    _, B, _, _, sc = full_AB(x, y, Pm)
                    return float(mp.log10(B/sc + mp.mpf(10)**(-60)))
                from scipy.optimize import minimize
                m = minimize(f, [float(L[a_]), float(Wg[b_])], method='Nelder-Mead',
                             options=dict(xatol=1e-12, fatol=1e-12, maxiter=3000))
                lr = min(max(m.x[0], -8.0), 3.0); yr = float(np.tanh(min(abs(m.x[1]), 6.0)))
                refined.append(dict(grid_value=mp.nstr(v, 4), log10_xm1=round(lr, 10), y=yr,
                                    log10_B_over_scale=round(m.fun, 3),
                                    zero=bool(m.fun < -25), at_box_edge=bool(lr in (-8.0, 3.0))))
            res[name] = dict(chi=c, q=q, n_strict_grid_local_minima_below_0p99=len(mins), refined=refined,
                             min_B_over_scale_on_grid_with_y_ge_0p01=mp.nstr(min_offeq_grid, 4),
                             zeros_found=[(d['log10_xm1'], d['y']) for d in refined if d['zero']])
    return res


def write_svg(chis, qs, cls, path):
    W, H, ml, mr, mt, mb = 760, 520, 70, 170, 30, 60
    pw, ph = W - ml - mr, H - mt - mb
    col = {0: '#ffffff', 1: '#d95f02', 2: '#1b9e77', 3: '#7570b3'}
    X = lambda c: ml + c*pw
    Y = lambda q: mt + (5 - q)/10*ph
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Helvetica,Arial,sans-serif" font-size="12">',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    ce = np.concatenate([[0], (chis[1:] + chis[:-1])/2, [1]])
    qe = np.concatenate([[qs[0]], (qs[1:] + qs[:-1])/2, [qs[-1]]])
    for i in range(len(chis)):
        j = 0
        while j < len(qs):
            k = j
            while k + 1 < len(qs) and cls[i, k + 1] == cls[i, j]: k += 1
            o.append(f'<rect x="{X(ce[i]):.2f}" y="{Y(qe[k+1]):.2f}" width="{X(ce[i+1])-X(ce[i])+0.3:.2f}" '
                     f'height="{Y(qe[j])-Y(qe[k+1])+0.3:.2f}" fill="{col[int(cls[i, j])]}"/>')
            j = k + 1
    cc = np.linspace(0.0005, 0.9995, 400); kap = np.sqrt(1 - cc**2)
    qst = kap**3*(np.arcsinh(2) + np.arctanh(kap))/(2*(1 + np.sqrt(5)))
    pts = ' '.join(f'{X(c):.2f},{Y(min(q, 5)):.2f}' for c, q in zip(cc, qst))
    o.append(f'<polyline points="{pts}" fill="none" stroke="#000" stroke-width="1.6" stroke-dasharray="5,3"/>')
    o.append(f'<line x1="{X(0)}" y1="{Y(0)}" x2="{X(1)}" y2="{Y(0)}" stroke="#000" stroke-width="1.2"/>')
    o.append(f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#000"/>')
    for c in np.linspace(0, 1, 11):
        o.append(f'<line x1="{X(c)}" y1="{mt+ph}" x2="{X(c)}" y2="{mt+ph+5}" stroke="#000"/>'
                 f'<text x="{X(c)}" y="{mt+ph+19}" text-anchor="middle">{c:.1f}</text>')
    for q in range(-5, 6):
        o.append(f'<line x1="{ml-5}" y1="{Y(q)}" x2="{ml}" y2="{Y(q)}" stroke="#000"/>'
                 f'<text x="{ml-9}" y="{Y(q)+4}" text-anchor="end">{q}</text>')
    o.append(f'<text x="{ml+pw/2}" y="{H-15}" text-anchor="middle">chi = a/M</text>')
    o.append(f'<text x="18" y="{mt+ph/2}" text-anchor="middle" transform="rotate(-90 18 {mt+ph/2})">q</text>')
    o.append(f'<text x="{ml}" y="18">Manko-Novikov B = 0 set (V15 Q1)</text>')
    leg = [(1, 'equatorial ring, R &lt; 0.5'), (2, 'equatorial ring, R &gt;= 0.5'), (3, 'off-equator ring pair (+-y)'),
           (0, 'none (q = 0, Kerr)')]
    for n, (k, t) in enumerate(leg):
        yy = mt + 10 + 24*n
        o.append(f'<rect x="{ml+pw+14}" y="{yy}" width="14" height="14" fill="{col[k]}" stroke="#000"/>'
                 f'<text x="{ml+pw+34}" y="{yy+11}">{t}</text>')
    yy = mt + 10 + 24*4
    o.append(f'<line x1="{ml+pw+14}" y1="{yy+7}" x2="{ml+pw+28}" y2="{yy+7}" stroke="#000" stroke-dasharray="5,3" stroke-width="1.6"/>'
             f'<text x="{ml+pw+34}" y="{yy+11}">q*(chi): R = 1/2</text>')
    o.append('</svg>')
    with open(path, 'w') as fh:
        fh.write('\n'.join(o))


def main():
    print('Part 1 symbolic ...', flush=True)
    OUT['part1_symbolic'] = part1()
    print('Part 2 PASS checks ...', flush=True)
    OUT['part2_pass_checks'] = part2()
    print(json.dumps({k: OUT['part2_pass_checks'][k].get('PASS') if isinstance(OUT['part2_pass_checks'][k], dict) else OUT['part2_pass_checks'][k]
                      for k in OUT['part2_pass_checks']}), flush=True)
    print('Part 4a global beta>0 scan ...', flush=True)
    OUT['part4_offequator_beta_pos'] = part4_global_beta_pos()
    print(OUT['part4_offequator_beta_pos']['verdict'], OUT['part4_offequator_beta_pos']['max_T_over_D>0_grid'],
          OUT['part4_offequator_beta_pos']['refined_max_T'], flush=True)
    print('Part 3/4 map ...', flush=True)
    chis, qs, nE, nO, cls, xs_eq, off_loc, mism, tdiag = part3_4_map()
    pos = qs > 0; neg = qs < 0; zer = qs == 0
    summary = dict(
        chi_grid=dict(n=len(chis), values=chis.tolist()), q_grid=dict(n=len(qs), lo=float(qs[0]), hi=float(qs[-1]), step=0.025),
        x_grid_for_equatorial_bracketing='6000 log points, x-1 in [1e-14, 1e6]', **mism,
        equatorial_roots_count_values_q_gt_0=sorted(set(nE[:, pos].ravel().tolist())),
        equatorial_roots_count_values_q_le_0=sorted(set(nE[:, ~pos].ravel().tolist())),
        offeq_roots_count_values_q_lt_0=sorted(set(nO[:, neg].ravel().tolist())),
        offeq_roots_count_values_q_ge_0=sorted(set(nO[:, ~neg].ravel().tolist())),
        offeq_points_from_asymptotic_tail=int(sum(1 for v in off_loc.values() if v[0] == 'tail')),
        class_counts={int(k): int(v) for k, v in zip(*np.unique(cls, return_counts=True))},
        class_legend={0: 'none', 1: 'equatorial ring R<0.5', 2: 'equatorial ring R>=0.5', 3: 'off-equator ring pair'},
        boundary_formulas=dict(existence_equatorial='q > 0 (all chi in (0,1)); exactly one root',
                               R_half='q*(chi) = (1-chi^2)^(3/2) [asinh 2 + artanh sqrt(1-chi^2)] / (2(1+sqrt5)); R<1/2 iff q<q*',
                               offequator='q < 0 (all chi in (0,1)); exactly one pair (x_off, +-y_off) per (chi,q) on the map',
                               kerr='q = 0: no B=0 point with x>1'),
        q_star_at=dict({f'{c}': float((1 - c*c)**1.5*(np.arcsinh(2) + np.arctanh(np.sqrt(1 - c*c)))/(2*(1 + np.sqrt(5)))) for c in [0.1, 0.3, 0.5, 0.6, 0.9, 80/89, 0.99]}),
        offeq_trace_diagnostics=tdiag)
    i09 = int(np.argmin(abs(chis - 0.9)))
    summary['chi0.9_equatorial_x_s_by_q'] = {str(float(q)): (float(xs_eq[i09, j]) if not np.isnan(xs_eq[i09, j]) else None)
                                            for j, q in enumerate(qs) if abs(q*4 - round(q*4)) < 1e-9}
    # sample of off-equator locations on the map
    smp = {}
    for c in [0.1, 0.3, 0.6, 0.9, 0.99]:
        i = int(np.argmin(abs(chis - c)))
        for q in [-0.025, -0.1, -0.5, -1, -2, -5]:
            j = int(np.argmin(abs(qs - q)))
            v = off_loc.get((i, j))
            if v is None: smp[f'chi{c}_q{q}'] = None
            elif v[0] == 'tail': smp[f'chi{c}_q{q}'] = 'x-1 < 1e-20000 (asymptotic tail)'
            else:
                yv = float(np.tanh(v[2])); xm = 10.0**float(v[1])
                smp[f'chi{c}_q{q}'] = dict(log10_x_minus_1=round(float(v[1]), 4), y=round(yv, 6),
                                           R=round(float(np.sqrt(xm*(2 + xm) + yv*yv)), 6))
        summary['offeq_location_samples'] = smp
    OUT['part3_map'] = summary
    np.savez_compressed(os.path.join(RES, 'v15_q1_map_arrays.npz'), chis=chis, qs=qs, nE=nE, nO=nO, cls=cls, xs_eq=xs_eq)
    print('Part 4b mp verification of off-equator roots ...', flush=True)
    OUT['part4_offequator_mp_verification'] = verify_offeq_mp(chis, qs, off_loc)
    print('Part 4c direct B minimisation ...', flush=True)
    OUT['part4_direct_B_minimisation'] = direct_B_minimisation()
    write_svg(chis, qs, cls, os.path.join(RES, 'v15_q1_map.svg'))
    OUT['runtime_s'] = round(time.time() - T_START, 1)
    with open(os.path.join(RES, 'v15_q1_existence.json'), 'w') as fh:
        json.dump(OUT, fh, indent=1, default=str)
    print('done', OUT['runtime_s'], 's', flush=True)


if __name__ == '__main__':
    main()
