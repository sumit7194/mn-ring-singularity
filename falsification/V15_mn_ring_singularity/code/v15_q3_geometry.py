"""V15 Q3 -- geometry of the Manko-Novikov equatorial ring B = A = 0 (ergosurface, CTC boundary, pole of g_phiphi).

Own analysis code (does not import the V9/V13 modules; formulas re-typed from GLM 2008 arXiv:0708.0628 Eqs 2-7 /
Brink 2008 arXiv:0807.1179 App. C, cross-checked against the transcription in V9 v9_nve.py and V13 v13_field.py).

Notation: small a, b are the MN functions a(x,y), b(x,y); capital A, B, C are the GLM polynomials in (x, y, a, b);
  f = e^{2psi} A/B,  omega = 2k e^{-2psi} C/A - w0,  w0 = 4k alpha/(1-alpha^2),  rho^2 = k^2 (x^2-1)(1-y^2),
  g_phiphi = -f omega^2 + rho^2/f = N/B,  N = e^{-2psi} k^2 Q + 4k w0 C - w0^2 e^{2psi} A,
  with Q the polynomial quotient (x^2-1)(1-y^2) B^2 - 4 C^2 = A Q (re-derived here by polynomial division).
Weyl: rho = k sqrt(x^2-1) sqrt(1-y^2), z = k x y.

Sections
  S1  symbolic, generic: equatorial factorisation of A, B, C, Q and N at y = 0, b = -a; parity; transversality.
  S2  per point (C0, F5-R): ring refined, second-order Taylor of A, B, C, Q, N (two implementations: mpmath numeric
      differentiation and sympy symbolic differentiation), compared with the S1 closed forms.
  S3  rays from the ring in Weyl (rho, z): scaling of g_phiphi, f, omega; regularised vs unregularised g_phiphi;
      g_phiphi / (-omega_ring^2 f) -> 1.
  S4  local curves A = 0 and N = 0 through the ring (exact roots at small y vs quadratic prediction).
  S5  C0 global map in (rho, z): sign classes, CTC components, extents (refined), boundedness checks, PNG/SVG.
  S6  Kerr positive control (beta = 0): same algebra at the Kerr ring x = -M/k, y = 0 vs Boyer-Lindquist closed form.
Run: OMP_NUM_THREADS=2 python v15_q3_geometry.py   (2 worker processes for S5)
"""
import os, sys, json, time, platform, zlib, struct, base64
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np, sympy as sp, mpmath as mp, scipy
from scipy import ndimage
from multiprocessing import get_context

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '../results'))
xs_, ys_, as_, bs_ = sp.symbols('x y a b', real=True)
R_ = sp.Rational

# ------------------------------------------------------------------------------------------------ polynomials
def polys():
    x, y, a, b = xs_, ys_, as_, bs_
    A = (x**2 - 1)*(1 + a*b)**2 - (1 - y**2)*(b - a)**2
    B = (x + 1 + (x - 1)*a*b)**2 + ((1 + y)*a + (1 - y)*b)**2
    C = (x**2 - 1)*(1 + a*b)*(b - a - y*(a + b)) + (1 - y**2)*(b - a)*(1 + a*b + x*(1 - a*b))
    num = sp.expand((x**2 - 1)*(1 - y**2)*B**2 - 4*C**2)
    Qp, rem = sp.div(sp.Poly(num, x, y, a, b), sp.Poly(A, x, y, a, b))
    assert rem.is_zero, 'A does not divide rho^2 B^2 - 4 C^2'
    return A, B, C, Qp.as_expr()

_POLY = None
def poly_funcs():
    global _POLY
    if _POLY is None:
        A, B, C, Q = polys()
        S1 = xs_ + 1 + (xs_ - 1)*as_*bs_; S2 = (1 + ys_)*as_ + (1 - ys_)*bs_
        _POLY = {nm: sp.lambdify((xs_, ys_, as_, bs_), e, 'mpmath') for nm, e in
                 dict(A=A, B=B, C=C, Q=Q, S1=S1, S2=S2).items()}
    return _POLY

# ------------------------------------------------------------------------------------------------ S1 symbolic
def section_S1():
    x, y, a = xs_, ys_, as_
    A, B, C, Q = polys()
    out = {'Q_division_remainder_zero': True, 'Q_n_terms': len(sp.Add.make_args(sp.expand(Q)))}
    eq = {ys_: 0, bs_: -a}
    F = a**2*(x - 1) - (x + 1)                         # = (x-1)(s - s0),  s = a^2,  s0 = (x+1)/(x-1)
    claims = {'A_eq': F*(a**2*(x + 1) - (x - 1)), 'B_eq': F**2, 'C_eq': 2*a*x*F, 'Q_eq': F*(a**2*(x - 1)**3 - (x + 1)**3)}
    got = {'A_eq': A, 'B_eq': B, 'C_eq': C, 'Q_eq': Q}
    for nm in claims:
        out[nm] = dict(form=str(sp.factor(sp.expand(got[nm].subs(eq)))),
                       equals_claim=sp.expand(got[nm].subs(eq) - claims[nm]) == 0, claim=str(claims[nm]))
    par = {ys_: -ys_, as_: -bs_, bs_: -as_}
    out['parity_even_under_y->-y,a->-b,b->-a'] = {nm: sp.expand(e.subs(par, simultaneous=True) - e) == 0
                                                  for nm, e in dict(A=A, C=C, Q=Q).items()}
    S1 = xs_ + 1 + (xs_ - 1)*as_*bs_; S2 = (1 + ys_)*as_ + (1 - ys_)*bs_
    out['B_is_S1^2+S2^2'] = sp.expand(B - S1**2 - S2**2) == 0
    out['S1_even_S2_odd'] = [sp.expand(S1.subs(par, simultaneous=True) - S1) == 0,
                             sp.expand(S2.subs(par, simultaneous=True) + S2) == 0]
    # N on the equator: E = e^{-psi} > 0, ring a = +sqrt((x+1)/(x-1)) (a > 0 for chi in (0,1) since alpha < 0)
    k, w0, E = sp.symbols('k w0 E', positive=True)
    Nexpr = E**2*k**2*Q + 4*k*w0*C - w0**2*A/E**2
    N_over_F = sp.cancel(sp.expand(Nexpr.subs(eq))/F)
    ring = {a: sp.sqrt((x + 1)/(x - 1))}
    Rr = sp.sqrt(x**2 - 1)
    omr = k*E**2*Rr - w0                              # omega at the ring (claim)
    chk = {}
    chk['N/F_at_ring == -(4x/(x-1)) (k R E - w0/E)^2'] = sp.simplify(
        (N_over_F.subs(ring) + 4*x/(x - 1)*(k*Rr*E - w0/E)**2).subs(x, sp.Symbol('X', positive=True) + 1)) == 0
    AF = sp.cancel(sp.expand(A.subs(eq))/F); CF = sp.cancel(sp.expand(C.subs(eq))/F); QF = sp.cancel(sp.expand(Q.subs(eq))/F)
    chk['A/F_at_ring == 4x/(x-1)'] = sp.simplify(AF.subs(ring) - 4*x/(x - 1)) == 0
    chk['Q/F_at_ring == -4x(x+1)'] = sp.simplify(QF.subs(ring) + 4*x*(x + 1)) == 0
    om_ring = 2*k*E**2*CF/AF - w0
    chk['omega_ring == k e^{-2psi} sqrt(x^2-1) - w0'] = sp.simplify((om_ring.subs(ring) - omr).subs(x, sp.Symbol('X', positive=True) + 1)) == 0
    # g/f at ring = (N/F) E^2 / (A/F) = -omega_ring^2
    chk['g_phiphi/f -> -omega_ring^2 (equator)'] = sp.simplify(((N_over_F*E**2/AF).subs(ring) + omr**2).subs(x, sp.Symbol('X', positive=True) + 1)) == 0
    out['equator_N_checks'] = chk
    # transversality: on y = 0, a^2 = alpha^2 exp(4 beta h(x)), h = 1 - x/R + x/(2R^3), R = sqrt(x^2-1)
    X = sp.Symbol('X', positive=True)
    Rx = sp.sqrt(X**2 - 1); h = 1 - X/Rx + X/(2*Rx**3)
    out['h_prime == -3/(2 R^5)'] = sp.simplify(sp.diff(h, X) + R_(3, 2)/Rx**5) == 0
    x_, y_ = xs_, ys_
    Rg = sp.sqrt(x_**2 + y_**2 - 1); ug = x_*y_/Rg; Pg = [1, ug, (3*ug**2 - 1)/2]
    g_a = -1 + sum((x_ - y_)*Pg[l]/Rg**(l + 1) for l in range(3))
    out['d_y g_a(x,0) == 3/(2 R^3)  (a = -alpha e^{-2 beta g_a})'] = sp.simplify((sp.diff(g_a, y_).subs(y_, 0) - R_(3, 2)/(x_**2 - 1)**R_(3, 2)).subs(x_, X + 1)) == 0
    out['consequence'] = ('at the ring (lambda = 1 - 3 beta/R^3): dS1/dx = -dF/dx = -2 lambda/(x-1), dS2/dy = 2 a lambda; '
                          'B ~ 4 lambda^2 [dx^2/(x-1)^2 + a^2 y^2] = K_B (drho^2 + dz^2), K_B = 4 lambda^2 (x+1)/((x-1) k^2 x^2): isotropic in Weyl coordinates')
    out['dF/dx_at_ring'] = '4 beta h\'(x)(x+1) + 2/(x-1) = 2/(x-1) - 6 beta (x+1)/R^5 ; zero iff R^3 = 3 beta (double root)'
    # y = 0 reduction of a, b (a = -b) -- symbolic check of the MN functions
    al, be = sp.symbols('alpha beta', real=True)
    a_f, b_f, psi_f = mn_ab_sym(al, be)
    out['a(x,0) == -b(x,0)'] = sp.simplify((a_f + b_f).subs(ys_, 0).subs(xs_, X)) == 0
    return out

def mn_ab_sym(al, be):
    x, y = xs_, ys_
    R = sp.sqrt(x**2 + y**2 - 1); u = x*y/R
    P = [sp.Integer(1), u, (3*u**2 - 1)/2]
    psi = be*P[2]/R**3
    a = -al*sp.exp(-2*be*(-1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))))
    b = al*sp.exp(2*be*(1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))))
    return a, b, psi

# ------------------------------------------------------------------------------------------------ mpmath field
class MN:
    def __init__(self, M, asp, beta):
        self.M, self.asp, self.beta = mp.mpf(M), mp.mpf(asp), mp.mpf(beta)
        chi = self.asp/self.M
        self.k = mp.sqrt(self.M**2 - self.asp**2)
        self.al = (-1 + mp.sqrt(1 - chi**2))/chi
        self.w0 = 4*self.k*self.al/(1 - self.al**2)
        self.P = poly_funcs()
    def ab(self, x, y):
        be = self.beta
        if be == 0:
            return -self.al, self.al, mp.mpf(0)
        R = mp.sqrt(x*x + y*y - 1); u = x*y/R
        P = [1, u, (3*u*u - 1)/2]
        psi = be*P[2]/R**3
        a = -self.al*mp.exp(-2*be*(-1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))))
        b = self.al*mp.exp(2*be*(1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))))
        return a, b, psi
    def fields(self, x, y):
        a, b, psi = self.ab(x, y); P = self.P; k, w0 = self.k, self.w0
        A, B, C, Q = P['A'](x, y, a, b), P['B'](x, y, a, b), P['C'](x, y, a, b), P['Q'](x, y, a, b)
        e2 = mp.exp(2*psi)
        N = k**2*Q/e2 + 4*k*w0*C - w0**2*e2*A
        rho2 = k**2*(x*x - 1)*(1 - y*y)
        d = dict(a=a, b=b, psi=psi, A=A, B=B, C=C, Q=Q, N=N, rho2=rho2, S1=P['S1'](x, y, a, b), S2=P['S2'](x, y, a, b))
        d['gpp'] = N/B
        if A != 0:
            f = e2*A/B; om = 2*k*C/(e2*A) - w0
            d.update(f=f, om=om, gpp_unreg=-f*om**2 + rho2/f)
        return d
    def xy_from_weyl(self, rho, z):
        k = self.k
        rp = mp.sqrt(rho**2 + (z + k)**2); rm = mp.sqrt(rho**2 + (z - k)**2)
        return (rp + rm)/(2*k), (rp - rm)/(2*k)
    def weyl(self, x, y):
        return self.k*mp.sqrt(x*x - 1)*mp.sqrt(1 - y*y), self.k*x*y

def sym_fields(M, asp, beta):
    """Exact-parameter sympy expressions (second implementation for S2)."""
    M, asp, beta = map(sp.nsimplify, (M, asp, beta))
    chi = asp/M; al = (-1 + sp.sqrt(1 - chi**2))/chi; k = sp.sqrt(M**2 - asp**2); w0 = 4*k*al/(1 - al**2)
    a, b, psi = mn_ab_sym(al, beta)
    A, B, C, Q = polys()
    sub = {as_: a, bs_: b}
    Ax, Bx, Cx, Qx = [e.subs(sub) for e in (A, B, C, Q)]
    N = sp.exp(-2*psi)*k**2*Qx + 4*k*w0*Cx - w0**2*sp.exp(2*psi)*Ax
    return dict(A=Ax, B=Bx, C=Cx, Q=Qx, N=N)

# ------------------------------------------------------------------------------------------------ points
POINTS = {
    'C0': dict(M=1, a='0.9', beta=(R_(19, 20)/(R_(19, 100))**R_(3, 2)), xs_given='1.92982825852', exact=(1, R_(9, 10), R_(19, 20)/(R_(19, 100))**R_(3, 2))),
    'F5R': dict(M=89, a=80, beta=R_(13394411, 1186380), xs_given='1.92192258792', exact=(89, 80, R_(13394411, 1186380))),
}

def mk(pt):
    p = POINTS[pt]
    return MN(p['M'], mp.mpf(p['a']) if isinstance(p['a'], str) else p['a'], mp.mpf(str(sp.N(p['beta'], mp.mp.dps + 10))))

def refine_ring(F, guess):
    g = lambda X: F.fields(X, mp.mpf(0))['S1']
    xs = mp.findroot(g, mp.mpf(guess))
    return xs

# ------------------------------------------------------------------------------------------------ S2
def taylor(F, xs, names, order=2):
    out = {}
    for nm in names:
        fn = lambda X, Y, nm=nm: F.fields(X, Y)[nm]
        c = {}
        for i in range(order + 1):
            for j in range(order + 1 - i):
                c[f'{i}{j}'] = mp.diff(fn, (xs, mp.mpf(0)), (i, j))
        out[nm] = c
    return out

def section_S2(pt):
    mp.mp.dps = 50
    F = mk(pt)
    xs = refine_ring(F, POINTS[pt]['xs_given'])
    d0 = F.fields(xs, mp.mpf(0))
    names = ['A', 'B', 'C', 'Q', 'N', 'S1', 'S2']
    T = taylor(F, xs, names)
    # second implementation: sympy derivatives of exact-parameter expressions, evaluated by mpmath lambdify
    t0 = time.time()
    S = sym_fields(*POINTS[pt]['exact'])
    T2 = {}
    for nm in ['A', 'B', 'C', 'N']:
        e = S[nm]; c = {}
        for (i, j) in [(1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]:
            de = sp.diff(e, xs_, i, ys_, j) if (i and j) else (sp.diff(e, xs_, i) if i else sp.diff(e, ys_, j))
            c[f'{i}{j}'] = sp.lambdify((xs_, ys_), de, 'mpmath')(xs, mp.mpf(0))
        T2[nm] = c
    t_sym = time.time() - t0
    agree = {}
    for nm in T2:
        for key, v in T2[nm].items():
            u = T[nm][key]; scale = max(abs(v), abs(u), mp.mpf('1e-40'))
            agree[f'{nm}_{key}'] = float(abs(u - v)/scale) if max(abs(u), abs(v)) > mp.mpf('1e-30') else float(abs(u - v))
    # closed-form predictions from S1
    x = xs; R = mp.sqrt(x*x - 1); psi = d0['psi']; k, w0, be = F.k, F.w0, F.beta
    Fp = 2/(x - 1) - 6*be*(x + 1)/R**5                       # dF/dx along the equator at the ring
    om_r = k*mp.exp(-2*psi)*d0['a']*(x - 1) - w0          # = k e^{-2psi} sqrt(x^2-1) - w0 since a_ring = +sqrt((x+1)/(x-1))
    assert abs(d0['a']*(x - 1) - R) < mp.mpf(10)**(-30)
    S2y = T['S2']['01']
    pred = {'A_x': 4*x/(x - 1)*Fp, 'B_xx/2': Fp**2, 'B_yy/2': S2y**2, 'C_x': 2*d0['a']*x*Fp,
            'Q_x': -4*x*(x + 1)*Fp, 'N_x': -4*x/(x - 1)*mp.exp(2*psi)*om_r**2*Fp}
    num = {'A_x': T['A']['10'], 'B_xx/2': T['B']['20']/2, 'B_yy/2': T['B']['02']/2, 'C_x': T['C']['10'],
           'Q_x': T['Q']['10'], 'N_x': T['N']['10']}
    lam = 1 - 3*be/R**3
    pred['B_yy/2'] = 4*d0['a']**2*lam**2                       # from d_y g_a = 3/(2R^3): dS2/dy = 2a(1 - 3beta/R^3)
    pred['dF/dx'] = 2/(x - 1)*lam; num['dF/dx'] = Fp
    predchk = {kk: float(abs(num[kk] - pred[kk])/abs(pred[kk])) for kk in pred}
    K_B = 4*lam**2*(x + 1)/((x - 1)*k**2*x**2)                    # B ~ K_B d^2,  d = Weyl distance to the ring
    K_f = mp.exp(2*psi)*2*k*x**2/(R*lam)                          # f ~ K_f cos(theta)/d, theta from the outward equator
    K_g = -om_r**2*K_f                                             # g_phiphi ~ K_g cos(theta)/d
    iso = dict(lambda_=mp.nstr(lam, 15), K_B=mp.nstr(K_B, 15), K_f=mp.nstr(K_f, 15), K_g=mp.nstr(K_g, 15))
    cA = -T['A']['02']/(2*T['A']['10']); cN = -T['N']['02']/(2*T['N']['10'])
    # curvature in Weyl: delta rho = (k x/R) dx - (k R/2) y^2, delta z = k x y  ->  delta rho = kappa delta z^2
    kap = lambda c: ((k*x/R)*c - k*R/2)/(k*x)**2
    res = dict(point=pt, M=str(F.M), a_spin=str(F.asp), beta=mp.nstr(be, 20), k=mp.nstr(k, 20), alpha=mp.nstr(F.al, 20),
               w0=mp.nstr(w0, 20), x_s=mp.nstr(xs, 30), x_s_given=POINTS[pt]['xs_given'],
               x_s_minus_given=float(xs - mp.mpf(POINTS[pt]['xs_given'])),
               rho_s=mp.nstr(F.weyl(xs, mp.mpf(0))[0], 20), R3_at_ring=mp.nstr(R**3, 15), three_beta=mp.nstr(3*be, 15),
               values_at_ring={nm: mp.nstr(d0[nm], 6) for nm in ['A', 'B', 'C', 'Q', 'N', 'S1', 'S2']},
               a_at_ring=mp.nstr(d0['a'], 20), b_at_ring=mp.nstr(d0['b'], 20), psi_at_ring=mp.nstr(psi, 20),
               taylor={nm: {kk: mp.nstr(v, 15) for kk, v in T[nm].items()} for nm in names},
               two_impl_max_rel_diff=max(agree.values()), two_impl_detail=agree, sympy_diff_seconds=round(t_sym, 1),
               closed_form_prediction={kk: mp.nstr(v, 15) for kk, v in pred.items()}, closed_form_rel_err=predchk,
               dF_dx=mp.nstr(Fp, 15), omega_ring=mp.nstr(om_r, 20), isotropic_leading_forms=iso,
               B_quadratic_form=dict(Bxx_half=mp.nstr(T['B']['20']/2, 15), Bxy=mp.nstr(T['B']['11'], 6),
                                     Byy_half=mp.nstr(T['B']['02']/2, 15),
                                     positive_definite=bool(T['B']['20'] > 0 and T['B']['02'] > 0 and T['B']['20']*T['B']['02'] > T['B']['11']**2)),
               A_leading=dict(A_x=mp.nstr(T['A']['10'], 15), A_y=mp.nstr(T['A']['01'], 6), A_yy_half=mp.nstr(T['A']['02']/2, 15)),
               N_leading=dict(N0=mp.nstr(d0['N'], 6), N_x=mp.nstr(T['N']['10'], 15), N_y=mp.nstr(T['N']['01'], 6),
                              N_yy_half=mp.nstr(T['N']['02']/2, 15)),
               curve_A0_dx_over_y2=mp.nstr(cA, 15), curve_N0_dx_over_y2=mp.nstr(cN, 15),
               curve_A0_drho_over_dz2=mp.nstr(kap(cA), 15), curve_N0_drho_over_dz2=mp.nstr(kap(cN), 15),
               ctc_side_on_equator='x < x_s (inner)' if T['N']['10'] > 0 else 'x > x_s (outer)',
               ergo_side_on_equator='x > x_s (outer)' if T['A']['10'] < 0 else 'x < x_s (inner)')
    return res, F, xs, T

# ------------------------------------------------------------------------------------------------ S3 rays
def section_S3(pt, F, xs, om_r, iso=None):
    mp.mp.dps = 50
    rs, _ = F.weyl(xs, mp.mpf(0))
    thetas = [0, 45, 80, 89, 90, 91, 100, 135, 180]
    ds = [mp.mpf(10)**(-e) for e in range(2, 7)]
    out = []
    for th in thetas:
        t = mp.radians(th); row = dict(theta_deg=th, d=[], gpp=[], f=[], omega=[], reg_vs_unreg_rel=[], gpp_over_minus_om2f=[])
        for d in ds:
            X, Y = F.xy_from_weyl(rs + d*mp.cos(t), d*mp.sin(t))
            q = F.fields(X, Y)
            row['d'].append(float(d)); row['gpp'].append(mp.nstr(q['gpp'], 12))
            row['f'].append(mp.nstr(q.get('f', mp.nan), 12)); row['omega'].append(mp.nstr(q.get('om', mp.nan), 12))
            row['reg_vs_unreg_rel'].append(float(abs(q['gpp'] - q['gpp_unreg'])/abs(q['gpp'])))
            row['gpp_over_minus_om2f'].append(mp.nstr(q['gpp']/(-om_r**2*q['f']), 12))
        lg = [mp.log10(abs(mp.mpf(g))) for g in row['gpp']]; ld = [mp.log10(d) for d in ds]
        row['slope_log_gpp_vs_log_d(last3)'] = float((lg[-1] - lg[-3])/(ld[-1] - ld[-3]))
        row['sign_gpp'] = int(mp.sign(mp.mpf(row['gpp'][-1])))
        if iso is not None and abs(mp.cos(t)) > mp.mpf('0.1'):
            d = ds[-1]; X, Y = F.xy_from_weyl(rs + d*mp.cos(t), d*mp.sin(t)); q = F.fields(X, Y)
            row['at_d=1e-6: f*d/(K_f cos) , gpp*d/(K_g cos), B/(K_B d^2)'] = [
                mp.nstr(q['f']*d/(mp.mpf(iso['K_f'])*mp.cos(t)), 10), mp.nstr(q['gpp']*d/(mp.mpf(iso['K_g'])*mp.cos(t)), 10),
                mp.nstr(q['B']/(mp.mpf(iso['K_B'])*d**2), 10)]
        out.append(row)
    return out

# ------------------------------------------------------------------------------------------------ S4 curves
def section_S4(F, xs, cA, cN):
    mp.mp.dps = 50
    rows = []
    for yv in ['1e-2', '3e-3', '1e-3', '3e-4']:
        Y = mp.mpf(yv); r = {'y': yv}
        for nm, c in [('A', cA), ('N', cN)]:
            g = lambda dx: F.fields(xs + dx, Y)[nm]
            dx = mp.findroot(g, c*Y**2)
            r[f'{nm}0_dx'] = mp.nstr(dx, 12); r[f'{nm}0_dx_over_y2'] = mp.nstr(dx/Y**2, 12)
        rows.append(r)
    return rows

# ------------------------------------------------------------------------------------------------ S5 global C0 map
def _worker(args):
    pt, rhos, z, dps = args
    mp.mp.dps = dps
    F = mk(pt)
    out = []
    for rho in rhos:
        X, Y = F.xy_from_weyl(mp.mpf(rho), mp.mpf(z))
        try:
            q = F.fields(X, Y)
            out.append((int(mp.sign(q['N'])), int(mp.sign(q['A']))))
        except (ZeroDivisionError, ValueError, OverflowError):
            out.append((0, 0))
    return out

def sign_grid(pt, rho, z, dps=30, procs=2):
    with get_context('spawn').Pool(procs) as P:
        rows = P.map(_worker, [(pt, list(map(str, rho)), str(zz), dps) for zz in z])
    arr = np.array(rows)       # (nz, nrho, 2)
    return arr[..., 0], arr[..., 1]

def write_png(path, rgb):
    h, w, _ = rgb.shape
    raw = b''.join(b'\x00' + rgb[i].astype(np.uint8).tobytes() for i in range(h))
    def chunk(t, d):
        c = struct.pack('>I', len(d)) + t + d
        return c + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    png = b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 2, 0, 0, 0)) + \
          chunk(b'IDAT', zlib.compress(raw, 9)) + chunk(b'IEND', b'')
    open(path, 'wb').write(png)
    return png

def colour(sN, sA):
    # classes: f>0 & g>0 white; f<0 & g>0 (ergo only) blue; f>0 & g<0 (CTC only) red; f<0 & g<0 purple; undefined grey
    img = np.full(sN.shape + (3,), 255, np.uint8)
    img[(sA < 0) & (sN > 0)] = (150, 190, 240)
    img[(sA > 0) & (sN < 0)] = (240, 150, 150)
    img[(sA < 0) & (sN < 0)] = (180, 120, 200)
    img[(sA == 0) | (sN == 0)] = (160, 160, 160)
    # zero curves: sign changes between neighbours
    for s, col in [(sA, (20, 60, 200)), (sN, (200, 20, 20))]:
        e = np.zeros(s.shape, bool)
        e[:, 1:] |= (s[:, 1:]*s[:, :-1] < 0); e[1:, :] |= (s[1:, :]*s[:-1, :] < 0)
        img[e] = col
    return img

def section_S5(F, xs, pt='C0'):
    k = float(F.k); rs = float(F.weyl(xs, mp.mpf(0))[0])
    # global grid (z >= 0, mirrored by the equatorial symmetry N, A even in y)
    h = 0.005
    rho = np.arange(h/2, 1.6, h); z = np.arange(0.0, 1.2 + 1e-12, h)
    t0 = time.time()
    sN, sA = sign_grid(pt, rho, z)
    tg = time.time() - t0
    sN_full = np.vstack([sN[:0:-1], sN]); sA_full = np.vstack([sA[:0:-1], sA]); z_full = np.concatenate([-z[:0:-1], z])
    ctc = sN_full < 0
    lab, nlab = ndimage.label(ctc)
    comps = []
    for L in range(1, nlab + 1):
        m = lab == L; iz, ir = np.nonzero(m)
        comps.append(dict(label=L, n_pixels=int(m.sum()), rho_min=float(rho[ir].min() - h/2), rho_max=float(rho[ir].max() + h/2),
                          z_min=float(z_full[iz].min() - h/2), z_max=float(z_full[iz].max() + h/2),
                          touches_rho0_column=bool((ir == 0).any()),
                          rho0_z_range=[float(z_full[iz[ir == 0]].min()), float(z_full[iz[ir == 0]].max())] if (ir == 0).any() else None,
                          contains_ring_pixel_nbhd=bool(m[np.abs(z_full) < 2*h][:, np.abs(rho - rs) < 2*h].any())))
    ergo = sA_full < 0; labE, nE = ndimage.label(ergo)
    ergo_comps = []
    for L in range(1, nE + 1):
        m = labE == L; iz, ir = np.nonzero(m)
        ergo_comps.append(dict(n_pixels=int(m.sum()), rho_min=float(rho[ir].min() - h/2), rho_max=float(rho[ir].max() + h/2),
                               z_min=float(z_full[iz].min() - h/2), z_max=float(z_full[iz].max() + h/2),
                               touches_rho0_column=bool((ir == 0).any())))
    overlap_pixels = int(((sN_full < 0) & (sA_full < 0)).sum())
    undefined = int(((sN_full == 0) | (sA_full == 0)).sum())
    # refine extents of each CTC component (z >= 0 representatives; the set is symmetric under z -> -z)
    mp.mp.dps = 40
    def Nw(r, zz):
        X, Y = F.xy_from_weyl(mp.mpf(r), mp.mpf(zz)); return F.fields(X, Y)['N']
    def edge(fixed, along, lo, hi, which):
        """root of N on a line: which='rho' -> N(rho, fixed) for rho in [lo, hi]; 'z' -> N(fixed, z). N(lo) < 0 < N(hi)."""
        g = (lambda t: Nw(t, fixed)) if which == 'rho' else (lambda t: Nw(fixed, t))
        a_, b_ = mp.mpf(lo), mp.mpf(hi)
        if not (g(a_) < 0 < g(b_)): return None
        for _ in range(70):
            m = (a_ + b_)/2
            if g(m) < 0: a_ = m
            else: b_ = m
        return (a_ + b_)/2
    def refine_max(mask, which):
        iz_, ir_ = np.nonzero(mask)
        if which == 'rho':            # outer rho edge per row, maximise over z
            cand = sorted(((rho[ir_[iz_ == i]].max(), z_full[i]) for i in np.unique(iz_)), reverse=True)[:1]
            r0, z0 = cand[0]
            fn = lambda zz: edge(zz, 'rho', r0 - 6*h, r0 + 3*h, 'rho') or mp.mpf(0)
            lo, hi = z0 - 2*h, z0 + 2*h
        else:                         # top z edge per column, maximise over rho
            cand = sorted(((z_full[iz_[ir_ == j]].max(), rho[j]) for j in np.unique(ir_)), reverse=True)[:1]
            z0, r0 = cand[0]
            fn = lambda rr: edge(rr, 'z', z0 - 6*h, z0 + 3*h, 'z') or mp.mpf(0)
            lo, hi = r0 - 2*h, r0 + 2*h
        lo, hi = lo - 8*h, hi + 8*h                     # scan a wide window first, then golden-section locally
        ts = [lo + (hi - lo)*i/48 for i in range(49)]; vs = [fn(t) for t in ts]; ib = max(range(49), key=lambda i: vs[i])
        lo, hi = ts[max(ib - 1, 0)], ts[min(ib + 1, 48)]
        gr = (mp.sqrt(5) - 1)/2
        for _ in range(45):
            m1 = hi - gr*(hi - lo); m2 = lo + gr*(hi - lo)
            if fn(m1) > fn(m2): hi = m2
            else: lo = m1
        arg = (lo + hi)/2; val = fn(arg)
        return dict(max=float(val), at=float(arg), window_scan_best=float(vs[ib]))
    refined = []
    for c in comps:
        m = lab == c['label']
        if c['z_min'] < 0 and c['z_max'] > 0:      # equatorial lobe: symmetric; use its z >= 0 half
            mh = m & (z_full[:, None] >= 0)
            r_eq = dict(label=c['label'], kind='equatorial lobe (contains the ring in its closure)',
                        rho_sup='rho_s (ring tip, not attained): ' + mp.nstr(F.weyl(xs, mp.mpf(0))[0], 15),
                        rho_edge_at_z=[dict(z=zz, rho=float(edge(zz, 'rho', rs - 0.01, rs + 0.002, 'rho')))
                                       for zz in [0.002, 0.005, 0.01, 0.02]],
                        z_max=refine_max(mh, 'z'))
        elif c['z_min'] >= 0:
            r_eq = dict(label=c['label'], kind='off-equator lobe z > 0 (mirror lobe at z < 0)',
                        rho_max=refine_max(m, 'rho'), z_max=refine_max(m, 'z'))
        else:
            continue
        refined.append(r_eq)
    # equatorial profile: sign changes of N (g_phiphi) and A (f) on y = 0 between x = 1 and x = 8
    def roots_eq(nm):
        xsamp = [1 + mp.mpf(10)**(-6 + 6*i/400.0) for i in range(401)] + [2 + 6*mp.mpf(i)/200 for i in range(1, 201)]
        vals = [F.fields(X, mp.mpf(0))[nm] for X in xsamp]; out_ = []
        for i in range(len(xsamp) - 1):
            if vals[i]*vals[i + 1] < 0:
                a_, b_ = xsamp[i], xsamp[i + 1]; ga = vals[i]
                for _ in range(80):
                    m = (a_ + b_)/2; gm = F.fields(m, mp.mpf(0))[nm]
                    if gm*ga > 0: a_, ga = m, gm
                    else: b_ = m
                xr = (a_ + b_)/2
                out_.append(dict(x=mp.nstr(xr, 15), rho=float(F.weyl(xr, mp.mpf(0))[0]),
                                 sign_left=int(mp.sign(vals[i])), sign_right=int(mp.sign(vals[i + 1]))))
        return out_
    equator = dict(N_sign_changes=roots_eq('N'), A_sign_changes=roots_eq('A'),
                   note='B > 0 off the ring, so sign(g_phiphi) = sign(N) and sign(f) = sign(A); at the ring both change sign (simple factor F)')
    # horizon rod x = 1, 0 < y < 1: where g_phiphi changes sign
    # (at x = 1 exactly, rho^2 = 0 and N suffers cancellation between terms ~ e^{+-2psi}; use x = 1 + eps, dps 80)
    rod_roots = []
    mp.mp.dps = 80
    for eps in ['1e-8', '1e-12', '1e-16']:
        X1 = 1 + mp.mpf(eps); gy = lambda Y: F.fields(X1, Y)['N']
        ys = [mp.mpf(i)/100 for i in range(1, 100)]; vy = [gy(Y) for Y in ys]
        for i in range(len(ys) - 1):
            if vy[i]*vy[i + 1] < 0:
                a_, b_ = ys[i], ys[i + 1]
                for _ in range(60):
                    m = (a_ + b_)/2
                    if gy(m)*vy[i] > 0: a_ = m
                    else: b_ = m
                rod_roots.append(dict(x_minus_1=eps, y=mp.nstr((a_ + b_)/2, 12), z=float(F.k*(a_ + b_)/2),
                                      sign_below=int(mp.sign(vy[i])), sign_above=int(mp.sign(vy[i + 1]))))
    # full re-check of the near-rod pixels (rho < 0.05) at dps 80
    near_bad = 0; near_n = 0
    for j in np.nonzero(rho < 0.05)[0]:
        for i in range(len(z)):
            X, Y = F.xy_from_weyl(mp.mpf(str(rho[j])), mp.mpf(str(z[i])))
            d = F.fields(X, Y); near_n += 1
            if int(mp.sign(d['N'])) != sN[i, j] or int(mp.sign(d['A'])) != sA[i, j]: near_bad += 1
    mp.mp.dps = 40
    # polar scan about the singular horizon-equator point (rho, z) = (0, 0), i.e. x = 1, y = 0: angular sign pattern
    origin = []
    for rW in ['0.05', '0.02', '0.01', '0.005']:
        pat_N = ''; pat_A = ''; ok = True
        for th in np.linspace(0.5, 89.5, 90):          # angle from the equator towards +z, degrees
            t = mp.radians(th); r_, z_ = mp.mpf(rW)*mp.cos(t), mp.mpf(rW)*mp.sin(t)
            vals = []
            for dps in (40, 80):
                mp.mp.dps = dps
                X, Y = F.xy_from_weyl(r_, z_); d = F.fields(X, Y); vals.append((int(mp.sign(d['N'])), int(mp.sign(d['A']))))
            ok &= vals[0] == vals[1]
            pat_N += '-' if vals[1][0] < 0 else '+'; pat_A += '-' if vals[1][1] < 0 else '+'
        mp.mp.dps = 40
        origin.append(dict(weyl_r=rW, angles_deg='0.5..89.5 step 1', sign_gpp=pat_N, sign_f=pat_A, dps40_eq_dps80=bool(ok)))
    # boundedness checks
    mp.mp.dps = 40
    far = []
    for rW in [2, 5, 20, 100]:
        signs = set()
        for th in np.linspace(0.01, np.pi - 0.01, 41):
            signs.add(int(mp.sign(Nw(rW*np.sin(th), rW*np.cos(th)))))
        far.append(dict(weyl_r=rW, signs_of_gpp=sorted(signs)))
    axis = []
    for zz in [0.45, 0.5, 0.7, 1.0, 2.0, 5.0]:
        q = Nw(1e-8, zz); X, Y = F.xy_from_weyl(mp.mpf('1e-8'), mp.mpf(zz)); d = F.fields(X, Y)
        axis.append(dict(z=zz, gpp_over_rho2=mp.nstr(d['gpp']/d['rho2'], 8)))
    rod = []   # horizon rod x -> 1, |y| < 1 (rho -> 0, |z| < k)
    for yv in [0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
        d = F.fields(1 + mp.mpf('1e-12'), mp.mpf(yv))
        rod.append(dict(y=yv, z=float(F.k*yv), gpp=mp.nstr(d['gpp'], 8), f=mp.nstr(d.get('f', mp.nan), 8)))
    # dps cross-check of the grid signs at random pixels
    rng = np.random.default_rng(15); bad = 0; nchk = 300
    for _ in range(nchk):
        i = rng.integers(len(z)); j = rng.integers(len(rho))
        mp.mp.dps = 70
        X, Y = F.xy_from_weyl(mp.mpf(str(rho[j])), mp.mpf(str(z[i])))
        try:
            d = F.fields(X, Y)
            if int(mp.sign(d['N'])) != sN[i, j] or int(mp.sign(d['A'])) != sA[i, j]: bad += 1
        except Exception:
            pass
    mp.mp.dps = 40
    img = colour(sN_full[::-1], sA_full[::-1])     # z up
    # ring marker and axis
    iz0 = np.argmin(np.abs(z_full[::-1])); ir0 = np.argmin(np.abs(rho - rs))
    img[max(iz0 - 2, 0):iz0 + 3, max(ir0 - 2, 0):ir0 + 3] = (0, 0, 0)
    kk = np.abs(z_full[::-1]) < k; img[kk, 0:2] = (0, 140, 0)          # horizon rod (green)
    png = write_png(os.path.join(RES, 'v15_q3_c0.png'), img)
    # zoom near the ring
    hz = 2e-4; rz = rs + np.arange(-0.02 + hz/2, 0.02, hz); zz_ = np.arange(0.0, 0.02 + 1e-12, hz)
    sNz, sAz = sign_grid(pt, rz, zz_, dps=40)
    sNz = np.vstack([sNz[:0:-1], sNz]); sAz = np.vstack([sAz[:0:-1], sAz])
    imz = colour(sNz[::-1], sAz[::-1]); c = imz.shape[0]//2; cc = imz.shape[1]//2
    imz[c - 1:c + 2, cc - 1:cc + 2] = (0, 0, 0)
    pngz = write_png(os.path.join(RES, 'v15_q3_c0_zoom.png'), imz)
    write_svg(png, img.shape, pngz, imz.shape, rho, z_full, rz, k, rs)
    return dict(grid=dict(h=h, rho_range=[float(rho[0]), float(rho[-1])], z_range=[-1.2, 1.2], pixels_full=list(sN_full.shape), seconds=round(tg, 1),
                          dps=30, sign_recheck_dps70_mismatches=f'{bad}/{nchk}', undefined_pixels=undefined),
                ctc_components=comps, ergo_components=ergo_comps, ctc_and_ergo_overlap_pixels=overlap_pixels,
                ctc_refined_extents=refined, origin_polar_scan=origin, equator_profile=equator, horizon_rod_gpp_sign_changes=rod_roots, near_rod_pixels_dps80_mismatches=f'{near_bad}/{near_n}',
                far_field_signs=far, axis_gpp_over_rho2=axis, horizon_rod=rod,
                zoom=dict(box_rho=[float(rz[0]), float(rz[-1])], box_z=[-0.02, 0.02], h=hz,
                          ctc_and_ergo_overlap_pixels=int(((sNz < 0) & (sAz < 0)).sum())))

def write_svg(png, shp, pngz, shpz, rho, z_full, rz, k, rs):
    H, W = shp[0], shp[1]; Hz, Wz = shpz[0], shpz[1]
    s = 2; sz = 2
    tw = 60 + W*s + 80 + Wz*sz + 40; th = 40 + max(H*s, Hz*sz) + 120
    b64 = base64.b64encode(png).decode(); b64z = base64.b64encode(pngz).decode()
    x2 = 60 + W*s + 80
    t = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{tw}" height="{th}" font-family="sans-serif" font-size="13">',
         '<rect width="100%" height="100%" fill="white"/>',
         f'<text x="60" y="24">MN C0 (chi 0.9, q 0.95): Weyl (rho, z), 0 &lt; rho &lt; 1.6, |z| &lt; 1.2</text>',
         f'<image x="60" y="40" width="{W*s}" height="{H*s}" style="image-rendering:pixelated" href="data:image/png;base64,{b64}"/>',
         f'<rect x="60" y="40" width="{W*s}" height="{H*s}" fill="none" stroke="black"/>']
    for r in [0, 0.4, 0.8, 1.2, 1.6]:
        X = 60 + r/1.6*W*s; t.append(f'<text x="{X-8:.0f}" y="{40+H*s+16}">{r}</text>')
    for zz in [-1.2, -0.6, 0, 0.6, 1.2]:
        Y = 40 + (1.2 - zz)/2.4*H*s; t.append(f'<text x="22" y="{Y+4:.0f}">{zz}</text>')
    t.append(f'<text x="{60+W*s/2-15}" y="{40+H*s+34}">rho</text><text x="4" y="{40+H*s/2}">z</text>')
    t += [f'<text x="{x2}" y="24">zoom: |rho - rho_s| &lt; 0.02, |z| &lt; 0.02 (ring at centre)</text>',
          f'<image x="{x2}" y="40" width="{Wz*sz}" height="{Hz*sz}" style="image-rendering:pixelated" href="data:image/png;base64,{b64z}"/>',
          f'<rect x="{x2}" y="40" width="{Wz*sz}" height="{Hz*sz}" fill="none" stroke="black"/>',
          f'<text x="{x2-10}" y="{40+Hz*sz+16}">{rs-0.02:.3f}</text><text x="{x2+Wz*sz-30}" y="{40+Hz*sz+16}">{rs+0.02:.3f}</text>']
    ly = 40 + max(H*s, Hz*sz) + 50
    leg = [((255, 255, 255), 'f &gt; 0, g_phiphi &gt; 0'), ((150, 190, 240), 'ergoregion only (f &lt; 0, g_phiphi &gt; 0)'),
           ((240, 150, 150), 'CTC only (f &gt; 0, g_phiphi &lt; 0)'), ((180, 120, 200), 'f &lt; 0 and g_phiphi &lt; 0'),
           ((20, 60, 200), 'A = 0 (ergosurface)'), ((200, 20, 20), 'g_phiphi = 0 (CTC boundary)'),
           ((0, 0, 0), 'ring B = A = 0'), ((0, 140, 0), 'horizon rod x = 1')]
    for i, (c, lbl) in enumerate(leg):
        X = 60 + (i % 4)*300; Y = ly + (i//4)*24
        t.append(f'<rect x="{X}" y="{Y-11}" width="14" height="14" fill="rgb{c}" stroke="black"/><text x="{X+20}" y="{Y}">{lbl}</text>')
    t.append('</svg>')
    open(os.path.join(RES, 'v15_q3_c0.svg'), 'w').write('\n'.join(t))

# ------------------------------------------------------------------------------------------------ S6 Kerr control
def section_S6():
    mp.mp.dps = 50
    M, asp = mp.mpf(1), mp.mpf('0.9')
    F = MN(M, asp, 0); k = F.k
    xs = -M/k
    out = dict(ring_x=mp.nstr(xs, 20), S1_at_ring=mp.nstr(F.fields(xs, mp.mpf(0))['S1'], 5))
    rows = []
    for (dx, yv) in [('1e-3', '0'), ('-1e-3', '0'), ('1e-3', '1e-3'), ('0', '1e-3'), ('2e-1', '3e-1')]:
        X = xs + mp.mpf(dx); Y = mp.mpf(yv)
        d = F.fields(X, Y)
        r = M + k*X; c2 = Y**2; s2 = 1 - c2; Sig = r*r + asp**2*c2
        g_bl = (r*r + asp**2 + 2*M*asp**2*r*s2/Sig)*s2
        f_bl = 1 - 2*M*r/Sig
        rows.append(dict(dx=dx, y=yv, gpp_MN_beta0=mp.nstr(d['gpp'], 15), gpp_BL=mp.nstr(g_bl, 15),
                         rel=float(abs(d['gpp'] - g_bl)/abs(g_bl)), f_rel=float(abs(d['f'] - f_bl)/abs(f_bl))))
    out['checks'] = rows
    # general ring formula omega_r = k e^{-2psi} a_ring (x_s - 1) - w0  (= k e^{-2psi} sqrt(x^2-1) - w0 when x_s > 1);
    # for Kerr x_s < -1 so a(x-1) = -sqrt(x^2-1) and omega_r = -a_spin + 2 a_spin = a_spin (BL value of g_tphi/f at Sigma -> 0)
    om_r = k*(-F.al)*(xs - 1) - F.w0
    out['omega_ring_formula'] = mp.nstr(om_r, 15); out['omega_ring_BL_expected'] = mp.nstr(asp, 15)
    out['gpp_over_minus_om2f_near_ring'] = []
    for dx in ['1e-4', '-1e-4', '1e-6']:
        d = F.fields(xs + mp.mpf(dx), mp.mpf(dx)/3)
        out['gpp_over_minus_om2f_near_ring'].append(dict(dx=dx, y=f'{dx}/3', ratio=mp.nstr(d['gpp']/(-om_r**2*d['f']), 12)))
    out['note'] = 'Kerr ring r = 0 is x = -M/k < -1; there a = -alpha, b = alpha so the S1 equatorial factorisation applies with s = alpha^2'
    return out

# ------------------------------------------------------------------------------------------------ main
def versions():
    return dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__, numpy=np.__version__,
                scipy=scipy.__version__, platform=platform.platform(), matplotlib='not installed in venv (PNG written by own zlib writer)')

if __name__ == '__main__':
    T0 = time.time(); out = dict(versions=versions(), date='2026-10-17', prereg='V15 PREREGISTRATION.md Q3')
    print('S1 ...', flush=True); out['S1_symbolic'] = section_S1(); print(json.dumps(out['S1_symbolic'], indent=1, default=str), flush=True)
    out['S6_kerr_control'] = section_S6(); print('S6', json.dumps(out['S6_kerr_control'], indent=1), flush=True)
    keep = {}
    for pt in ['C0', 'F5R']:
        print(pt, 'S2 ...', flush=True)
        res, F, xs, T = section_S2(pt); keep[pt] = (F, xs, T)
        print(json.dumps({kk: v for kk, v in res.items() if kk not in ('taylor', 'two_impl_detail')}, indent=1), flush=True)
        om_r = mp.mpf(res['omega_ring'])
        res['S3_rays'] = section_S3(pt, F, xs, om_r, res['isotropic_leading_forms'])
        for r in res['S3_rays']:
            print(pt, 'theta', r['theta_deg'], 'slope', round(r['slope_log_gpp_vs_log_d(last3)'], 4), 'sign', r['sign_gpp'],
                  'g/(-om^2 f)', r['gpp_over_minus_om2f'][-1], 'reg-unreg', max(r['reg_vs_unreg_rel']),
                  r.get('at_d=1e-6: f*d/(K_f cos) , gpp*d/(K_g cos), B/(K_B d^2)'), flush=True)
        res['S4_curves'] = section_S4(F, xs, mp.mpf(res['curve_A0_dx_over_y2']), mp.mpf(res['curve_N0_dx_over_y2']))
        print(json.dumps(res['S4_curves'], indent=1), flush=True)
        out[pt] = res
    print('S5 ...', flush=True)
    mp.mp.dps = 40
    F, xs, T = keep['C0']
    out['S5_C0_global'] = section_S5(F, xs)
    print(json.dumps(out['S5_C0_global'], indent=1, default=str), flush=True)
    out['runtime_s'] = round(time.time() - T0, 1)
    json.dump(out, open(os.path.join(RES, 'v15_q3_geometry.json'), 'w'), indent=1, default=str)
    print('done', out['runtime_s'], 's')
