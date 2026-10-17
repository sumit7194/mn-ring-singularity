"""V15 Q2b: derive (not assume) the leading curvature invariant at a simple pole of the Ernst potential.

Part A  (symbolic, sympy): perturbative solution of the Ernst equation near a simple pole E ~ c/W, W = (rho-rho0) + i z,
         through degree 1, and the Laurent expansion of the exact identity
             I = (4/H^2) (Q_WW Q_WbWb + 3 Q_WWb^2),     H = e^{2 gamma}/f,
         giving the leading term  I -> 3 (c/cbar)^2 / (rho0^2 Lambda^2 W^6)   (free Laurent data drop out).
         Also: wider monomial ansatz (robustness), and the exact Kerr check of the coefficient.
Part B  (numeric, mpmath): the exact identity against (i) the Kerr closed forms at generic points,
         (ii) the Q2 direct-Riemann K values stored in results/v15_q2_curvature.json (MN, C0 and F5-R, near the ring),
         (iii) the Q2 negative-control point; plus the unsubstituted form (metric's own d ln H, no Ernst/gamma eq. used).
Part C  (numeric): Laurent data at the 4 rings (C0, F5-R, (chi,q)=(0.6,0.1024), Kerr's ring), closed-form
         K delta^6 and K y^6 limits, against Richardson-extrapolated numerical limits of K.
Metric transcription: V9_mn_obstruction_check/code/v9_nve.py :: mn_functions (imported, read-only) for f and e^{2gamma};
the Ernst potential E = e^{2psi} A_-/A_+ (Brink 2008 App. C form, as in Q2 Sec. B) is re-typed here with the same a, b, psi,
and Re E = f(v9) is checked numerically.
Run: nice -n 10 python v15_q2b_check.py   (single core; < 1 GB)
"""
import sys, os, json, time, platform, random
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.abspath(os.path.join(ROOT, '..', 'V9_mn_obstruction_check', 'code')))
from v9_nve import mn_functions  # noqa: E402  (verified transcription; read-only)

DPS = 90
mp.mp.dps = DPS
T0 = time.time()
OUT = dict(meta=dict(python=sys.version, sympy=sp.__version__, mpmath=mp.__version__, platform=platform.platform(),
                     dps=DPS, script='code/v15_q2b_check.py'))


def s(v, n=20):
    if isinstance(v, mp.mpc) or isinstance(v, complex):
        return [mp.nstr(mp.re(v), n), mp.nstr(mp.im(v), n)]
    return mp.nstr(v, n)


# =====================================================================================================================
# Part A: symbolic Laurent expansion
# =====================================================================================================================
_rng = random.Random(7)


def is_zero(v):
    v = sp.sympify(v)
    fs = sorted(v.free_symbols, key=str)
    for _ in range(2):
        sub = {q: sp.Rational(_rng.randint(1, 97), _rng.randint(1, 89)) + sp.I*sp.Rational(_rng.randint(1, 50), 71)
               for q in fs}
        val = complex(sp.N(v.subs(sub), 40))
        if abs(val) > 1e-25:
            return False
    return True


class S:
    """Truncated Laurent series in the scaling parameter eps (W -> eps W, Wb -> eps Wb).
    terms {p: coefficient(W, Wb)}; the series is exact through power `top`."""
    def __init__(s_, terms, top):
        s_.t = {p: v for p, v in terms.items() if p <= top}; s_.top = top

    def low(s_):
        nz = [p for p, v in s_.t.items() if not is_zero(v)]
        return min(nz) if nz else s_.top

    def __add__(a, b):
        if not isinstance(b, S): b = S({0: sp.sympify(b)}, 10**6)
        t = dict(a.t)
        for p, v in b.t.items(): t[p] = t.get(p, 0) + v
        return S(t, min(a.top, b.top))
    __radd__ = __add__

    def __neg__(a): return S({p: -v for p, v in a.t.items()}, a.top)
    def __sub__(a, b): return a + (-b if isinstance(b, S) else -sp.sympify(b))

    def __mul__(a, b):
        if not isinstance(b, S): return S({p: v*b for p, v in a.t.items()}, a.top)
        la, lb = a.low(), b.low(); top = min(a.top + lb, b.top + la); t = {}
        for p, v in a.t.items():
            for q, w in b.t.items():
                if p + q <= top: t[p+q] = t.get(p+q, 0) + v*w
        return S(t, top)
    __rmul__ = __mul__

    def inv(a):
        p = a.low(); ap = a.t[p]; n = a.top - p
        u = S({q - p: v/ap for q, v in a.t.items() if q > p}, n)
        r = S({0: sp.Integer(1)}, n); term = S({0: sp.Integer(1)}, n)
        for _ in range(n):
            term = term*(-u); r = r + term
        return S({q - p: v/ap for q, v in r.t.items()}, n - p)

    def d(a, x):
        return S({p-1: sp.diff(v, x) for p, v in a.t.items()}, a.top - 1)

    def simp(a):
        return S({p: sp.factor(sp.cancel(sp.together(v))) for p, v in a.t.items()}, a.top)

    def dump(a):
        return {f'eps^{p}': str(sp.factor(sp.cancel(sp.together(a.t[p])))) for p in sorted(a.t)}


def part_A(kextra):
    W, Wb = sp.symbols('W Wb'); c, cb, r0 = sp.symbols('c cb rho0')
    co, cob = {}, {}
    for n in (0, 1):
        for k in range(0, 2*n + 3 + kextra):
            co[(n, k)] = sp.Symbol(f'a{n}_{k}'); cob[(n, k)] = sp.Symbol(f'ab{n}_{k}')

    def mk(cs, X, Y, cc):
        return S({-1: cc/X, 0: sum(cs[(0, k)]*X**(0-k)*Y**k for k in range(3 + kextra)),
                  1: sum(cs[(1, k)]*X**(1-k)*Y**k for k in range(5 + kextra))}, 1)
    E = mk(co, W, Wb, c); Eb = mk(cob, Wb, W, cb)     # Eb = conj(E): conj(W^p Wb^q) = Wb^p W^q
    rho = S({0: r0, 1: (W + Wb)/2}, 4)
    half = sp.Rational(1, 2)
    sol = {}; rec = {}
    # Ernst equation times rho:  f (4 rho E_WWb + E_W + E_Wb) - 4 rho E_W E_Wb = 0, order by order
    for lev, n in ((-4, None), (-3, 0), (-2, 1)):
        eqs = []
        for EE, other in ((E, Eb), (Eb, E)):
            f = (EE + other)*half; EW, EWb = EE.d(W), EE.d(Wb)
            r = f*(4*rho*EW.d(Wb) + EW + EWb) - 4*rho*EW*EWb
            num = sp.numer(sp.together(sp.expand(r.t.get(lev, 0)).subs(sol)))
            if num != 0:
                eqs += sp.Poly(sp.expand(num), W, Wb).coeffs()
        if n is None:
            rec['order_eps^-4_identically_zero'] = (len(eqs) == 0)   # E0 = c/W holomorphic: planar Ernst eq. exact
            continue
        unk = [co[(n, k)] for k in range(2*n + 3 + kextra)] + [cob[(n, k)] for k in range(2*n + 3 + kextra)]
        sl = sp.solve(eqs, unk, dict=True)
        assert len(sl) == 1
        sol.update({k_: sp.simplify(v.subs(sol)) for k_, v in sl[0].items()})
        rec[f'degree{n}_free'] = [str(u) for u in unk if u not in sl[0]]
        rec[f'degree{n}_solution'] = {str(k_): str(sp.factor(v)) for k_, v in sl[0].items()}
    Es = S({p: v.subs(sol) for p, v in E.t.items()}, 1); Ebs = S({p: v.subs(sol) for p, v in Eb.t.items()}, 1)
    f = (Es + Ebs)*half; finv = f.inv()
    EW, EWb, EbW, EbWb = Es.d(W), Es.d(Wb), Ebs.d(W), Ebs.d(Wb)
    dlnH = rho*EW*EbW*finv*finv - (EW + EbW)*finv*half          # d_W ln H from the vacuum gamma equation
    QWW = ((EW.d(W) + EW*EW*finv*half - rho*EW*EW*EbW*finv*finv)*finv*half).simp()
    QbB = ((EWb.d(Wb) + EWb*EWb*finv*half - rho*EWb*EWb*EbWb*finv*finv)*finv*half).simp()
    QWWb = ((EW*EWb*finv*half - (EW + EWb)*rho.inv()*sp.Rational(1, 4))*finv*half).simp()
    X = QWW*QbB + 3*QWWb*QWWb
    rec['dlnH'] = dlnH.simp().dump(); rec['Q_WW'] = QWW.dump(); rec['Q_WbWb'] = QbB.dump()
    rec['Q_WWb'] = QWWb.dump(); rec['X'] = X.simp().dump()
    lead = sp.factor(sp.cancel(sp.together(X.t[-2])))
    target = 3*Wb**2*c**2/(4*W**4*cb**2*r0**2)
    rec['X_leading_eps^-2'] = str(lead)
    rec['X_leading_equals_3Wb^2c^2/(4W^4cb^2rho0^2)'] = bool(sp.simplify(lead - target) == 0)
    rec['X_orders_below_-2_all_zero'] = all(is_zero(X.t[p]) for p in X.t if p < -2)
    rec['free_symbols_in_X_leading'] = sorted(str(q) for q in lead.free_symbols)
    return rec


def kerr_symbolic():
    """Exact Kerr: local data at the ring r = 0, cos(theta) = 0 and the check I_lead = 3 Psi2^2 exactly."""
    M, a, eps, r1, c1 = sp.symbols('M a epsilon r1 c1', positive=True)
    r, cth = eps*r1, eps*c1
    sth = sp.sqrt(1 - cth**2)
    rho = sp.sqrt(r**2 - 2*M*r + a**2)*sth
    z = (r - M)*cth
    Wser = sp.series(rho - a + sp.I*z, eps, 0, 2).removeO()
    Wlin = sp.simplify(Wser.coeff(eps, 1))                      # W = Wlin*eps + O(eps^2)
    Wlin_expected = -(M/a)*(r1 + sp.I*a*c1)
    E = 1 - 2*M/(r - sp.I*a*cth)                                 # Kerr Ernst potential (BL)
    # residue in Wb: E ~ c/Wb with Wb = conj(W)
    Wb_lin = sp.conjugate(Wlin).subs({sp.conjugate(r1): r1, sp.conjugate(c1): c1, sp.conjugate(M): M,
                                     sp.conjugate(a): a})
    cres = sp.simplify(sp.limit(E*eps*Wb_lin, eps, 0))
    Sigma = r**2 + a**2*cth**2
    H = Sigma/((r - M)**2 - (M**2 - a**2)*cth**2)               # H = e^{2gamma}/f (Kerr, Weyl-Papapetrou)
    Lam = sp.simplify(sp.limit(H/(eps**2*Wlin*sp.conjugate(Wlin).subs({sp.conjugate(r1): r1,
                                                                        sp.conjugate(c1): c1,
                                                                        sp.conjugate(M): M,
                                                                        sp.conjugate(a): a})), eps, 0))
    rho0 = a
    Ilead = 3*(cres/cres)**2/(rho0**2*Lam**2*(Wb_lin*eps)**6)      # c real: phase factor 1; pole is in Wb (conj case)
    Iexact = 3*M**2/(r - sp.I*a*cth)**6                             # 3 Psi2^2, Psi2 = -M/(r - i a cos)^3
    ratio = sp.simplify(Ilead/Iexact)
    return dict(W_linear=str(Wlin), W_linear_is_minus_M_over_a_times_r_plus_i_a_c=bool(sp.simplify(Wlin - Wlin_expected) == 0),
                residue_c_in_Wbar=str(cres), Lambda=str(Lam), rho0=str(rho0),
                I_lead_over_3Psi2sq=str(ratio), exact=bool(ratio == 1),
                Kdelta6_equator_closed_form='48 M^2/k^6 (exact for all delta: K = 48 M^2/r^6, r = k delta)')


# =====================================================================================================================
# Part B/C: numerics
# =====================================================================================================================
x, y = sp.symbols('x y')


def build(M, a_, beta):
    """sympy E(x,y), H(x,y) and coordinate maps; lambdified to mpmath."""
    F = mn_functions(M, a_, beta)
    k, al = F['k'], F['al']
    beta = sp.nsimplify(beta)
    R = x*sp.sqrt(1 + (y**2 - 1)/x**2); u = x*y/R
    P = [sp.Integer(1), u, (3*u**2 - 1)/2, (5*u**3 - 3*u)/2]
    psi = beta*P[2]/R**3
    aa = -al*sp.exp(-2*beta*(-1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))))
    bb = al*sp.exp(2*beta*(1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))))
    Wc = (1 - sp.I*aa)*(1 - sp.I*bb)
    Ap = x*(1 + aa*bb) + sp.I*y*(bb - aa) + Wc
    Am = x*(1 + aa*bb) + sp.I*y*(bb - aa) - Wc
    E = sp.exp(2*psi)*Am/Ap
    H = F['e2g']/F['f']
    rho = k*sp.sqrt((x**2 - 1)*(1 - y**2)); z = k*x*y
    L = lambda e: sp.lambdify((x, y), e, modules='mpmath', cse=True)
    d = dict(k=mp.mpf(sp.N(k, DPS + 10)), al=mp.mpf(sp.N(al, DPS + 10)), beta=mp.mpf(sp.N(beta, DPS + 10)))
    d['E'] = [L(E), L(sp.diff(E, x)), L(sp.diff(E, y)), L(sp.diff(E, x, 2)), L(sp.diff(E, x, y)), L(sp.diff(E, y, 2))]
    d['H'] = [L(H), L(sp.diff(H, x)), L(sp.diff(H, y))]
    d['f9'] = L(F['f'])
    J = [[rho.diff(x), rho.diff(y)], [z.diff(x), z.diff(y)]]
    d['J'] = [[L(J[i][j]) for j in range(2)] for i in range(2)]
    d['dJ'] = [[[L(J[i][j].diff(v)) for v in (x, y)] for j in range(2)] for i in range(2)]
    d['rho'] = L(rho); d['z'] = L(z)
    d['Ap'] = [L(Ap), L(sp.diff(Ap, x)), L(sp.diff(Ap, y))]
    d['Am'] = L(Am); d['psi'] = L(psi)
    return d


def mat_inv(J):
    det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    return [[J[1][1]/det, -J[0][1]/det], [-J[1][0]/det, J[0][0]/det]]


def wderivs(d, X, Y):
    """E, E_W, E_Wb, E_WW, E_WbWb, E_WWb at (X, Y) via the exact chain rule (x,y) -> (rho, z)."""
    X, Y = mp.mpf(X), mp.mpf(Y)
    E0, Ex, Ey, Exx, Exy, Eyy = [g(X, Y) for g in d['E']]
    J = [[d['J'][i][j](X, Y) for j in range(2)] for i in range(2)]
    Mi = mat_inv(J)                                           # Mi[i][a] = d u_i / d v_a,  u = (x,y), v = (rho,z)
    gu = [Ex, Ey]; Hu = [[Exx, Exy], [Exy, Eyy]]
    gv = [sum(gu[i]*Mi[i][a_] for i in range(2)) for a_ in range(2)]
    dJ = [[[d['dJ'][i][j][e](X, Y) for e in range(2)] for j in range(2)] for i in range(2)]
    # d(Mi)/du_e = -Mi (dJ/du_e) Mi
    dMi_du = []
    for e in range(2):
        A = [[dJ[i][j][e] for j in range(2)] for i in range(2)]
        MA = [[sum(Mi[i][p]*A[p][q] for p in range(2)) for q in range(2)] for i in range(2)]
        dMi_du.append([[-sum(MA[i][q]*Mi[q][j] for q in range(2)) for j in range(2)] for i in range(2)])
    Hv = [[None, None], [None, None]]
    for a_ in range(2):
        for b_ in range(2):
            t1 = sum(Hu[i][j]*Mi[i][a_]*Mi[j][b_] for i in range(2) for j in range(2))
            t2 = sum(gu[i]*sum(dMi_du[e][i][a_]*Mi[e][b_] for e in range(2)) for i in range(2))
            Hv[a_][b_] = t1 + t2
    Er, Ez = gv; Err, Erz, Ezz = Hv[0][0], Hv[0][1], Hv[1][1]
    EW = (Er - 1j*Ez)/2; EWb = (Er + 1j*Ez)/2
    EWW = (Err - 2j*Erz - Ezz)/4; EWbWb = (Err + 2j*Erz - Ezz)/4; EWWb = (Err + Ezz)/4
    # d_W ln H from the metric (for the unsubstituted check)
    H0, Hx, Hy = [g(X, Y) for g in d['H']]
    Hr = Hx*Mi[0][0] + Hy*Mi[1][0]; Hz = Hx*Mi[0][1] + Hy*Mi[1][1]
    dWlnH = (Hr - 1j*Hz)/(2*H0)
    return dict(E=E0, EW=EW, EWb=EWb, EWW=EWW, EWbWb=EWbWb, EWWb=EWWb, H=H0, dWlnH=dWlnH,
                rho=d['rho'](X, Y), z=d['z'](X, Y), f9=d['f9'](X, Y))


def invariant(d, X, Y, form='substituted'):
    q = wderivs(d, X, Y)
    E, EW, EWb, EWW, EWbWb, EWWb, H, rho = (q[k_] for k_ in ('E', 'EW', 'EWb', 'EWW', 'EWbWb', 'EWWb', 'H', 'rho'))
    f = mp.re(E)
    EbW = mp.conj(EWb); EbWb = mp.conj(EW)                    # d_W conj(E) = conj(d_Wb E)
    if form == 'substituted':   # Ernst equation and vacuum gamma equation used (the form expanded in Part A)
        QWW = (EWW + EW**2/(2*f) - rho*EW**2*EbW/f**2)/(2*f)
        QbB = (EWbWb + EWb**2/(2*f) - rho*EWb**2*EbWb/f**2)/(2*f)
        QWWb = (EW*EWb/(2*f) - (EW + EWb)/(4*rho))/(2*f)
    else:                       # raw identity: metric's own d ln H and E_WWb; no field equation substituted
        dWlnH = q['dWlnH']; dWblnH = mp.conj(dWlnH)
        QWW = (EWW - dWlnH*EW - EW*EbW/(2*f))/(2*f)
        QbB = (EWbWb - dWblnH*EWb - EWb*EbWb/(2*f))/(2*f)
        QWWb = (EWWb - EW*EWb/(2*f))/(2*f)
    I = 4*(QWW*QbB + 3*QWWb**2)/H**2
    return I, q, dict(QWW=QWW, QbB=QbB, QWWb=QWWb)


def ernst_resid(d, X, Y):
    q = wderivs(d, X, Y)
    f = mp.re(q['E'])
    lhs = f*(4*q['EWWb'] + (q['EW'] + q['EWb'])/q['rho']); rhs = 4*q['EW']*q['EWb']
    gam_lhs = q['dWlnH'] + (q['EW'] + mp.conj(q['EWb']))/(2*f)       # = 2 gamma_W  (from the metric)
    gam_rhs = q['rho']*q['EW']*mp.conj(q['EWb'])/f**2                 # vacuum gamma equation
    return dict(ernst_rel=s(abs(lhs - rhs)/abs(rhs), 4), gamma_eq_rel=s(abs(gam_lhs - gam_rhs)/abs(gam_rhs), 4),
                ReE_minus_f9_rel=s(abs(f - q['f9'])/abs(f), 4))


def part_B():
    res = {}
    # (i) Kerr M = 1, a = 0.9 at generic points (both sheets)
    M, a_ = sp.Integer(1), sp.Rational(9, 10)
    d = build(M, a_, 0)
    k = d['k']; Mm = mp.mpf(1); am = mp.mpf('0.9')
    rows = []
    for (X, Y) in [('1.5', '0.3'), ('3.0', '0.7'), ('1.1', '-0.55'), ('-1.7', '0.4'), ('-3.2', '-0.25')]:
        X, Y = mp.mpf(X), mp.mpf(Y)
        I, q, _ = invariant(d, X, Y)
        I2, _, _ = invariant(d, X, Y, 'raw')
        r = Mm + k*X; c = Y; Sig = r**2 + am**2*c**2
        Kc = 48*Mm**2*(r**6 - 15*am**2*r**4*c**2 + 15*am**4*r**2*c**4 - am**6*c**6)/Sig**6
        Iex = 3*Mm**2/(r - 1j*am*c)**6
        rows.append(dict(x=s(X, 6), y=s(Y, 6), K_identity=s(16*mp.re(I), 25), K_closed=s(Kc, 25),
                         K_rel=s(abs(16*mp.re(I) - Kc)/abs(Kc), 4),
                         I_vs_3Psi2sq_rel=s(min(abs(I - Iex), abs(I - mp.conj(Iex)))/abs(Iex), 4),
                         raw_vs_substituted_rel=s(abs(I2 - I)/abs(I), 4), **ernst_resid(d, X, Y)))
    res['kerr_generic'] = rows
    # (ii),(iii) MN against Q2's stored direct-Riemann values
    q2 = json.load(open(os.path.join(ROOT, 'results', 'v15_q2_curvature.json')))
    pars = {'C0': (sp.Integer(1), sp.Rational(9, 10), sp.Rational(95, 100)/sp.Rational(19, 100)**sp.Rational(3, 2)),
            'F5-R': (sp.Integer(89), sp.Integer(80), sp.Rational(13394411, 1186380))}
    builds = {}
    for lab, (M, a_, b_) in pars.items():
        t = time.time(); d = build(M, a_, b_); builds[lab] = d
        xs = mp.mpf(q2['ring'][lab]['x_s'])
        xs = mp.findroot(lambda X: mp.re(d['Ap'][0](X, 0)), xs)
        rows = []
        for name, (sx, sy) in {'y=0, x>x_s': (1, 0), 'y=0, x<x_s': (-1, 0)}.items():
            for row in q2['ring'][lab]['approaches'][name]['rows']:
                dd = mp.mpf(row['d'])
                I, _, _ = invariant(d, xs + sx*dd, sy*dd)
                Kq2 = mp.mpf(row['K'])
                rows.append(dict(ray=name, d=row['d'], K_identity=s(16*mp.re(I), 30), K_Q2_direct=row['K'][:32],
                                 rel=s(abs(16*mp.re(I) - Kq2)/abs(Kq2), 4)))
        extra = {}
        for name in q2['ring'][lab]['approaches']:
            if name.startswith('y=0'): continue
            extra[name] = 'present in Q2 json (direction parametrisation not re-derived here; not compared)'
        res[f'MN_{lab}_vs_Q2'] = dict(build_s=round(time.time() - t, 1), rows=rows,
                                      max_rel=s(max(mp.mpf(r_['rel']) for r_ in rows), 4))
        if lab == 'C0':
            X, Y = mp.mpf(5), mp.mpf('0.3')
            I, _, _ = invariant(d, X, Y); I2, _, _ = invariant(d, X, Y, 'raw')
            Kq = mp.mpf('0.08357926111805873116582765')
            res['MN_C0_negative_control_x5_y0.3'] = dict(K_identity=s(16*mp.re(I), 28), K_Q2='0.08357926111805873116582765',
                                                         rel=s(abs(16*mp.re(I) - Kq)/Kq, 4),
                                                         raw_vs_substituted_rel=s(abs(I2 - I)/abs(I), 4),
                                                         **ernst_resid(d, X, Y))
            res['MN_C0_resid_near_ring_d1e-3'] = ernst_resid(d, xs + mp.mpf('1e-3'), mp.mpf('1e-3'))
            I, _, _ = invariant(d, xs + mp.mpf('1e-3'), mp.mpf('1e-3')); I2, _, _ = invariant(d, xs + mp.mpf('1e-3'), mp.mpf('1e-3'), 'raw')
            res['MN_C0_raw_vs_substituted_near_ring_d1e-3'] = s(abs(I2 - I)/abs(I), 4)
    return res, builds


def ring_data(d, xs, sheet=1):
    """Laurent data at an equatorial ring (y = 0) from the exact functions."""
    k, al, beta = d['k'], d['al'], d['beta']
    Ap, Apx, Apy = [g(xs, mp.mpf(0)) for g in d['Ap']]
    J = [[d['J'][i][j](xs, mp.mpf(0)) for j in range(2)] for i in range(2)]
    Mi = mat_inv(J)
    Apr = Apx*Mi[0][0] + Apy*Mi[1][0]; Apz = Apx*Mi[0][1] + Apy*Mi[1][1]
    lamW = (Apr - 1j*Apz)/2; lamWb = (Apr + 1j*Apz)/2              # A_+ ~ lamW W + lamWb Wb
    orient = 'holomorphic in W' if abs(lamWb) < abs(lamW)*mp.mpf(10)**(-DPS//2) else (
        'holomorphic in Wbar' if abs(lamW) < abs(lamWb)*mp.mpf(10)**(-DPS//2) else 'neither')
    lam = lamW if orient == 'holomorphic in W' else lamWb
    psi0 = d['psi'](xs, mp.mpf(0)) if beta != 0 else mp.mpf(0)
    Am0 = d['Am'](xs, mp.mpf(0))
    cres = mp.exp(2*psi0)*Am0/lam                                   # E ~ c/W (or c/Wbar)
    rho0 = d['rho'](xs, mp.mpf(0))
    # Lambda = lim H/|W|^2.  H = e^{Gamma-2psi} |A_+|^2/((x^2-y^2)(1-al^2)^2) exactly (Q2 Sec. B); at the ring:
    Hfun = d['H'][0]
    # closed form via the gamma exponent evaluated at the ring: Gamma0 = ln(H (x^2-y^2)(1-al^2)^2/|A_+|^2) + 2 psi;
    # evaluate the bracket at a nearby point and extrapolate is unnecessary: use the exact closed form of Q2 for MN,
    # and an independent numerical limit.
    Rs = mp.sqrt(xs**2 - 1)
    if beta != 0:
        Gam0 = -3*beta**2/(4*Rs**6) + 4*beta*xs/Rs - 2*beta*xs/Rs**3 - 4*beta
    else:
        Gam0 = mp.mpf(0)
    Lam_closed = mp.exp(Gam0 - 2*psi0)*abs(lam)**2/(xs**2*(1 - al**2)**2)
    Lam_num = []
    for (dx, dy) in [('1e-30', '0'), ('0', '1e-30'), ('-7e-31', '7e-31')]:
        X, Y = xs + mp.mpf(dx), mp.mpf(dy)
        Wv = (d['rho'](X, Y) - rho0) + 1j*d['z'](X, Y)
        Lam_num.append(Hfun(X, Y)/abs(Wv)**2)
    # the same Lambda written through the Ernst residue c = e^{2 psi0} A_-(p)/lam:
    Lam_res = mp.exp(Gam0 + 2*psi0)*abs(Am0)**2/(abs(cres)**2*xs**2*(1 - al**2)**2)
    return dict(lamW=lamW, lamWb=lamWb, orient=orient, lam=lam, psi0=psi0, Am0=Am0, c=cres, rho0=rho0, Rs=Rs,
                Gam0=Gam0, Lam=Lam_closed, Lam_num=Lam_num, k=k, Lam_res=Lam_res)


def richardson(fun, deltas):
    """fun(delta) = L + a delta^2 + b delta^4 + ...; deltas decreasing by factor 10."""
    vals = [fun(dl) for dl in deltas]
    R1 = [(100*vals[i+1] - vals[i])/99 for i in range(len(vals) - 1)]
    R2 = [(10000*R1[i+1] - R1[i])/9999 for i in range(len(R1) - 1)]
    return vals, R1, R2


def part_C(builds):
    res = {}
    rings = {}
    # MN rings (C0, F5-R from Part B builds) + third ring + Kerr ring
    rings['C0'] = (builds['C0'], mp.mpf('1.92982825852115000840276900535422580128'))
    rings['F5-R'] = (builds['F5-R'], mp.mpf('1.92192258791675977224705941968613828922'))
    d3 = build(sp.Integer(1), sp.Rational(3, 5), sp.Rational(1, 5))       # chi 0.6, q 0.1024 -> beta = 0.2
    rings['chi0.6_q0.1024'] = (d3, mp.mpf('1.0759622011947267412'))
    dk = build(sp.Integer(1), sp.Rational(9, 10), 0)                      # Kerr M = 1, a = 0.9; ring at x = -M/k
    rings['Kerr_ring_M1_a0.9'] = (dk, -1/dk['k'])
    for lab, (d, x0) in rings.items():
        xs = mp.findroot(lambda X: mp.re(d['Ap'][0](X, 0)), x0)
        L = ring_data(d, xs)
        k, rho0, Lam, c, Rs = L['k'], L['rho0'], L['Lam'], L['c'], L['Rs']
        phase2 = (c/mp.conj(c))**2
        coeff = 3*phase2/(rho0**2*Lam**2)                                  # I ~ coeff / W^6 (or conj orientation)
        dWdx = d['J'][0][0](xs, mp.mpf(0))                                 # on y = 0: W = (d rho/dx) delta + ...
        dWdy = d['J'][1][1](xs, mp.mpf(0))                                 # on x = x_s: W = i (dz/dy) y + ...
        Kd6_pred = 16*mp.re(coeff)/dWdx**6
        Ky6_pred = 16*mp.re(coeff/(1j)**6)/dWdy**6
        # numerical limits (symmetric average kills odd orders; Richardson in delta^2)
        Kx = lambda dl: (16*mp.re(invariant(d, xs + dl, 0)[0]) + 16*mp.re(invariant(d, xs - dl, 0)[0]))/2*dl**6
        Ky = lambda dl: 16*mp.re(invariant(d, xs, dl)[0])*dl**6
        dls = [mp.mpf(10)**(-e) for e in (5, 6, 7, 8)]
        vx, R1x, R2x = richardson(Kx, dls)
        vy, R1y, R2y = richardson(Ky, dls)
        # direct I*W^6 (or I*Wb^6) at the closest point: check of the complex coefficient and orientation
        X, Y = xs + mp.mpf('1e-9'), mp.mpf('1e-9')
        I, q, Qs = invariant(d, X, Y)
        Wv = (q['rho'] - rho0) + 1j*q['z']
        IW6 = I*Wv**6; IWb6 = I*mp.conj(Wv)**6
        # ingredient checks against the Part A expansion (real c): Q_WWb -> c Wb/(4 cbar rho0 W^2),
        # Q_WbWb -> 3 c Wb (cbar W + c Wb)/(16 cbar^2 rho0^2 W^2), Q_WW -> 3 c Wb/(W^2 (cbar W + c Wb))  [W-orientation]
        Wo = Wv if L['orient'] == 'holomorphic in W' else mp.conj(Wv)
        cc = c   # mirrored coordinates W' = Wbar keep the residue c (E = c/W')
        # in the conj orientation the roles of Q_WW and Q_WbWb swap (and c -> cbar); compare moduli for both
        pred_QWWb = cc*mp.conj(Wo)/(4*mp.conj(cc)*rho0*Wo**2)
        pred_QbB = 3*cc*mp.conj(Wo)*(mp.conj(cc)*Wo + cc*mp.conj(Wo))/(16*mp.conj(cc)**2*rho0**2*Wo**2)
        pred_QWW = 3*cc*mp.conj(Wo)/(Wo**2*(mp.conj(cc)*Wo + cc*mp.conj(Wo)))
        if L['orient'] == 'holomorphic in W':
            gotWW, gotbB = Qs['QWW'], Qs['QbB']
        else:
            gotWW, gotbB = Qs['QbB'], Qs['QWW']
        res[lab] = {
            'x_s': s(xs, 40), 'R_s': s(Rs, 30), 'k': s(k, 30), 'rho0': s(rho0, 30), 'pole_orientation': L['orient'],
            'dAplus_dW': s(L['lamW'], 25), 'dAplus_dWbar': s(L['lamWb'], 5), 'A_minus_at_ring': s(L['Am0'], 25),
            'psi0': s(L['psi0'], 25), 'Gamma0': s(L['Gam0'], 25), 'residue_c': s(c, 25), 'arg_c': s(mp.arg(c), 5),
            'residue_check_W_times_E_at_d1e-9': s(Wo*q['E'], 15),
            'Lambda_closed': s(Lam, 30), 'Lambda_via_residue_e^(G0+2psi0)Am^2/(|c|^2 xs^2 (1-al^2)^2)': s(L['Lam_res'], 30), 'Lambda_numeric_limit_3dirs_d1e-30': [s(v, 25) for v in L['Lam_num']],
            'coefficient_3(c/cbar)^2/(rho0^2 Lambda^2)': s(coeff, 25),
            'I_times_W6_at_d1e-9': s(IW6, 15), 'I_times_Wbar6_at_d1e-9': s(IWb6, 15),
            'Kdelta6_closed': s(Kd6_pred, 25), 'Kdelta6_numeric_symavg_d1e-5..1e-8': [s(v, 25) for v in vx],
            'Kdelta6_richardson1': [s(v, 25) for v in R1x], 'Kdelta6_richardson2': [s(v, 25) for v in R2x],
            'Kdelta6_rel_err': s(abs(R2x[-1] - Kd6_pred)/abs(Kd6_pred), 4),
            'Kdelta6_rel_err_raw_d1e-8_no_extrapolation': s(abs(vx[-1] - Kd6_pred)/abs(Kd6_pred), 4),
            'Ky6_closed': s(Ky6_pred, 25), 'Ky6_numeric_d1e-5..1e-8': [s(v, 25) for v in vy],
            'Ky6_richardson2': [s(v, 25) for v in R2y], 'Ky6_rel_err': s(abs(R2y[-1] - Ky6_pred)/abs(Ky6_pred), 4),
            'ingredient_rel_err_d1e-9': {'Q_WWb': s(abs(Qs['QWWb'] - pred_QWWb)/abs(pred_QWWb), 3),
                                         'Q_lead_pp': s(abs(gotWW - pred_QWW)/abs(pred_QWW), 3),
                                         'Q_second_order': s(abs(gotbB - pred_QbB)/abs(pred_QbB), 3)},
            'Kdelta6_closed_alt_48Rs4_over_Lam2k8xs6': s(48*Rs**4/(Lam**2*k**8*xs**6), 25),
        }
        if lab.startswith('Kerr'):
            res[lab]['Kerr_exact_Kdelta6_48M2_over_k6'] = s(48/k**6, 25)
            res[lab]['Kerr_exact_Lambda_a2_over_M4'] = s(mp.mpf('0.81'), 25)
            res[lab]['Kerr_exact_coefficient_3M8_over_a6'] = s(3/mp.mpf('0.9')**6, 25)
            res[lab]['Kerr_rel_err_closed_vs_exact'] = s(abs(Kd6_pred - 48/k**6)/(48/k**6), 4)
        print(lab, 'Kd6 closed', s(Kd6_pred, 18), 'numeric', s(R2x[-1], 18), 'rel', res[lab]['Kdelta6_rel_err'],
              '| Ky6 rel', res[lab]['Ky6_rel_err'], flush=True)
    # comparison with the numbers Q2 stored (finite-d, Q2's own Riemann code)
    q2 = json.load(open(os.path.join(ROOT, 'results', 'v15_q2_curvature.json')))
    cmp_ = {}
    for lab in ('C0', 'F5-R'):
        lc = q2['ring'][lab]['leading_coefficient']
        mine = mp.mpf(res[lab]['Kdelta6_closed'])
        for key in ('Kdelta6_numeric_sym_avg_d1em5', 'Kdelta6_predicted_y0'):
            cmp_[f'{lab}:{key}'] = dict(Q2=lc[key], rel_to_derived=s(abs(mp.mpf(lc[key]) - mine)/mine, 4))
        mine_y = mp.mpf(res[lab]['Ky6_closed'])
        cmp_[f'{lab}:Kd6_along_y_numeric_d1em5'] = dict(Q2=lc['Kd6_along_y_numeric_d1em5'],
                                                       rel_to_derived=s(abs(mp.mpf(lc['Kd6_along_y_numeric_d1em5']) - mine_y)/abs(mine_y), 4))
    t3 = q2['third_point_law_check']
    rows = t3['rows']
    sym = (mp.mpf(rows[4]['K']) + mp.mpf(rows[5]['K']))/2*mp.mpf('1e-6')**6
    mine = mp.mpf(res['chi0.6_q0.1024']['Kdelta6_closed'])
    cmp_['chi0.6_q0.1024:Q2_K_symavg_d1e-6_times_d6'] = dict(Q2=s(sym, 15), rel_to_derived=s(abs(sym - mine)/mine, 4))
    cmp_['chi0.6_q0.1024:Q2_N0_law'] = dict(Q2=t3['N0_law'],
                                            rel_to_derived=s(abs(mp.mpf(t3['N0_law']) - abs(mp.mpf(res['chi0.6_q0.1024']['coefficient_3(c/cbar)^2/(rho0^2 Lambda^2)'][0])))/mp.mpf(t3['N0_law']), 4))
    res['comparison_with_Q2_stored_numbers'] = cmp_
    return res


def part_D():
    """Off-equator pole (q < 0, Q1: chi 0.9, q -1 has a +-y pair of B = A = 0 rings near x 1.69, y 0.80).
    No reflection symmetry at the point, so arg c is generic: tests the phase factor (c/cbar)^2 and the
    angular law K |W|^6 -> 48 cos(6 theta - 4 arg c)/(rho0^2 Lambda^2). Not part of the registered PASS criterion."""
    M, a_ = sp.Integer(1), sp.Rational(9, 10)
    beta = -1/sp.Rational(19, 100)**sp.Rational(3, 2)                    # q = -1: beta = q M^3/k^3
    d = build(M, a_, beta)
    k = d['k']
    F = lambda X, Y: [mp.re(d['Ap'][0](X, Y)), mp.im(d['Ap'][0](X, Y))]
    x0, y0 = mp.findroot(F, (mp.mpf('1.69'), mp.mpf('0.80')))
    Ap, Apx, Apy = [g(x0, y0) for g in d['Ap']]
    J = [[d['J'][i][j](x0, y0) for j in range(2)] for i in range(2)]
    Mi = mat_inv(J)
    Apr = Apx*Mi[0][0] + Apy*Mi[1][0]; Apz = Apx*Mi[0][1] + Apy*Mi[1][1]
    lamW, lamWb = (Apr - 1j*Apz)/2, (Apr + 1j*Apz)/2
    holo = abs(lamWb) < abs(lamW)
    lam = lamW if holo else lamWb
    psi0 = d['psi'](x0, y0); Am0 = d['Am'](x0, y0)
    c = mp.exp(2*psi0)*Am0/lam
    rho0, z0 = d['rho'](x0, y0), d['z'](x0, y0)

    def xy(rho, z):
        rp, rm = mp.sqrt(rho**2 + (z + k)**2), mp.sqrt(rho**2 + (z - k)**2)
        return (rp + rm)/(2*k), (rp - rm)/(2*k)
    Lam_num = []
    for ang in (0, 1, 2):
        Wv = mp.mpf('1e-30')*mp.expj(ang)
        X, Y = xy(rho0 + mp.re(Wv), z0 + mp.im(Wv))
        Lam_num.append(d['H'][0](X, Y)/abs(Wv)**2)
    Lam = Lam_num[0]
    coeff = 3*(c/mp.conj(c))**2/(rho0**2*Lam**2)
    rows = []
    for th in [mp.mpf(j)*mp.pi/7 for j in range(7)]:
        vals = []
        for t in (mp.mpf('1e-6'), mp.mpf('1e-7'), mp.mpf('1e-8')):
            Wv = t*mp.expj(th)
            X, Y = xy(rho0 + mp.re(Wv), z0 + mp.im(Wv))
            I, q, _ = invariant(d, X, Y)
            Wo = Wv if holo else mp.conj(Wv)
            vals.append((I*Wo**6, 16*mp.re(I)*t**6))
        IW6 = (10*vals[2][0] - vals[1][0])/9; Kt6 = (10*vals[2][1] - vals[1][1])/9   # linear Richardson (O(t) residual)
        thW = th if holo else -th
        Kpred = 48*mp.cos(6*thW - 4*mp.arg(c))/(rho0**2*Lam**2)
        rows.append(dict(theta=s(th, 6), I_times_W6_extrap=s(IW6, 15), rel_to_coeff=s(abs(IW6 - coeff)/abs(coeff), 3),
                         K_t6_extrap=s(Kt6, 15), K_pred_48cos=s(Kpred, 15),
                         K_abs_err_over_scale=s(abs(Kt6 - Kpred)/(48/(rho0**2*Lam**2)), 3)))
    return dict(params='M=1, a=0.9, q=-1 (beta=-1/0.19^1.5)', x0=s(x0, 30), y0=s(y0, 30), rho0=s(rho0, 25), z0=s(z0, 25),
                A_plus_at_point=s(abs(Ap), 3), pole_orientation='holomorphic in W' if holo else 'holomorphic in Wbar',
                dAplus_dW=s(lamW, 20), dAplus_dWbar=s(lamWb, 5), A_minus_at_point=s(Am0, 20), residue_c=s(c, 20),
                arg_c=s(mp.arg(c), 20), Lambda_numeric_3dirs=[s(v, 20) for v in Lam_num],
                rows=rows, coefficient=s(coeff, 20),
                max_rel_IW6=s(max(mp.mpf(r_['rel_to_coeff']) for r_ in rows), 3),
                max_K_err_over_scale=s(max(mp.mpf(r_['K_abs_err_over_scale']) for r_ in rows), 3))


def main():
    print('Part A ...', flush=True)
    t = time.time()
    OUT['A_series_ansatz_k_le_2n+2'] = part_A(0)
    OUT['A_series_ansatz_k_le_2n+6'] = part_A(4)
    OUT['A_kerr_exact'] = kerr_symbolic()
    OUT['A_time_s'] = round(time.time() - t, 1)
    print(json.dumps({k_: OUT['A_series_ansatz_k_le_2n+2'][k_] for k_ in ('X_leading_eps^-2', 'degree0_free', 'degree1_free')}),
          OUT['A_kerr_exact']['exact'], flush=True)
    print('Part B ...', flush=True)
    t = time.time()
    OUT['B_identity_checks'], builds = part_B()
    OUT['B_time_s'] = round(time.time() - t, 1)
    print('Part C ...', flush=True)
    t = time.time()
    OUT['C_rings'] = part_C(builds)
    OUT['C_time_s'] = round(time.time() - t, 1)
    print('Part D ...', flush=True)
    t = time.time()
    OUT['D_offequator_phase_test'] = part_D()
    OUT['D_time_s'] = round(time.time() - t, 1)
    print('D', OUT['D_offequator_phase_test']['arg_c'], OUT['D_offequator_phase_test']['max_rel_IW6'],
          OUT['D_offequator_phase_test']['max_K_err_over_scale'], flush=True)
    rels = [mp.mpf(OUT['C_rings'][lab]['Kdelta6_rel_err']) for lab in ('C0', 'F5-R', 'chi0.6_q0.1024', 'Kerr_ring_M1_a0.9')]
    OUT['PASS_criterion'] = dict(rule='derived coefficient reproduces numeric K delta^6 at the 4 rings to <= 1e-6 relative',
                                 max_rel=s(max(rels), 4), PASS=bool(max(rels) <= mp.mpf('1e-6')),
                                 kerr_symbolic_exact=OUT['A_kerr_exact']['exact'])
    OUT['runtime_s'] = round(time.time() - T0, 1)
    with open(os.path.join(ROOT, 'results', 'v15_q2b_check.json'), 'w') as fh:
        json.dump(OUT, fh, indent=1)
    print('PASS' if OUT['PASS_criterion']['PASS'] else 'FAIL', OUT['PASS_criterion'], flush=True)


if __name__ == '__main__':
    main()
