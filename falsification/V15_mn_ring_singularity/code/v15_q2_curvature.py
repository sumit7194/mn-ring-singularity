#!/usr/bin/env python
"""V15 Q2: curvature of the Manko-Novikov (MN) metric near the equatorial B = A = 0 point.

NEW ROUTE (independent of V13 v13_kretschmann.py and of tabula's script, neither of which was read):
  * The metric is treated as an orthogonally transitive 2+2 split
        ds^2 = G_AB(x,y) dy^A dy^B + h_ij(x,y) dx^i dx^j ,   y^A = (t, phi),  x^i = (x, y),
    with G = [[-f, f w], [f w, rho^2/f - f w^2]],  h = diag(H, H (x^2-1)/(1-y^2)),
    H = e^{2 gamma} k^2 (x^2-y^2) / (f (x^2-1)).
  * The Riemann tensor is assembled directly from closed-form block formulas (derived in
    results/v15_q2_analytic.md, Sec. A), not from a generic Christoffel -> Riemann loop:
        R_ijkl = K_h (h_ik h_jl - h_il h_jk)                 (K_h = Gauss curvature of h, Brioschi)
        R_ABCD = -1/4 h^ij (dG_i,AC dG_j,BD - dG_i,AD dG_j,BC)   (Gauss equation for the torus orbits)
        R_ABij = -1/4 (dG_i G^-1 dG_j - dG_j G^-1 dG_i)_AB     (normal curvature)
        R_iAjB = -1/2 (D_i D_j G)_AB + 1/4 (dG_j G^-1 dG_i)_AB  (D = Levi-Civita of h on scalars)
    all odd-fibre-index components vanish identically.
  * Invariants: K = R_abcd R^abcd, Pontryagin K* = *R_abcd R^abcd, the complex Weyl invariants
    I = (K - i K*)/16, J = (C C C - i C C *C)/96 and the speciality index S = 27 J^2 / I^3 (=1 for
    type D); Ricci tensor (vacuum check) and first-Bianchi residual (self-consistency check).
  * Exact symbolic derivatives (sympy) of G and h to 2nd order, lambdified to mpmath; dps >= 100.
Metric functions: GLM 2008 arXiv:0708.0628 Eqs 2-7 / Brink 2008 arXiv:0807.1179 App. C, as transcribed
in V9_mn_obstruction_check/code/v9_nve.py (re-typed here from that transcription).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))   # snapshot: paths relative to this file
import json, sys, time, platform, itertools
import sympy as sp
import mpmath as mp

OUT = _os.path.join(_HERE, '..', 'results', 'v15_q2_curvature.json')
DPS = 100
mp.mp.dps = DPS
x, y = sp.symbols('x y', real=True)


# ---------------------------------------------------------------- metric (MN, prolate x, y)
def mn_symbolic(M, a, beta):
    M, a, beta = map(sp.nsimplify, (M, a, beta))
    chi = a / M
    al = (-1 + sp.sqrt(1 - chi**2)) / chi
    k = sp.sqrt(M**2 - a**2)
    R = sp.sqrt(x**2 + y**2 - 1)
    u = x * y / R
    P = [sp.Integer(1), u, (3*u**2 - 1)/2, (5*u**3 - 3*u)/2]
    psi = beta * P[2] / R**3
    aa = -al * sp.exp(-2*beta*(-1 + sum((x - y)*P[l]/R**(l+1) for l in range(3))))
    bb = al * sp.exp(2*beta*(1 + sum((-1)**(3-l)*(x + y)*P[l]/R**(l+1) for l in range(3))))
    A = (x**2 - 1)*(1 + aa*bb)**2 - (1 - y**2)*(bb - aa)**2
    B = (x + 1 + (x - 1)*aa*bb)**2 + ((1 + y)*aa + (1 - y)*bb)**2
    C = (x**2 - 1)*(1 + aa*bb)*(bb - aa - y*(aa + bb)) + (1 - y**2)*(bb - aa)*(1 + aa*bb + x*(1 - aa*bb))
    f = sp.exp(2*psi) * A / B
    w = 2*k*sp.exp(-2*psi)*C/A - 4*k*al/(1 - al**2)
    # e^{2 gamma}; the "-2" sits inside the l-sum (3 terms -> -6), which makes e^{2gamma} -> 1 at infinity
    gexp = 3*beta**2*(P[3]**2 - P[2]**2)/R**6 \
        + 2*beta*sum((x - y + (-1)**l*(x + y))*P[l]/R**(l+1) - 2 for l in range(3)) + 8*beta
    e2g = sp.exp(gexp) * A / ((x**2 - 1)*(1 - al**2)**2) * (x**2 - 1)/(x**2 - y**2)
    return dict(f=f, w=w, e2g=e2g, k=k, al=al, A=A, B=B, C=C, a=aa, b=bb, psi=psi)


class Spacetime:
    """Lambdified G_AB, h_ii and their first and second (x, y) derivatives."""

    def __init__(self, M, a, beta, label):
        t0 = time.time()
        self.label = label
        F = mn_symbolic(M, a, beta)
        k = F['k']; f = F['f']; w = F['w']; e2g = F['e2g']
        rho2 = k**2 * (x**2 - 1) * (1 - y**2)
        Gs = [[-f, f*w], [f*w, rho2/f - f*w**2]]
        H = e2g * k**2 * (x**2 - y**2) / (f * (x**2 - 1))
        hs = [H, H * (x**2 - 1) / (1 - y**2)]
        V = [x, y]
        exprs = []
        # layout: G00,G01,G11 ; h0,h1 ; then d/dx, d/dy of each ; then second derivs xx, xy, yy
        base = [Gs[0][0], Gs[0][1], Gs[1][1], hs[0], hs[1]]
        exprs += base
        d1 = [[sp.diff(e, v) for v in V] for e in base]
        d2 = [[sp.diff(d1[n][0], x), sp.diff(d1[n][0], y), sp.diff(d1[n][1], y)] for n in range(5)]
        for n in range(5):
            exprs += d1[n]
        for n in range(5):
            exprs += d2[n]
        self.fun = sp.lambdify((x, y), exprs, modules='mpmath', cse=True)
        aa, bb = F['a'], F['b']
        Wc = (1 - sp.I*aa)*(1 - sp.I*bb)
        Ap = x*(1 + aa*bb) + sp.I*y*(bb - aa) + Wc        # B = |A+|^2, A = Re(A- conj(A+))
        Am_ = x*(1 + aa*bb) + sp.I*y*(bb - aa) - Wc
        Ern = sp.exp(2*F['psi'])*Am_/Ap                    # Ernst potential candidate, Re = f
        aux = [F['A'], F['B'], F['C'], aa, bb, f, w, e2g, x + 1 + (x - 1)*aa*bb,
               (1 + y)*aa + (1 - y)*bb, Ap, Am_, F['psi'], sp.sqrt(rho2), k*x*y]
        self.auxfun = sp.lambdify((x, y), aux, modules='mpmath')
        # Ernst equation (prolate form) residual, normalised
        Ex_, Ey_ = sp.diff(Ern, x), sp.diff(Ern, y)
        lhs = (Ern + sp.conjugate(Ern))/2*(sp.diff((x**2 - 1)*Ex_, x) + sp.diff((1 - y**2)*Ey_, y))
        rhs = (x**2 - 1)*Ex_**2 + (1 - y**2)*Ey_**2
        self.ernst_res = sp.lambdify((x, y), [lhs, rhs], modules='mpmath')
        self.beta = mp.mpf(sp.N(sp.nsimplify(beta), DPS + 20))
        self.k = mp.mpf(sp.N(k, DPS + 20)); self.al = mp.mpf(sp.N(F['al'], DPS + 20))
        self.Hfun = sp.lambdify((x, y), sp.sqrt(H), modules='mpmath')
        self.build_s = time.time() - t0

    def aux(self, X, Y):
        v = self.auxfun(mp.mpf(X), mp.mpf(Y))
        dct = dict(zip(['A', 'B', 'C', 'a', 'b', 'f', 'w', 'e2g', 'B1', 'B2', 'Ap', 'Am', 'psi', 'rho', 'z'], v))
        dct['Ern'] = mp.exp(2*dct['psi'])*dct['Am']/dct['Ap'] if dct['Ap'] != 0 else mp.inf
        return dct

    def jets(self, X, Y):
        v = self.fun(mp.mpf(X), mp.mpf(Y))
        G = [v[0], v[1], v[2]]; h = [v[3], v[4]]
        d = [[v[5 + 2*n], v[6 + 2*n]] for n in range(5)]
        dd = [[v[15 + 3*n], v[16 + 3*n], v[17 + 3*n]] for n in range(5)]
        return G, h, d, dd


# ---------------------------------------------------------------- 2+2 block curvature
def riemann_2p2(G, h, d, dd):
    """Covariant Riemann R[a][b][c][d], index order 0=t, 1=phi (fibre), 2=x, 3=y (base)."""
    mat = lambda c00, c01, c11: mp.matrix([[c00, c01], [c01, c11]])
    Gm = mat(*G)
    Gi = Gm**-1
    dG = [mat(d[0][i], d[1][i], d[2][i]) for i in range(2)]
    def d2(n, i, j):
        idx = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 2}[(i, j)]
        return dd[n][idx]
    ddG = [[mat(d2(0, i, j), d2(1, i, j), d2(2, i, j)) for j in range(2)] for i in range(2)]
    E, Gb = h
    Ex, Ey = d[3]; Gx, Gy = d[4]
    Eyy = dd[3][2]; Gxx = dd[4][0]
    hinv = [1/E, 1/Gb]
    # Christoffels of diagonal 2D metric h = diag(E, Gb), Gam[k][i][j]
    Gam = [[[Ex/(2*E), Ey/(2*E)], [Ey/(2*E), -Gx/(2*E)]],
           [[-Ey/(2*Gb), Gx/(2*Gb)], [Gx/(2*Gb), Gy/(2*Gb)]]]
    # Gauss curvature (orthogonal-coordinates formula)
    W = E*Gb; Wx = Ex*Gb + E*Gx; Wy = Ey*Gb + E*Gy; sW = mp.sqrt(W)
    Kh = -(Gxx/sW - Gx*Wx/(2*W*sW) + Eyy/sW - Ey*Wy/(2*W*sW)) / (2*sW)
    hm = [[E, 0], [0, Gb]]
    R = [[[[mp.mpf(0)]*4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    # base-base
    for i, j, kk, l in itertools.product(range(2), repeat=4):
        R[2+i][2+j][2+kk][2+l] = Kh*(hm[i][kk]*hm[j][l] - hm[i][l]*hm[j][kk])
    # fibre-fibre
    for A, B, C, D in itertools.product(range(2), repeat=4):
        s = 0
        for i in range(2):
            s += hinv[i]*(dG[i][A, C]*dG[i][B, D] - dG[i][A, D]*dG[i][B, C])
        R[A][B][C][D] = -s/4
    # R_ABij
    for i, j in itertools.product(range(2), repeat=2):
        Nm = (dG[i]*Gi*dG[j] - dG[j]*Gi*dG[i]) * (-mp.mpf(1)/4)
        for A, B in itertools.product(range(2), repeat=2):
            R[A][B][2+i][2+j] = Nm[A, B]
            R[2+i][2+j][A][B] = Nm[A, B]
    # R_iAjB
    for i, j in itertools.product(range(2), repeat=2):
        DD = ddG[i][j] - Gam[0][i][j]*dG[0] - Gam[1][i][j]*dG[1]
        X = DD*(-mp.mpf(1)/2) + (dG[j]*Gi*dG[i])/4
        for A, B in itertools.product(range(2), repeat=2):
            v = X[A, B]
            R[2+i][A][2+j][B] = v
            R[A][2+i][B][2+j] = v
            R[2+i][A][B][2+j] = -v
            R[A][2+i][2+j][B] = -v
    g = [[G[0], G[1], 0, 0], [G[1], G[2], 0, 0], [0, 0, E, 0], [0, 0, 0, Gb]]
    gi = [[Gi[0, 0], Gi[0, 1], 0, 0], [Gi[1, 0], Gi[1, 1], 0, 0], [0, 0, 1/E, 0], [0, 0, 0, 1/Gb]]
    detg = (G[0]*G[2] - G[1]**2)*E*Gb
    return R, g, gi, detg


def raise_all(R, gi):
    """R^{ab}_{cd} (mixed, first pair raised) and R^{abcd}."""
    n = range(4)
    def raise_idx(T, pos):
        out = [[[[mp.mpf(0)]*4 for _ in n] for _ in n] for _ in n]
        for a, b, c, d in itertools.product(n, repeat=4):
            idx = [a, b, c, d]
            s = 0
            for e in n:
                if gi[idx[pos]][e] == 0:
                    continue
                j = list(idx); j[pos] = e
                s += gi[idx[pos]][e]*T[j[0]][j[1]][j[2]][j[3]]
            out[a][b][c][d] = s
        return out
    Rm = raise_idx(raise_idx(R, 0), 1)       # R^{ab}_{cd}
    Ru = raise_idx(raise_idx(Rm, 2), 3)      # R^{abcd}
    return Rm, Ru


def levi(a, b, c, d):
    p = [a, b, c, d]
    if len(set(p)) < 4:
        return 0
    s = 1
    for i in range(4):
        for j in range(i+1, 4):
            if p[i] > p[j]:
                s = -s
    return s


def invariants(st, X, Y):
    G, h, d, dd = st.jets(X, Y)
    R, g, gi, detg = riemann_2p2(G, h, d, dd)
    Rm, Ru = raise_all(R, gi)
    n = range(4)
    K = mp.fsum(R[a][b][c][e]*Ru[a][b][c][e] for a, b, c, e in itertools.product(n, repeat=4))
    # Ricci R_bd = g^{ac} R_abcd ; scalar invariants of Ricci
    Ric = [[mp.fsum(gi[a][c]*R[a][b][c][e] for a in n for c in n) for e in n] for b in n]
    RicU = [[mp.fsum(gi[b][p]*gi[e][q]*Ric[p][q] for p in n for q in n) for e in n] for b in n]
    RicSq = mp.fsum(Ric[b][e]*RicU[b][e] for b in n for e in n)
    # first Bianchi residual (scaled by sqrt K)
    bian = max(abs(R[a][b][c][e] + R[a][c][e][b] + R[a][e][b][c]) for a, b, c, e in itertools.product(n, repeat=4))
    # dual: *R_{abcd} = 1/2 eps_{abef} R^{ef}_{cd}, eps_{0123} = sqrt(-g)
    sg = mp.sqrt(-detg)
    Rd = [[[[mp.mpf(0)]*4 for _ in n] for _ in n] for _ in n]
    for a, b in itertools.product(n, repeat=2):
        if a == b:
            continue
        pairs = [(e, f_) for e in n for f_ in n if levi(a, b, e, f_) != 0]
        for c, e2 in itertools.product(n, repeat=2):
            Rd[a][b][c][e2] = sg/2*mp.fsum(levi(a, b, e, f_)*Rm[e][f_][c][e2] for e, f_ in pairs)
    Kst = mp.fsum(Rd[a][b][c][e]*Ru[a][b][c][e] for a, b, c, e in itertools.product(n, repeat=4))
    # cubic: R_ab^cd R_cd^ef R_ef^ab and with one dual. Need R_{ab}^{cd} = lower first pair of Rm... use
    # Q[(ab)][(cd)] = R^{ab}_{cd}; Qd[(ab)][(cd)] = (*R)^{ab}_{cd}
    Rdm = [[[[mp.fsum(gi[a][p]*gi[b][q]*Rd[p][q][c][e] for p in n for q in n if gi[a][p] != 0 and gi[b][q] != 0)
              for e in n] for c in n] for b in n] for a in n]
    pairs = [(a, b) for a in n for b in n]
    cub = mp.fsum(Rm[a][b][c][e]*Rm[c][e][p][q]*Rm[p][q][a][b]
                  for (a, b) in pairs for (c, e) in pairs for (p, q) in pairs)
    cubd = mp.fsum(Rm[a][b][c][e]*Rm[c][e][p][q]*Rdm[p][q][a][b]
                   for (a, b) in pairs for (c, e) in pairs for (p, q) in pairs)
    I = mp.mpc(K, -Kst)/16
    J = mp.mpc(cub, -cubd)/96
    S = 27*J**2/I**3
    return dict(K=K, Kstar=Kst, I=I, J=J, S=S, Psi2_typeD=-3*J/I, RicSq=RicSq,
                Ric_max=max(abs(Ric[b][e]) for b in n for e in n), bianchi_max=bian,
                g=g, h=h)


# ---------------------------------------------------------------- helpers
def s(v, n=25):
    if isinstance(v, mp.mpc):
        return [mp.nstr(v.real, n), mp.nstr(v.imag, n)]
    return mp.nstr(v, n)


def kerr_K_closed(M, a, k, X, Y):
    r = M + k*X; c = Y
    Sig = r**2 + a**2*c**2
    K = 48*M**2*(r**6 - 15*a**2*r**4*c**2 + 15*a**4*r**2*c**4 - a**6*c**6)/Sig**6
    Kst = 96*M**2*a*r*c*(3*r**2 - a**2*c**2)*(r**2 - 3*a**2*c**2)/Sig**6
    Psi2 = -M/(r - 1j*a*c)**3
    return K, Kst, Psi2


def digits(a_, b_):
    if a_ == b_:
        return mp.mpf(DPS)
    return -mp.log10(abs(a_ - b_)/abs(b_))


def fit_slope(ds, Ks):
    X = [mp.log(dd) for dd in ds]; Yv = [mp.log(abs(K)) for K in Ks]
    n = len(X); mx = sum(X)/n; my = sum(Yv)/n
    return sum((X[i]-mx)*(Yv[i]-my) for i in range(n))/sum((X[i]-mx)**2 for i in range(n))


def main():
    res = dict(meta=dict(python=sys.version, sympy=sp.__version__, mpmath=mp.__version__,
                         platform=platform.platform(), dps=DPS,
                         route='2+2 orthogonally-transitive block Riemann (closed-form block formulas), '
                               'prolate (x,y) base, exact sympy jets -> mpmath'))
    T0 = time.time()
    # ---------------- positive control: Kerr (beta = 0)
    Mk, ak = 1, sp.Rational(9, 10)
    kerr = Spacetime(Mk, ak, 0, 'Kerr M=1 a=0.9')
    kpts = [(mp.mpf('1.5'), mp.mpf('0.3')), (mp.mpf('3.0'), mp.mpf('0.7')), (mp.mpf('1.1'), mp.mpf('-0.55'))]
    kres = []
    for X, Y in kpts:
        inv = invariants(kerr, X, Y)
        Kc, Ksc, P2c = kerr_K_closed(mp.mpf(1), mp.mpf('0.9'), kerr.k, X, Y)
        kres.append(dict(x=s(X, 6), y=s(Y, 6), K=s(inv['K']), K_closed=s(Kc), K_digits=s(digits(inv['K'], Kc), 5),
                         Kstar=s(inv['Kstar']), Kstar_closed=s(Ksc),
                         Kstar_digits_absval=s(digits(abs(inv['Kstar']), abs(Ksc)), 5),
                         Kstar_sign_rel=s(mp.sign(inv['Kstar']*Ksc), 3),
                         S_minus_1=s(abs(inv['S'] - 1), 5),
                         absPsi2=s(abs(inv['Psi2_typeD'])), absPsi2_closed=s(abs(P2c)),
                         absPsi2_digits=s(digits(abs(inv['Psi2_typeD']), abs(P2c)), 5),
                         RicSq_over_K=s(abs(inv['RicSq'])/abs(inv['K']), 5),
                         bianchi_max=s(inv['bianchi_max'], 5)))
        print('Kerr', kres[-1], flush=True)
    res['positive_control_kerr'] = dict(params='M=1, a=0.9, beta=0', points=kres, build_s=kerr.build_s)
    # Kerr's own ring: B = 0 on y = 0 at x = -M/k (r = 0); same local law as used for MN (analytic note Sec. D)
    xr = mp.findroot(lambda X: kerr.aux(X, 0)['B1'], -1/kerr.k + mp.mpf('1e-3'))
    rho0k = kerr.aux(xr, 0)['rho']
    kring = dict(x_ring=s(xr, 30), x_ring_closed_minus_M_over_k=s(xr + 1/kerr.k, 5), rows=[])
    for dd in [mp.mpf('1e-3'), mp.mpf('1e-5')]:
        for (X, Y) in [(xr + dd, 0), (xr, dd), (xr - dd/mp.sqrt(2), dd/mp.sqrt(2))]:
            inv = invariants(kerr, X, Y); q = kerr.aux(X, Y)
            Wh = (q['rho'] - rho0k) + 1j*q['z']
            Kc, _, _ = kerr_K_closed(mp.mpf(1), mp.mpf('0.9'), kerr.k, X, Y)
            kring['rows'].append(dict(x=s(X, 12), y=s(Y, 12), K=s(inv['K'], 20), K_closed_digits=s(digits(inv['K'], Kc), 5),
                                      I_times_conjW6=s(inv['I']*mp.conj(Wh)**6, 12), I_times_W6=s(inv['I']*Wh**6, 12),
                                      absI_times_absW6=s(abs(inv['I'])*abs(Wh)**6, 12),
                                      Lambda=s(q['e2g']/(q['f']*abs(Wh)**2), 12)))
    c1k = mp.diff(lambda X: kerr.aux(X, 0)['Ap'], xr); c2k = mp.diff(lambda Y: kerr.aux(xr, Y)['Ap'], 0)
    kring['Aplus_coeffs'] = dict(dAp_dx=s(c1k, 15), dAp_dy=s(c2k, 15), Rs=s(mp.sqrt(xr**2 - 1), 15),
                                 dz_dy_sign=s(mp.sign(xr), 2), drho_dx_sign=s(mp.sign(xr), 2))
    kring['N0_exact_3M8_over_a6'] = s(3/mp.mpf('0.9')**6, 12)
    kring['Lambda_exact_a2_over_M4'] = s(mp.mpf('0.81'), 12)
    res['positive_control_kerr_ring'] = kring
    print('kerr ring', kring, flush=True)

    # ---------------- MN points
    pars = {'C0': (1, sp.Rational(9, 10), sp.Rational(95, 100)/sp.Rational(19, 100)**sp.Rational(3, 2),
                   '1.92982825852'),
            'F5-R': (89, 80, sp.Rational(13394411, 1186380), '1.92192258792')}
    sts = {}
    for lab, (M, a, beta, xs_given) in pars.items():
        sts[lab] = Spacetime(M, a, beta, lab)
        print(lab, 'built in', sts[lab].build_s, flush=True)

    # ---------------- negative control: C0 at (5, 0.3) and small displacements
    st = sts['C0']
    neg = []
    base = invariants(st, mp.mpf(5), mp.mpf('0.3'))
    for dx, dy in [(0, 0), (mp.mpf('1e-6'), 0), (-mp.mpf('1e-6'), 0), (0, mp.mpf('1e-6')), (0, -mp.mpf('1e-6')),
                   (mp.mpf('1e-3'), mp.mpf('1e-3'))]:
        inv = base if (dx, dy) == (0, 0) else invariants(st, 5 + dx, mp.mpf('0.3') + dy)
        neg.append(dict(dx=s(dx, 3), dy=s(dy, 3), K=s(inv['K']), relchange=s((inv['K'] - base['K'])/base['K'], 6),
                        RicSq_over_K=s(abs(inv['RicSq'])/abs(inv['K']), 5), S=s(inv['S'], 12)))
    # smoothness: symmetric second difference ratio should be O(h^2)
    Kp = mp.mpf(neg[1]['K']); Km = mp.mpf(neg[2]['K']); K0 = base['K']
    res['negative_control_C0_x5_y0.3'] = dict(points=neg,
                                               second_diff_x_rel=s((Kp - 2*K0 + Km)/K0, 6),
                                               first_diff_x_rel=s((Kp - Km)/K0, 6))
    print('neg', res['negative_control_C0_x5_y0.3'], flush=True)

    # ---------------- ring
    ds = [mp.mpf(10)**(-e) for e in (2, 3, 4, 5)]
    rings = {}
    for lab, (M, a, beta, xs_given) in pars.items():
        st = sts[lab]
        B1 = lambda X: st.aux(X, 0)['B1']
        xs = mp.findroot(B1, mp.mpf(xs_given))
        au = st.aux(xs, 0)
        info = dict(x_s=s(xs, 40), x_s_given=xs_given, x_s_minus_given=s(xs - mp.mpf(xs_given), 5),
                    at_ring=dict(A=s(au['A'], 5), B=s(au['B'], 5), a_plus_b=s(au['a'] + au['b'], 5),
                                 C=s(au['C'], 8)))
        # scaling of A, B, f, e2g, g_phiphi near ring (y = 0, right side)
        scal = []
        for dd in ds:
            aup = st.aux(xs + dd, 0)
            G, hh, _, _ = st.jets(xs + dd, 0)
            scal.append(dict(d=s(dd, 3), A=s(aup['A'], 8), B=s(aup['B'], 8), f=s(aup['f'], 8), C=s(aup['C'], 8),
                             g_tt=s(G[0], 8), g_tphi=s(G[1], 8), g_phiphi=s(G[2], 8), h_xx=s(hh[0], 8)))
        info['scaling_y0_right'] = scal
        # exact-structure checks (Sec. B of the analytic note)
        chk = []
        for (X, Y) in [(xs + mp.mpf('0.05'), mp.mpf('0.2')), (mp.mpf(3), mp.mpf('0.6')), (xs - mp.mpf('0.01'), -mp.mpf('0.1'))]:
            q = st.aux(X, Y); er = st.ernst_res(X, Y)
            chk.append(dict(x=s(X, 8), y=s(Y, 8),
                            B_minus_absAp2_rel=s(abs(q['B'] - abs(q['Ap'])**2)/q['B'], 5),
                            A_minus_ReAmconjAp_rel=s(abs(q['A'] - (q['Am']*mp.conj(q['Ap'])).real)/abs(q['A']), 5),
                            f_minus_ReErnst_rel=s(abs(q['f'] - q['Ern'].real)/abs(q['f']), 5),
                            ernst_eq_residual_rel=s(abs(er[0] - er[1])/abs(er[1]), 5)))
        info['ernst_structure_checks'] = chk
        Apf = lambda X, Y: st.aux(X, Y)['Ap']
        c1 = mp.diff(lambda X: Apf(X, 0), xs)
        c2 = mp.diff(lambda Y: Apf(xs, Y), 0)
        Rs = mp.sqrt(xs**2 - 1)
        c1_an = -2/(xs - 1) + 6*st.beta*(xs + 1)/Rs**5
        info['Aplus_linear_coeffs'] = dict(dAp_dx=s(c1, 25), dAp_dx_analytic=s(c1_an, 25), dAp_dy=s(c2, 25),
                                           ratio_abs_c2_over_c1=s(abs(c2)/abs(c1), 25),
                                           predicted_R_s=s(Rs, 25),
                                           Am_at_ring=s(st.aux(xs, 0)['Am'], 20))
        rho0 = st.aux(xs, 0)['rho']
        # e^{2gamma} normalisation (Ricci-invisible constant): axis and infinity values
        info['e2gamma_axis_and_infinity'] = dict(
            axis_x3_yplus=s(st.aux(3, 1 - mp.mpf('1e-30'))['e2g'], 20),
            axis_x1p5_yminus=s(st.aux(mp.mpf('1.5'), -1 + mp.mpf('1e-30'))['e2g'], 20),
            far_x1e12_y0p3=s(st.aux(mp.mpf('1e12'), mp.mpf('0.3'))['e2g'], 20),
            note='-2 inside the l-sum; the other reading multiplies e^{2gamma} by exp(8 beta) (Ricci-flat too, '
                 'but conical axis, and I scales by exp(-16 beta))')
        dirs = {'y=0, x>x_s': lambda dd: (xs + dd, 0),
                'y=0, x<x_s': lambda dd: (xs - dd, 0),
                'x=x_s, y>0': lambda dd: (xs, dd),
                'x=x_s, y<0': lambda dd: (xs, -dd),
                'diag (+,+)': lambda dd: (xs + dd/mp.sqrt(2), dd/mp.sqrt(2)),
                'diag (-,+)': lambda dd: (xs - dd/mp.sqrt(2), dd/mp.sqrt(2))}
        table = {}; raw = {}
        for name, pt in dirs.items():
            rows = []; Ks = []
            for dd in ds:
                X, Y = pt(dd)
                inv = invariants(st, X, Y)
                q = st.aux(X, Y)
                Wh = (q['rho'] - rho0) + 1j*q['z']
                Ks.append(inv['K'])
                rows.append(dict(d=s(dd, 3), K=s(inv['K'], 40), Kstar=s(inv['Kstar'], 20), absI=s(abs(inv['I']), 20),
                                 absI_times_B3=s(abs(inv['I'])*q['B']**3, 15),
                                 I_times_W6=s(inv['I']*Wh**6, 15), I_times_conjW6=s(inv['I']*mp.conj(Wh)**6, 15),
                                 S=s(inv['S'], 12), RicSq_over_K=s(abs(inv['RicSq'])/abs(inv['K']), 5),
                                 Ric_max=s(inv['Ric_max'], 5),
                                 bianchi_over_sqrtK=s(inv['bianchi_max']/mp.sqrt(abs(inv['K'])), 5)))
                print(lab, name, rows[-1]['d'], rows[-1]['K'][:20], rows[-1]['RicSq_over_K'], flush=True)
            loc = [s((mp.log(abs(Ks[i+1])) - mp.log(abs(Ks[i])))/(mp.log(ds[i+1]) - mp.log(ds[i])), 8)
                   for i in range(len(ds) - 1)]
            table[name] = dict(rows=rows, fit_slope=s(fit_slope(ds, Ks), 8), local_slopes=loc)
            raw[name] = Ks
            print(lab, name, 'slope', table[name]['fit_slope'], loc, flush=True)
        info['approaches'] = table
        info['fit_slope_y0_both_sides_combined_8pts'] = s(fit_slope(ds + ds, raw['y=0, x>x_s'] + raw['y=0, x<x_s']), 8)
        info['fit_slope_y0_d1e-3_to_1e-5_right_left'] = [s(fit_slope(ds[1:], raw['y=0, x>x_s'][1:]), 8),
                                                        s(fit_slope(ds[1:], raw['y=0, x<x_s'][1:]), 8)]
        info['K_C0_or_F5R_d1e-3_right_left_40dig'] = [s(raw['y=0, x>x_s'][1], 40), s(raw['y=0, x<x_s'][1], 40)]
        # leading coefficient (analytic note Sec. D): Lambda = lim e^{2gamma}/(f |W|^2), closed form at the ring,
        # and the observed law I ~ 3 / (rho0^2 Lambda^2 conj(W)^6)  =>  K delta^6 -> 48 R^4/(Lambda^2 k^8 x_s^6) on y=0
        b_ = st.beta; al_ = st.al; k_ = st.k
        gexp0 = -3*b_**2/(4*Rs**6) + 4*b_*xs/Rs - 2*b_*xs/Rs**3 - 4*b_
        psi0 = -b_/(2*Rs**3)
        Lam_an = mp.exp(gexp0 - 2*psi0)*c1.real**2*Rs**2/(k_**2*xs**4*(1 - al_**2)**2)
        Lam_num = []
        for (X, Y) in [(xs + mp.mpf('1e-9'), 0), (xs, mp.mpf('1e-9')), (xs - mp.mpf('7e-10'), mp.mpf('7e-10'))]:
            q = st.aux(X, Y); Wh = (q['rho'] - rho0) + 1j*q['z']
            Lam_num.append(s(q['e2g']/(q['f']*abs(Wh)**2), 15))
        Klead = 48*Rs**4/(Lam_an**2*k_**8*xs**6)
        Ksym = (raw['y=0, x>x_s'][3] + raw['y=0, x<x_s'][3])/2*ds[3]**6
        Ksym4 = (raw['y=0, x>x_s'][2] + raw['y=0, x<x_s'][2])/2*ds[2]**6
        Kyl = raw['x=x_s, y>0'][3]*ds[3]**6
        info['leading_coefficient'] = dict(
            Lambda_closed_form=s(Lam_an, 20), Lambda_numeric_3dirs_d1em9=Lam_num,
            N0_predicted_3_over_rho0sq_Lambdasq=s(3/(rho0**2*Lam_an**2), 15),
            N0_numeric_ItimesconjW6_d1em5_y0=table['y=0, x>x_s']['rows'][3]['I_times_conjW6'],
            Kdelta6_predicted_y0=s(Klead, 15),
            Kdelta6_numeric_sym_avg_d1em5=s(Ksym, 15), Kdelta6_numeric_sym_avg_d1em4=s(Ksym4, 15),
            Kd6_along_y_predicted=s(-Klead/Rs**6, 15), Kd6_along_y_numeric_d1em5=s(Kyl, 15))
        # precision check: repeat at dps 150 and compare raw (unrounded) values
        old = mp.mp.dps
        mp.mp.dps = 150
        xs150 = mp.findroot(B1, xs)
        K150 = invariants(st, xs150 + mp.mpf('1e-5'), 0)['K']
        K150m = invariants(st, xs150 + mp.mpf('1e-3'), 0)['K']
        Kr = raw['y=0, x>x_s'][3]; Kr3 = raw['y=0, x>x_s'][1]
        dig5 = digits(Kr, K150); dig3 = digits(Kr3, K150m)
        mp.mp.dps = old
        info['precision_check'] = dict(xs_dps150_minus_xs_dps100=s(xs150 - xs, 5),
                                       K_d1em5_dps150=s(K150, 40), digits_vs_dps100=s(dig5, 5),
                                       K_d1em3_dps150=s(K150m, 40), digits_d1em3=s(dig3, 5))
        # proper distance along y=0 from ring: s(d) = int sqrt(h_xx) dx
        prop = []
        for dd in ds:
            sd = mp.quad(lambda X: st.Hfun(X, mp.mpf(0)), [xs, xs + dd])
            prop.append(dict(d=s(dd, 3), proper=s(sd, 12)))
        info['proper_distance_y0_right'] = prop
        info['proper_slope_logs_vs_logd'] = s(fit_slope(ds, [mp.mpf(p['proper']) for p in prop]), 8)
        rings[lab] = info
    res['ring'] = rings
    res['third_point_law_check'] = third_point_check()
    print('third point', res['third_point_law_check'], flush=True)
    res['local_model_ppwave_check'] = local_model_check()
    print('local model', res['local_model_ppwave_check'], flush=True)
    res['runtime_s'] = time.time() - T0
    with open(OUT, 'w') as fh:
        json.dump(res, fh, indent=1)
    print('wrote', OUT, 'in', res['runtime_s'], 's')


def local_model_check():
    """Leading-order near-ring model (analytic note Sec. E): ds^2 = -Re(c/zeta) dT^2 - 2 rho0 dT dphi
    + Lam |zeta|^2 (dxi^2 + deta^2), zeta = xi + i eta.  Claimed: vacuum pp-wave, all invariants zero."""
    xi, eta = sp.symbols('xi eta', real=True)
    c, r0, Lam = sp.Rational(-23, 1000), sp.Rational(72, 100), sp.Rational(697, 100)
    f = c*xi/(xi**2 + eta**2)
    base = [-f, -r0 + 0*xi, sp.Integer(0) + 0*xi, Lam*(xi**2 + eta**2), Lam*(xi**2 + eta**2)]
    V = [xi, eta]
    exprs = list(base)
    d1 = [[sp.diff(e, v) for v in V] for e in base]
    for n in range(5):
        exprs += d1[n]
    for n in range(5):
        exprs += [sp.diff(d1[n][0], xi), sp.diff(d1[n][0], eta), sp.diff(d1[n][1], eta)]
    fun = sp.lambdify((xi, eta), exprs, modules='mpmath')
    out = []
    for (X, Y) in [(mp.mpf('0.3'), mp.mpf('0.1')), (mp.mpf('-0.02'), mp.mpf('0.05')), (mp.mpf('1e-4'), mp.mpf('-3e-4'))]:
        v = fun(X, Y)
        G = [v[0], v[1], v[2]]; h = [v[3], v[4]]
        d = [[v[5 + 2*n], v[6 + 2*n]] for n in range(5)]
        dd = [[v[15 + 3*n], v[16 + 3*n], v[17 + 3*n]] for n in range(5)]
        R, g, gi, detg = riemann_2p2(G, h, d, dd)
        Rm, Ru = raise_all(R, gi)
        n4 = range(4)
        K = mp.fsum(R[a][b][c_][e]*Ru[a][b][c_][e] for a, b, c_, e in itertools.product(n4, repeat=4))
        Ric = max(abs(mp.fsum(gi[a][c_]*R[a][b][c_][e] for a in n4 for c_ in n4)) for b in n4 for e in n4)
        Rmax = max(abs(Rm[a][b][c_][e]) for a, b, c_, e in itertools.product(n4, repeat=4))
        out.append(dict(xi=s(X, 5), eta=s(Y, 5), K=s(K, 5), max_Ricci=s(Ric, 5), max_mixed_Riemann=s(Rmax, 8)))
    return out


def third_point_check():
    """Universality test of the leading-coefficient law at the prereg's (chi, q) = (0.6, 0.1024) root
    (M = 1, a = 0.6, k = 0.8, beta = q/k^3 = 0.2), x_s ~ 1.0759622 (inside R < 0.5)."""
    st = Spacetime(1, sp.Rational(3, 5), sp.Rational(1, 5), 'chi0.6_q0.1024')
    xs = mp.findroot(lambda X: st.aux(X, 0)['B1'], mp.mpf('1.0759622'))
    q0 = st.aux(xs, 0); rho0 = q0['rho']
    Rs = mp.sqrt(xs**2 - 1)
    c1 = mp.diff(lambda X: st.aux(X, 0)['Ap'], xs).real
    c2 = mp.diff(lambda Y: st.aux(xs, Y)['Ap'], 0)
    b_, al_, k_ = st.beta, st.al, st.k
    gexp0 = -3*b_**2/(4*Rs**6) + 4*b_*xs/Rs - 2*b_*xs/Rs**3 - 4*b_
    Lam = mp.exp(gexp0 + b_/Rs**3)*c1**2*Rs**2/(k_**2*xs**4*(1 - al_**2)**2)
    out = dict(x_s=s(xs, 20), R_s=s(Rs, 12), c2_over_c1=s(c2/c1, 15), Lambda_closed=s(Lam, 15),
               N0_law=s(3/(rho0**2*Lam**2), 15), rows=[])
    for dd in [mp.mpf('1e-4'), mp.mpf('1e-6')]:
        for (X, Y) in [(xs + dd, 0), (xs - dd, 0), (xs, dd), (xs + dd/mp.sqrt(2), dd/mp.sqrt(2))]:
            inv = invariants(st, X, Y); q = st.aux(X, Y)
            Wh = (q['rho'] - rho0) + 1j*q['z']
            out['rows'].append(dict(x=s(X, 12), y=s(Y, 8), K=s(inv['K'], 15), S=s(inv['S'], 8),
                                    I_times_conjW6=s(inv['I']*mp.conj(Wh)**6, 12), I_times_W6=s(inv['I']*Wh**6, 12),
                                    RicSq_over_K=s(abs(inv['RicSq'])/abs(inv['K']), 4)))
    K4 = mp.mpf(out['rows'][0]['K']); K6 = mp.mpf(out['rows'][4]['K'])
    out['slope_y0_right_1e-4_to_1e-6'] = s(mp.log(abs(K6/K4))/mp.log(mp.mpf('1e-2')), 8)
    return out


def compare_v13():
    """Run ONLY after main() has written our numbers: compare with the V13 K table (read-only)."""
    V13 = _os.path.join(_HERE, '..', '..', 'V13_mn_chaotic_sea', 'results', 'v13_kretschmann.json')
    mp.mp.dps = 50
    ours = json.load(open(OUT)); v13 = json.load(open(V13))
    out = dict(note='our numbers were computed and written before V13 file was opened', rows=[])
    for lab, vlab in [('C0', 'C0'), ('F5-R', 'F5R')]:
        tab = ours['ring'][lab]['approaches']
        for r in v13[vlab]['rows']:
            side = 'y=0, x>x_s' if r['side'] == 1 else 'y=0, x<x_s'
            mine = [row['K'] for row in tab[side]['rows'] if mp.mpf(row['d']) == mp.mpf(r['d'])][0]
            rel = abs(mp.mpf(mine) - mp.mpf(r['K']))/abs(mp.mpf(mine))
            out['rows'].append(dict(point=lab, side=r['side'], d=r['d'], K_v15=mine[:22], K_v13=r['K'], rel=s(rel, 4)))
        out[lab + '_xs_v13'] = v13[vlab]['x_s']
        out[lab + '_xs_v15_minus_v13'] = s(mp.mpf(ours['ring'][lab]['x_s']) - mp.mpf(v13[vlab]['x_s']), 4)
    c0 = [r for r in out['rows'] if r['point'] == 'C0' and r['d'] == '1e-3']
    out['criterion_C0_d1e-3_rel_le_1e-6'] = all(mp.mpf(r['rel']) <= mp.mpf('1e-6') for r in c0)
    ours['comparison_v13'] = out
    json.dump(ours, open(OUT, 'w'), indent=1)
    for r in out['rows']:
        print(r)
    print('criterion C0 d=1e-3 <= 1e-6:', out['criterion_C0_d1e-3_rel_le_1e-6'])


if __name__ == '__main__':
    if '--compare-v13' in sys.argv:
        compare_v13()
    else:
        main()
