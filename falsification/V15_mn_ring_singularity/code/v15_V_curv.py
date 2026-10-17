"""V15 verifier: plain Christoffel -> Riemann curvature engine in coordinates (t, x, y, phi).

Independent of Q2's 2+2 block formulation: the metric components and their first and second x, y derivatives are
taken exactly by sympy (from v15_V_mn.build), evaluated in mpmath, and the full 4D Christoffel symbols, Riemann
tensor, Ricci tensor, Kretschmann K = R_abcd R^abcd and Chern-Pontryagin P = *R_abcd R^abcd are formed by brute-force
index sums.  A second derivative route (mpmath numerical differentiation, mp.diff) is available as a cross-check.
"""
import itertools
import sympy as sp
import mpmath as mp
from v15_V_mn import x, y, al, be, kk, build

IDX = {'tt': (0, 0), 'tp': (0, 3), 'pp': (3, 3), 'xx': (1, 1), 'yy': (2, 2)}
DERS = [(), (x,), (y,), (x, x), (x, y), (y, y)]


class Metric:
    def __init__(self, **build_kw):
        F = build(**build_kw)
        self.F = F
        self.funcs = {}
        for name, expr in F['g'].items():
            for d in DERS:
                e = sp.diff(expr, *d) if d else expr
                self.funcs[(name, d)] = sp.lambdify((x, y, al, be, kk), e, modules='mpmath', cse=True)
        self.scalar = {n: sp.lambdify((x, y, al, be, kk), F[n], modules='mpmath', cse=True)
                       for n in ('f', 'omega', 'e2g', 'A', 'B', 'C', 'a', 'b', 'psi', 'gp')}

    def jets(self, X, Y, P):
        """Return g[a][b], dg[c][a][b], ddg[c][d][a][b] at (X, Y) as mpmath matrices/lists."""
        args = (X, Y, P['alpha'], P['beta'], P['k'])
        g = [[mp.mpf(0)] * 4 for _ in range(4)]
        dg = [[[mp.mpf(0)] * 4 for _ in range(4)] for _ in range(4)]
        ddg = [[[[mp.mpf(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
        cidx = {x: 1, y: 2}
        for (name, d), fn in self.funcs.items():
            v = fn(*args)
            if isinstance(v, mp.mpc):
                assert abs(v.imag) <= mp.mpf(10)**(-mp.mp.dps + 5) * (1 + abs(v.real)), v
                v = v.real
            v = mp.mpf(v)
            i, j = IDX[name]
            if len(d) == 0:
                g[i][j] = g[j][i] = v
            elif len(d) == 1:
                c = cidx[d[0]]
                dg[c][i][j] = dg[c][j][i] = v
            else:
                c, e = cidx[d[0]], cidx[d[1]]
                for (p, q) in ((c, e), (e, c)):
                    ddg[p][q][i][j] = ddg[p][q][j][i] = v
        return g, dg, ddg


def curvature(g, dg, ddg):
    n = 4
    G = mp.matrix(g)
    Gi = G ** -1
    gi = [[Gi[i, j] for j in range(n)] for i in range(n)]
    # d gi^{ab}/dx^c = - gi^{ae} dg_{ef,c} gi^{fb}
    dgi = [[[-sum(gi[a][e] * dg[c][e][f] * gi[f][b] for e in range(n) for f in range(n))
             for b in range(n)] for a in range(n)] for c in range(n)]
    # Gamma^a_{bc}
    low = [[[(dg[b][d][c] + dg[c][d][b] - dg[d][b][c]) / 2 for c in range(n)] for b in range(n)] for d in range(n)]
    Gam = [[[sum(gi[a][d] * low[d][b][c] for d in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
    # d_e Gamma^a_{bc}
    dlow = [[[[(ddg[e][b][d][c] + ddg[e][c][d][b] - ddg[e][d][b][c]) / 2 for c in range(n)] for b in range(n)]
             for d in range(n)] for e in range(n)]
    dGam = [[[[sum(dgi[e][a][d] * low[d][b][c] + gi[a][d] * dlow[e][d][b][c] for d in range(n))
               for c in range(n)] for b in range(n)] for a in range(n)] for e in range(n)]
    # R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}
    Rm = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        if c >= d:
            continue
        v = (dGam[c][a][d][b] - dGam[d][a][c][b]
             + sum(Gam[a][c][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][c][b] for e in range(n)))
        Rm[(a, b, c, d)] = v
        Rm[(a, b, d, c)] = -v
    for a, b, c in itertools.product(range(n), repeat=3):
        Rm[(a, b, c, c)] = mp.mpf(0)
    Ric = [[sum(Rm[(a, b, a, d)] for a in range(n)) for d in range(n)] for b in range(n)]
    Ricmix = [[sum(gi[a][c] * Ric[c][b] for c in range(n)) for b in range(n)] for a in range(n)]
    # fully lowered and fully raised
    Rl = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        Rl[(a, b, c, d)] = sum(g[a][e] * Rm[(e, b, c, d)] for e in range(n))
    # raise: R^{abcd} = gi.. ; metric is block (t,phi) + diagonal x,y, so do it generally but sparsely
    def raise_all(T):
        out = {}
        # raise index by index
        cur = T
        for pos in range(4):
            new = {}
            for key in itertools.product(range(n), repeat=4):
                s = mp.mpf(0)
                for e in range(n):
                    if gi[key[pos]][e] == 0:
                        continue
                    k2 = list(key); k2[pos] = e
                    s += gi[key[pos]][e] * cur[tuple(k2)]
                new[key] = s
            cur = new
        return cur
    Ru = raise_all(Rl)
    K = sum(Rl[k] * Ru[k] for k in Rl)
    sq = mp.sqrt(-mp.det(G))
    perm_sign = {}
    for p in itertools.permutations(range(4)):
        s = 1
        pl = list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if pl[i] > pl[j]:
                    s = -s
        perm_sign[p] = s
    # P = (1/2) eps_{abef} R^{ef}_{cd} R^{abcd};  R^{ef}_{cd} = g_{cg} g_{dh} R^{efgh}
    Rud = {}
    for e, f, c, d in itertools.product(range(n), repeat=4):
        Rud[(e, f, c, d)] = sum(g[c][gg] * g[d][hh] * Ru[(e, f, gg, hh)]
                                for gg in range(n) for hh in range(n) if g[c][gg] != 0 and g[d][hh] != 0)
    Pinv = mp.mpf(0)
    for (a, b, e, f), s in perm_sign.items():
        for c, d in itertools.product(range(n), repeat=2):
            Pinv += s * sq * Rud[(e, f, c, d)] * Ru[(a, b, c, d)]
    Pinv /= 2
    RicSq = sum(Ricmix[a][b] * Ricmix[b][a] for a in range(n) for b in range(n))
    maxRic = max(abs(Ricmix[a][b]) for a in range(n) for b in range(n))
    return dict(K=K, P=Pinv, RicSq=RicSq, maxRicmix=maxRic, Ricmix=Ricmix)
