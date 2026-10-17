"""V15 verifier test 2: compare the verifier's own transcription (v15_V_mn) with the fleet's v9_nve.mn_functions.

Numerical comparison at 12 random points per parameter set (C0, F5-R) at 80 digits, plus an exact symbolic check of the
gamma exponent.  Exact parameters as in the V13 POINTS tables (rationals).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))   # snapshot: paths relative to this file
import sys, json, random, platform
import sympy as sp
import mpmath as mp
sys.path.insert(0, _os.path.join(_HERE, '..', '..', 'V9_mn_obstruction_check', 'code'))
import v9_nve
from v15_V_mn import build, x, y, al, be, kk, params

mp.mp.dps = 80
Rt = sp.Rational
SETS = {'C0': (1, Rt(9, 10), Rt(19, 20) / Rt(19, 100)**Rt(3, 2)),
        'F5R': (89, 80, Rt(13394411, 1186380))}
out = dict(versions=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__),
           dps=mp.mp.dps, sets={})
mine = {r: build(gamma_reading=r) for r in ('outside', 'inside')}
random.seed(20261017)
for name, (M, a_, beta) in SETS.items():
    V = v9_nve.mn_functions(M, a_, beta)
    # parameters, my route (GLM Eq. 5)
    chi = sp.nsimplify(a_) / M
    alpha_s = (-1 + sp.sqrt(1 - chi**2)) / chi
    k_s = M * (1 - alpha_s**2) / (1 + alpha_s**2)
    beta_mine = Rt(19, 20) * M**3 / k_s**3
    sub = {al: alpha_s, be: beta, kk: k_s}
    rec = dict(param_check=dict(
        k_mine_minus_v9=str(sp.nsimplify(sp.simplify(k_s - V['k']))),
        alpha_mine_minus_v9=str(sp.simplify(alpha_s - V['al'])),
        beta_from_q_minus_registered=mp.nstr(mp.mpf(sp.N(beta_mine - beta, 90)), 5)), points=[])
    fn_m = {r: {q: sp.lambdify((x, y), mine[r][q].subs(sub), 'mpmath') for q in ('f', 'omega', 'e2g')}
            for r in mine}
    fn_v = {q: sp.lambdify((x, y), V[q2], 'mpmath') for q, q2 in (('f', 'f'), ('omega', 'om'), ('e2g', 'e2g'))}
    worst = {'outside': 0, 'inside': 0}
    for i in range(12):
        X = mp.mpf(1) + mp.mpf(random.uniform(0.02, 6.0))
        Y = mp.mpf(random.uniform(-0.98, 0.98))
        row = dict(x=mp.nstr(X, 20), y=mp.nstr(Y, 20))
        for r in ('outside', 'inside'):
            for q in ('f', 'omega', 'e2g'):
                vm = fn_m[r][q](X, Y); vv = fn_v[q](X, Y)
                rel = abs(vm - vv) / max(abs(vv), mp.mpf(10)**-70)
                row[f'{r}:{q}'] = mp.nstr(rel, 5)
                worst[r] = max(worst[r], rel)
        if i == 0:
            row['e2g_ratio_inside_over_v9'] = mp.nstr(fn_m['inside']['e2g'](X, Y) / fn_v['e2g'](X, Y), 30)
            row['exp(-8beta)'] = mp.nstr(mp.exp(-8 * mp.mpf(sp.N(beta, 90))), 30)
        rec['points'].append(row)
    rec['worst_rel_outside_reading'] = mp.nstr(worst['outside'], 5)
    rec['worst_rel_inside_reading'] = mp.nstr(worst['inside'], 5)
    # exact check of the gamma exponent: replace both spellings of R by one symbol and expand
    Rs = sp.symbols('Rs', positive=True)
    E9 = [a_ for a_ in sp.Mul.make_args(V['e2g']) if isinstance(a_, sp.exp)]
    assert len(E9) == 1, E9
    E9 = sp.Add(*[e.args[0] for e in E9])
    mine_exp = 2 * (mine['outside']['gp'] - sp.Rational(1, 2) * sp.log((x**2 - 1) / (x**2 - y**2)))
    rep = lambda e: e.subs({sp.sqrt(x**2 + y**2 - 1): Rs}).subs({x * sp.sqrt(1 + (y**2 - 1) / x**2): Rs})
    dl = sp.expand(rep(mine_exp.subs(sub)) - rep(E9))
    rec['symbolic_gamma_exponent_difference_outside'] = str(sp.nsimplify(sp.simplify(dl)))
    out['sets'][name] = rec
    print(name, rec['param_check'], 'worst outside', rec['worst_rel_outside_reading'],
          'worst inside', rec['worst_rel_inside_reading'], 'symbolic', rec['symbolic_gamma_exponent_difference_outside'][:80])
json.dump(out, open(_os.path.join(_HERE, '..', 'results', 'v15_V_t2_compare.json'), 'w'),
          indent=1)
