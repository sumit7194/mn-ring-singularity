"""V15 verifier tests 6-7 (and the curvature half of test 4): is the C0 ring a genuine, reachable curvature
singularity in the verifier's own transcription?

6a  K = R_abcd R^abcd and the Chern-Pontryagin invariant P = *R_abcd R^abcd at d = 1e-2 ... 1e-5 on four rays,
    plain Christoffel->Riemann, 60 digits (and a 90-digit repeat at d = 1e-4); for BOTH gamma readings.
6b  Uniform bound: on a ring of 36 angles at each d, min and max of |I| d^6, where |I| = sqrt(K^2 + P^2)
    (= 48 |Psi_2|^2 for type D, so it cannot vanish on K's nodal directions).
6c  Proper distance to the ring along y = 0 from both sides; equatorial timelike and null geodesics from x = 1e3
    reaching the ring at finite proper time / affine parameter.
7   Ring outside the horizon rod (x = 1); unique equatorial zero of B for x > 1; reflection symmetry (so y = 0 is
    totally geodesic); no CTCs (g_phiphi > 0) on y = 0, x > x_s; Q1 closed form at the verifier's x_s.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))   # snapshot: paths relative to this file
import json, platform, time
import sympy as sp
import mpmath as mp
from v15_V_mn import build, x, y, al, be, kk, C0, F5R
from v15_V_curv import Metric, curvature

RES = _os.path.join(_HERE, '..', 'results', '')
out = dict(versions=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__))
mp.mp.dps = 110
PC = C0(mp)
Fo = build('outside')
lam = lambda e: sp.lambdify((x, y, al, be, kk), e, 'mpmath', cse=True)
fa, fb = lam(Fo['a']), lam(Fo['b'])
Pa = (PC['alpha'], PC['beta'], PC['k'])
Gq = lambda X: X + 1 + (X - 1) * fa(X, 0, *Pa) * fb(X, 0, *Pa)
xs = mp.findroot(Gq, mp.mpf('1.93'))
out['x_s'] = mp.nstr(xs, 40)

# --- 7: closed form, uniqueness, horizon, symmetry
s_ = mp.acoth(xs); kap = mp.sqrt(1 - mp.mpf('0.81'))
out['Q1_closed_form_residual'] = mp.nstr(PC['beta'] * (mp.cosh(s_) - 1)**2 * (mp.cosh(s_) + 2) - (s_ + mp.atanh(kap)), 5)
mp.mp.dps = 40
grid = [1 + mp.mpf(10)**(-6 + 12 * i / 4000.0) for i in range(4001)]  # x - 1 from 1e-6 to 1e6, log grid
signs = [mp.sign(Gq(X)) for X in grid]
changes = [mp.nstr(grid[i], 12) for i in range(len(grid) - 1) if signs[i] != signs[i + 1]]
out['equatorial_sign_changes_of_x+1+(x-1)ab_for_x_in_(1+1e-6,1+1e6)'] = changes
apb = [abs(fa(X, 0, *Pa) + fb(X, 0, *Pa)) for X in grid[::200]]
out['max|a+b| on y=0 (sampled)'] = mp.nstr(max(apb), 5)
out['ring_outside_horizon_rod'] = dict(x_s_minus_1=mp.nstr(xs - 1, 15),
                                       rho_s=mp.nstr(PC['k'] * mp.sqrt(xs**2 - 1), 15), z_s=0)
# reflection symmetry and d_y g = 0 on y = 0
g = Fo['g']
sym = {}
for nm, e in g.items():
    fe = lam(e); fey = lam(sp.diff(e, y))
    X0, Y0 = mp.mpf('2.7'), mp.mpf('0.37')
    sym[nm] = dict(rel_g_minus_g_reflected=mp.nstr((fe(X0, Y0, *Pa) - fe(X0, -Y0, *Pa)) / fe(X0, Y0, *Pa), 5),
                   dy_g_at_y0_x2p7=mp.nstr(fey(X0, 0, *Pa), 5))
out['reflection_symmetry'] = sym
fpp = lam(g['pp'])
gpp_vals = [fpp(xs + mp.mpf(10)**e_, 0, *Pa) for e_ in [-8 + 0.1 * i for i in range(120)]]
out['min g_phiphi on y=0, x in (x_s+1e-8, x_s+1e4)'] = mp.nstr(min(gpp_vals), 8)
print(out, flush=True)

# --- 6c proper distance and geodesics (outside reading; gamma reading also checked)
mp.mp.dps = 50
res6c = {}
for reading in ('outside', 'inside'):
    F = build(reading)
    gxx = lam(F['g']['xx']); gtt = lam(F['g']['tt']); gtp = lam(F['g']['tp']); gpp_ = lam(F['g']['pp'])
    r = {}
    r['proper_distance_y0_from_x_s+0.5'] = mp.nstr(mp.quad(lambda X: mp.sqrt(gxx(X, 0, *Pa)), [xs, xs + mp.mpf('0.25'), xs + mp.mpf('0.5')]), 20)
    r['proper_distance_y0_from_x_s-0.5'] = mp.nstr(mp.quad(lambda X: mp.sqrt(gxx(X, 0, *Pa)), [xs - mp.mpf('0.5'), xs - mp.mpf('0.25'), xs]), 20)
    # leading behaviour: sqrt(g_xx) ~ c |x - x_s|  ->  distance ~ c d^2 / 2
    r['sqrt(gxx)/d at d=1e-6,1e-9 (+side)'] = [mp.nstr(mp.sqrt(gxx(xs + dd, 0, *Pa)) / dd, 15) for dd in (mp.mpf('1e-6'), mp.mpf('1e-9'))]

    def Vrad(X, E, L, mu2):
        a_, b_, c_ = gtt(X, 0, *Pa), gtp(X, 0, *Pa), gpp_(X, 0, *Pa)
        det = a_ * c_ - b_**2
        Q = (c_ * E**2 + 2 * b_ * E * L + a_ * L**2) / det   # g^{ab} p_a p_b with p_t = -E, p_phi = L
        return -mu2 - Q            # = g_xx (dx/dtau)^2
    geo = []
    for (E, L, mu2) in ((mp.mpf('1.2'), mp.mpf(0), 1), (mp.mpf('1.05'), mp.mpf('0.5'), 1), (mp.mpf(1), mp.mpf(0), 0)):
        xsamp = [xs + mp.mpf(10)**(-10 + 13 * i / 600.0) for i in range(601)]
        Vmin = min(Vrad(X, E, L, mu2) for X in xsamp)
        tau = mp.quad(lambda X: mp.sqrt(gxx(X, 0, *Pa) / Vrad(X, E, L, mu2)),
                      [xs, xs + mp.mpf('1e-6'), xs + mp.mpf('1e-3'), xs + mp.mpf('0.1'), xs + 1, xs + 10, xs + 100, xs + 1000])
        tau_last = mp.quad(lambda X: mp.sqrt(gxx(X, 0, *Pa) / Vrad(X, E, L, mu2)), [xs, xs + mp.mpf('0.1')])
        geo.append({'E': mp.nstr(E, 5), 'L': mp.nstr(L, 5), 'mu2': mu2,
                    'min g_xx xdot^2 on path (x_s+1e-10 .. x_s+1e3)': mp.nstr(Vmin, 8),
                    'proper time (mu2=1) or affine param (mu2=0) from x_s+1e3 to ring': mp.nstr(tau, 15),
                    'same, last 0.1 in x': mp.nstr(tau_last, 15)})
    r['equatorial_geodesics'] = geo
    res6c[reading] = r
    print(reading, r, flush=True)
out['reachability'] = res6c
json.dump(out, open(RES + 'v15_V_t6_singularity.json', 'w'), indent=1)

# --- 6a/6b curvature at the ring, both readings
mp.mp.dps = 60
P60 = {k_: mp.mpf(v) for k_, v in PC.items()}
rays = {'+x (outside, ergo side)': 0, '+y': mp.pi / 2, '-x (inside, CTC side)': mp.pi, 'diag(-,+)': 3 * mp.pi / 4}
ds = [mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')]
curv = {}
for reading in ('outside', 'inside'):
    t0 = time.time(); Mt = Metric(gamma_reading=reading); print('built', reading, time.time() - t0, flush=True)
    rec = {'rays': {}, 'uniform': []}
    for rn, th in rays.items():
        vals = []
        for d in ds:
            c = curvature(*Mt.jets(xs + d * mp.cos(th), d * mp.sin(th), P60))
            vals.append((d, c['K'], c['P'], c['maxRicmix'] / mp.sqrt(abs(c['K']))))
        slopes = [mp.nstr(mp.log(abs(vals[i + 1][1] / vals[i][1])) / mp.log(vals[i + 1][0] / vals[i][0]), 6)
                  for i in range(len(ds) - 1)]
        Islopes = [mp.nstr(mp.log(mp.hypot(vals[i + 1][1], vals[i + 1][2]) / mp.hypot(vals[i][1], vals[i][2]))
                           / mp.log(vals[i + 1][0] / vals[i][0]), 6) for i in range(len(ds) - 1)]
        rec['rays'][rn] = dict(d=[mp.nstr(v[0], 2) for v in vals], K=[mp.nstr(v[1], 15) for v in vals],
                               P=[mp.nstr(v[2], 15) for v in vals], maxRic_over_sqrtK=[mp.nstr(v[3], 3) for v in vals],
                               K_local_slopes=slopes, I_local_slopes=Islopes)
        print(reading, rn, rec['rays'][rn], flush=True)
    # 90-digit repeat at d = 1e-4 on +x
    with mp.workdps(90):
        P90 = {k_: mp.mpf(v) for k_, v in C0(mp).items()}
        xs90 = mp.findroot(lambda X: X + 1 + (X - 1) * fa(X, 0, P90['alpha'], P90['beta'], P90['k']) * fb(X, 0, P90['alpha'], P90['beta'], P90['k']), xs)
        c90 = curvature(*Mt.jets(xs90 + mp.mpf('1e-4'), mp.mpf(0), P90))
        rec['K_+x_d1e-4_at_90_digits'] = mp.nstr(c90['K'], 30)
    for d in ds:
        Is, Ks = [], []
        for j in range(36):
            th = 2 * mp.pi * j / 36 + mp.mpf('0.01')
            c = curvature(*Mt.jets(xs + d * mp.cos(th), d * mp.sin(th), P60))
            Is.append(mp.hypot(c['K'], c['P']) * d**6); Ks.append(c['K'])
        nK = sum(1 for j in range(36) if mp.sign(Ks[j]) != mp.sign(Ks[(j + 1) % 36]))
        rec['uniform'].append(dict(d=mp.nstr(d, 2), min_I_d6=mp.nstr(min(Is), 10), max_I_d6=mp.nstr(max(Is), 10),
                                   K_sign_changes_around_circle=nK))
        print(reading, rec['uniform'][-1], flush=True)
    curv[reading] = rec
    out['curvature'] = curv
    json.dump(out, open(RES + 'v15_V_t6_singularity.json', 'w'), indent=1)
print('done')
