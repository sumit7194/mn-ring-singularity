"""V15 verifier test 3 (vacuum) and test 4 (gamma readings / axis regularity), with the verifier's own transcription.

3a  Full 4D Ricci (plain Christoffel->Riemann in (t,x,y,phi)) at 10 generic points + 4 points at d = 1e-3 from the
    C0 ring + 3 F5-R points, at 60 and 100 digits (Ricci should fall with precision if it is round-off).
3b  The same for the variant transcriptions (GLM-literal '-2 inside', Brink's 9/4 beta^2 coefficient, Brink's
    P2/r^2 in b) -- a test with teeth must reject the non-vacuum variants.
3c  Ernst equation for Brink's E = e^{2psi}A_-/A_+ (exact sympy derivatives), Re E = GLM f, and the twist relation
    between Im E and GLM's omega (Weyl rho, z chain rule).
4   e^{2gamma} on the axis y = +-1 and at spatial infinity for both readings; elementary flatness ratio
    (circumference / (2 pi proper radius)) for small circles about the axis.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))   # snapshot: paths relative to this file
import sys, json, random, platform, time
import sympy as sp
import mpmath as mp
from v15_V_mn import build, x, y, al, be, kk, C0, F5R
from v15_V_curv import Metric, curvature

RES = _os.path.join(_HERE, '..', 'results', '')
out = dict(versions=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__))


def ring_x(F, P, guess):
    """Own root of x+1+(x-1)ab = 0 on y = 0 (B = (that)^2 + (a+b)^2 there)."""
    fa = sp.lambdify((x, y, al, be, kk), F['a'], 'mpmath')
    fb = sp.lambdify((x, y, al, be, kk), F['b'], 'mpmath')
    G = lambda X: X + 1 + (X - 1) * fa(X, 0, P['alpha'], P['beta'], P['k']) * fb(X, 0, P['alpha'], P['beta'], P['k'])
    xs = mp.findroot(G, mp.mpf(guess))
    apb = fa(xs, 0, P['alpha'], P['beta'], P['k']) + fb(xs, 0, P['alpha'], P['beta'], P['k'])
    return xs, apb


def run_points(Mt, pts, dps_list=(60, 100)):
    rows = []
    for (label, PP, X, Y) in pts:
        row = dict(label=label, x=mp.nstr(X, 25), y=mp.nstr(Y, 25))
        for dps in dps_list:
            with mp.workdps(dps):
                P = {k_: mp.mpf(v) for k_, v in PP.items()}
                c = curvature(*Mt.jets(mp.mpf(X), mp.mpf(Y), P))
                row[f'K@{dps}'] = mp.nstr(c['K'], 20)
                row[f'P@{dps}'] = mp.nstr(c['P'], 20)
                row[f'maxRicmix/sqrt|K|@{dps}'] = mp.nstr(c['maxRicmix'] / mp.sqrt(abs(c['K'])), 5)
                row[f'RicSq/K@{dps}'] = mp.nstr(abs(c['RicSq'] / c['K']), 5)
        rows.append(row)
        print(label, row['x'][:12], row['y'][:12], row['K@60'][:14],
              row['maxRicmix/sqrt|K|@60'], row.get('maxRicmix/sqrt|K|@100'), flush=True)
    return rows


mp.mp.dps = 110
PC, PF = C0(mp), F5R(mp)
Fo = build('outside')
xs_C0, apb_C0 = ring_x(Fo, PC, '1.93')
xs_F5, apb_F5 = ring_x(Fo, PF, '1.92')
out['ring'] = dict(C0=dict(x_s=mp.nstr(xs_C0, 40), a_plus_b=mp.nstr(apb_C0, 5),
                           dev_from_registered=mp.nstr(xs_C0 - mp.mpf('1.92982825852'), 5)),
                   F5R=dict(x_s=mp.nstr(xs_F5, 40), a_plus_b=mp.nstr(apb_F5, 5),
                            dev_from_registered=mp.nstr(xs_F5 - mp.mpf('1.92192258792'), 5)))
print(out['ring'], flush=True)

random.seed(15)
pts = []
for i in range(10):
    pts.append((f'C0 generic {i}', PC, mp.mpf(1) + mp.mpf(random.uniform(0.05, 5.0)), mp.mpf(random.uniform(-0.95, 0.95))))
d = mp.mpf('1e-3')
for lab, (dx, dy) in {'+x': (1, 0), '-x': (-1, 0), '+y': (0, 1), 'diag(-,+)': (-1 / mp.sqrt(2), 1 / mp.sqrt(2))}.items():
    pts.append((f'C0 ring d=1e-3 {lab}', PC, xs_C0 + d * dx, d * dy))
for (X, Y) in ((mp.mpf(3), mp.mpf('0.4')), (xs_F5 + d, mp.mpf(0)), (xs_F5, d)):
    pts.append(('F5R', PF, X, Y))

import os
ONLY_ERNST_AXIS = bool(os.environ.get('V15_ONLY_ERNST_AXIS'))   # rerun 3c and 4 only, keep the stored 3a/3b
if ONLY_ERNST_AXIS:
    prev = json.load(open(RES + 'v15_V_t3_vacuum.json'))
    out.update({k_: v for k_, v in prev.items() if k_.startswith('vacuum_')})
t0 = time.time()
Mo = Metric(gamma_reading='outside') if not ONLY_ERNST_AXIS else None
print('built outside', time.time() - t0, flush=True)
if not ONLY_ERNST_AXIS:
    out['vacuum_outside'] = run_points(Mo, pts)
json.dump(out, open(RES + 'v15_V_t3_vacuum.json', 'w'), indent=1)

# 3b variants (fewer points: 3 generic + 2 ring, 60 digits only)
sub = [pts[0], pts[1], pts[2], pts[10], pts[12]]
for tag, kw in {'inside_GLM_literal': dict(gamma_reading='inside'),
                'beta2_9over4_Brink': dict(beta2_coeff=sp.Rational(9, 4)),
                'b_P2_over_r2_Brink': dict(b_variant='Brink')}.items():
    if ONLY_ERNST_AXIS:
        break
    t0 = time.time()
    Mv = Metric(**kw)
    print('built', tag, time.time() - t0, flush=True)
    out['vacuum_' + tag] = run_points(Mv, sub, dps_list=(60,))
    json.dump(out, open(RES + 'v15_V_t3_vacuum.json', 'w'), indent=1)

# 3c Ernst equation, Re E = f, twist relation
E = Fo['ernst']
Ex, Ey = sp.diff(E, x), sp.diff(E, y)
lam = lambda e: sp.lambdify((x, y, al, be, kk), e, 'mpmath', cse=True)
fE, fEx, fEy = lam(E), lam(Ex), lam(Ey)
fExx, fEyy = lam(sp.diff(E, x, 2)), lam(sp.diff(E, y, 2))
ff, fom = lam(Fo['f']), lam(Fo['omega'])
fomx, fomy = lam(sp.diff(Fo['omega'], x)), lam(sp.diff(Fo['omega'], y))
ern_rows = []
mp.mp.dps = 60
for (label, PP, X, Y) in pts:
    P = (PP['alpha'], PP['beta'], PP['k'])
    e, ex, ey, exx, eyy = fE(X, Y, *P), fEx(X, Y, *P), fEy(X, Y, *P), fExx(X, Y, *P), fEyy(X, Y, *P)
    lap = (X**2 - 1) * exx + 2 * X * ex + (1 - Y**2) * eyy - 2 * Y * ey
    grad2 = (X**2 - 1) * ex**2 + (1 - Y**2) * ey**2
    res = (mp.re(e) * lap - grad2) / max(abs(grad2), abs(mp.re(e) * lap))
    fv = ff(X, Y, *P)
    # twist: rho = k sqrt((x^2-1)(1-y^2)), z = k x y; d/drho, d/dz via inverse Jacobian
    k_ = PP['k']
    rho = k_ * mp.sqrt((X**2 - 1) * (1 - Y**2))
    J = mp.matrix([[k_ * X * mp.sqrt((1 - Y**2) / (X**2 - 1)), -k_ * Y * mp.sqrt((X**2 - 1) / (1 - Y**2))],
                   [k_ * Y, k_ * X]])   # rows: rho, z ; cols: d/dx, d/dy
    JiT = (J ** -1).T                  # [d_rho; d_z] = JiT [d_x; d_y]
    Phx, Phy = mp.im(ex), mp.im(ey)
    Phr = JiT[0, 0] * Phx + JiT[0, 1] * Phy
    Phz = JiT[1, 0] * Phx + JiT[1, 1] * Phy
    omx, omy = fomx(X, Y, *P), fomy(X, Y, *P)
    omr = JiT[0, 0] * omx + JiT[0, 1] * omy
    omz = JiT[1, 0] * omx + JiT[1, 1] * omy
    den1, den2 = -(fv**2 / rho) * omz, (fv**2 / rho) * omr
    # on y = 0 both d_z Phi and d_z omega vanish by reflection symmetry: report the numerator instead
    r1 = Phr / den1 if den1 != 0 else ('den=0, num=' + mp.nstr(Phr, 5))
    r2 = Phz / den2 if den2 != 0 else ('den=0, num=' + mp.nstr(Phz, 5))
    ern_rows.append(dict(label=label, ernst_rel_residual=mp.nstr(res, 5),
                         ReE_minus_f_rel=mp.nstr((mp.re(e) - fv) / fv, 5),
                         twist_ratio_rho=r1 if isinstance(r1, str) else mp.nstr(r1, 25),
                         twist_ratio_z=r2 if isinstance(r2, str) else mp.nstr(r2, 25)))
    print(ern_rows[-1], flush=True)
out['ernst'] = ern_rows
json.dump(out, open(RES + 'v15_V_t3_vacuum.json', 'w'), indent=1)

# 4 axis regularity
xp = sp.symbols('xp', positive=True)
ax = {}
for r in ('outside', 'inside'):
    F = build(r)
    rec = {}
    for ysgn in (1, -1):
        # approach the axis: y = ysgn*(1 - eps); exact limit via substitution with x > 1 positive
        # y = +-1 is not a 0/0 point of e^{2gamma} (A ~ (x^2-1)(1+ab)^2 there), so direct substitution is exact
        e2g_axis = sp.simplify(F['e2g'].subs(y, ysgn).subs(x, xp + 1))
        rec[f'e2g(y={ysgn})'] = str(e2g_axis)
    ax[r] = rec
    print(r, rec, flush=True)
# numeric: at infinity and the elementary-flatness ratio
for r in ('outside', 'inside'):
    F = build(r)
    fe2g = lam(F['e2g'])
    gpp = lam(F['g']['pp'])
    Hyy = lam(kk**2 * F['e2g'] * (x**2 - y**2) / F['f'])   # = g_yy (1 - y^2), regular on the axis
    P = (PC['alpha'], PC['beta'], PC['k'])
    rec = ax[r]
    rec['e2g_at_x=1e6_y=0'] = mp.nstr(fe2g(mp.mpf(10)**6, mp.mpf(0), *P), 20)
    rec['e2g_at_x=1e6_y=0.7'] = mp.nstr(fe2g(mp.mpf(10)**6, mp.mpf('0.7'), *P), 20)
    for X in (mp.mpf('1.5'), mp.mpf(3), mp.mpf(10)):
        for ysgn in (1, -1):
            vals = []
            for eps in (mp.mpf('1e-6'), mp.mpf('1e-10')):
                Y1 = ysgn * (1 - eps)
                circ = 2 * mp.pi * mp.sqrt(gpp(X, Y1, *P))
                # proper radius int sqrt(g_yy) dy from the axis; with 1 - |y| = t^2 the integrand is
                # 2 sqrt(H/(2 - t^2)), H = g_yy (1 - y^2), regular at t = 0
                rad = mp.quad(lambda t: 2 * mp.sqrt(Hyy(X, ysgn * (1 - t**2), *P) / (2 - t**2)), [0, mp.sqrt(eps)])
                vals.append(mp.nstr(circ / (2 * mp.pi * rad), 20))
            rec[f'circ/(2pi radius) x={mp.nstr(X, 3)} y={ysgn}(1-eps) eps=1e-6,1e-10'] = vals
    rec['exp(-4beta) C0'] = mp.nstr(mp.exp(-4 * PC['beta']), 20)
    print(r, rec, flush=True)
out['axis'] = ax
json.dump(out, open(RES + 'v15_V_t3_vacuum.json', 'w'), indent=1)
print('done')
