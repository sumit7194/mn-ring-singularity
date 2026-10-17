"""V15 verifier test 5: asymptotics / multipole moments from the verifier's own transcription.

(a) Exact: Brink's Ernst potential on the upper axis (rho = 0, z = k x), xi = (1 - E)/(1 + E) expanded in 1/z;
    Fodor-Hoenselaers-Perjes: xi = sum m_n z^{-n-1} and the Geroch-Hansen moments equal m_n for n <= 3.
    Compared with GLM Eq. 4 (M0, S1, M2, S3) and with Q = M2 = -M(a^2 + q M^2) (a = chi M).
(b) Numeric far field of the GLM metric functions on the equator and off it: M from f = 1 - 2M/r + ...,
    J from g_tphi = f omega ~ -/+ 2 J sin^2(theta)/r, and e^{2gamma} -> 1.
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))   # snapshot: paths relative to this file
import json, platform
import sympy as sp
import mpmath as mp
from v15_V_mn import build, x, y, al, be, kk, C0, F5R

RES = _os.path.join(_HERE, '..', 'results', '')
out = dict(versions=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__))

w = sp.symbols('w', positive=True)
xp = sp.symbols('xp', positive=True)
F = build('outside')
E_axis = F['ernst'].subs(y, 1).subs(x, xp)
E_axis = sp.simplify(E_axis.subs(sp.sqrt(xp**2), xp))
E_w = E_axis.subs(xp, 1 / (kk * w))        # z = k x  ->  x = 1/(k w), w = 1/z
xi = (1 - E_w) / (1 + E_w)
ser = sp.series(xi, w, 0, 5).removeO()
m = [sp.simplify(sp.factor(ser.coeff(w, n + 1))) for n in range(4)]
glm = dict(M0=kk * (1 + al**2) / (1 - al**2),
           S1=-2 * al * kk**2 * (1 + al**2) / (1 - al**2)**2,
           M2=-kk**3 * (be + 4 * al**2 * (1 + al**2) / (1 - al**2)**3),
           S3=4 * al * kk**4 * (be + 2 * al**2 * (1 + al**2) / (1 - al**2)**3) / (1 - al**2))
chk = dict(
    m0_minus_M0=sp.simplify(m[0] - glm['M0']),
    m1=sp.simplify(m[1]),
    m1_over_i_minus_S1=sp.simplify(m[1] / sp.I - glm['S1']),
    m1_over_i_plus_S1=sp.simplify(m[1] / sp.I + glm['S1']),
    m2_minus_M2=sp.simplify(m[2] - glm['M2']),
    m3_over_i_minus_S3=sp.simplify(m[3] / sp.I - glm['S3']),
    m3_over_i_plus_S3=sp.simplify(m[3] / sp.I + glm['S3']),
)
# in terms of M, chi, q (GLM Eq. 5): M2 = -M^3 (chi^2 + q)
Mq, chi, q = sp.symbols('M chi q', positive=True)
alpha_e = (-1 + sp.sqrt(1 - chi**2)) / chi
k_e = Mq * (1 - alpha_e**2) / (1 + alpha_e**2)
sub5 = {al: alpha_e, kk: k_e, be: q * Mq**3 / k_e**3}
chk['M0(M,chi,q)'] = sp.simplify(m[0].subs(sub5))
chk['S1(M,chi,q) = |m1|'] = sp.simplify((m[1] / sp.I).subs(sub5))
chk['M2(M,chi,q)'] = sp.simplify(m[2].subs(sub5))
chk['M2 + M^3(chi^2+q)'] = sp.simplify(m[2].subs(sub5) + Mq**3 * (chi**2 + q))
out['axis_series'] = {f'm{n}': str(m[n]) for n in range(4)}
out['checks'] = {k_: str(v) for k_, v in chk.items()}
for k_, v in out['checks'].items():
    print(k_, ':', v)

# (b) numeric far field
mp.mp.dps = 60
lam = lambda e: sp.lambdify((x, y, al, be, kk), e, 'mpmath')
ff, fom, fe2g = lam(F['f']), lam(F['omega']), lam(F['e2g'])
num = {}
for name, P in (('C0', C0(mp)), ('F5R', F5R(mp))):
    Pa = (P['alpha'], P['beta'], P['k'])
    rows = []
    for X in (mp.mpf(10)**4, mp.mpf(10)**6, mp.mpf(10)**8):
        for Y in (mp.mpf(0), mp.mpf('0.6')):
            rho = P['k'] * mp.sqrt((X**2 - 1) * (1 - Y**2)); z = P['k'] * X * Y
            r = mp.sqrt(rho**2 + z**2); s2 = rho**2 / r**2
            fv, omv = ff(X, Y, *Pa), fom(X, Y, *Pa)
            rows.append(dict(x=mp.nstr(X, 3), y=mp.nstr(Y, 3),
                             M_from_f=mp.nstr(r * (1 - fv) / 2, 15),
                             J_from_gtphi=mp.nstr(-(fv * omv) * r / (2 * s2), 15),
                             e2g_minus_1=mp.nstr(fe2g(X, Y, *Pa) - 1, 5)))
    num[name] = dict(M=mp.nstr(P['M'], 10), J_expected_aM=mp.nstr(P['M'] * P['M'] * (mp.mpf('0.9') if name == 'C0' else mp.mpf(80) / 89), 15),
                     rows=rows)
    print(name, num[name])
out['far_field_numeric'] = num
json.dump(out, open(RES + 'v15_V_t5_multipoles.json', 'w'), indent=1)
