"""V15 adversarial verifier: an INDEPENDENT transcription of the Manko-Novikov (MN) metric.

Written from the arXiv LaTeX sources (read 2026-10-17, before opening V9's v9_nve.py):
  GLM  = Gair, Li, Mandel 2008, arXiv:0708.0628 (e-print dated 2008-01-11), Eq. 2 (line element, prolate),
         Eq. 3a-3l (metric functions), Eq. 4 (moments), Eq. 5 (k, alpha, beta), Eq. 6-7 (Weyl rho, z).
  Brink 2008, arXiv:0807.1179, Appendix C ("Manko-Novikov Metric"): Ernst potential E = e^{2 psi~} A_-/A_+,
         A_-+ = x(1+ab) + i y(b-a) -+ (1-ia)(1-ib); metric functions; a, b, psi~, gamma~.
  LGAC = Lukes-Gerakopoulos, Apostolatos, Contopoulos 2010, arXiv:1003.3120, Eqs. (ffunc)-(fLeg), (freepar).

Paper-vs-paper differences recorded (see results/v15_V_verifier.md):
  D1  gamma' linear-in-beta term: GLM 3h writes beta * sum_{l=0}^{2} ( [x-y+(-1)^{2-l}(x+y)] P_l/R^{l+1} - 2 )
      ("-2 inside", total -6 beta); LGAC writes beta(-2 + sum ...) and says it moved the -2 outside "to avoid
      any confusion"; Brink's form -1/2(ln(a/-alpha)+ln(b/alpha)) expands to the "-2 outside" reading.
  D2  beta^2 coefficient in gamma': GLM 9 beta^2/(6 R^6) = 3/2; LGAC 3 beta^2/(2 R^6) = 3/2;
      Brink alpha_2^2 (3/2)^2 = 9/4.
  D3  b: GLM 3j and LGAC (fb) have P_2/R^3; Brink's ln(b/alpha) has P_2/r^2 (an r^2 where r^3 is expected).
  D4  Brink's line element k^2 e^{-2psi}[e^{2gamma}(drho^2+dz^2)+R^2 dphi^2] with R^2 = k^2(x^2-1)(1-y^2)
      (k^4 in g_phiphi, unless k=1, which Brink's comment says the original MN paper uses); his e^{2gamma}
      carries an extra (x^2-y^2) because his "rho, z" are x = cosh rho, y = cos z (not Weyl rho, z).
All variants are implemented so that the vacuum and axis tests can adjudicate between them.
"""
import sympy as sp

x, y = sp.symbols('x y', real=True)
al, be, kk = sp.symbols('alpha beta k', real=True)
u_ = sp.symbols('u', real=True)


def legendre_rodrigues(n):
    """P_n(u) from the Rodrigues formula quoted under GLM Eq. 3l / LGAC (fLeg)."""
    return sp.expand(sp.diff((u_**2 - 1)**n, u_, n) / (2**n * sp.factorial(n)))


R = sp.sqrt(x**2 + y**2 - 1)                       # GLM 3k
U = x * y / R                                      # argument of P_n, GLM 3l
P = [legendre_rodrigues(n).subs(u_, U) for n in range(4)]


def build(gamma_reading='outside', beta2_coeff=sp.Rational(3, 2), b_variant='GLM'):
    """Return dict of sympy expressions for the MN functions.

    gamma_reading: 'outside' (LGAC/Brink: beta(-2+sum)) or 'inside' (GLM 3h literal: beta*sum(...-2)).
    beta2_coeff:  3/2 (GLM 9/6, LGAC 3/2) or 9/4 (Brink (3/2)^2).
    b_variant:    'GLM' (P_2/R^3) or 'Brink' (P_2/r^2 as printed in Brink App. C).
    """
    # GLM 3i (= LGAC fa = Brink ln(a/-alpha))
    Sa = sum((x - y) * P[l] / R**(l + 1) for l in range(3))
    a = -al * sp.exp(-2 * be * (-1 + Sa))
    # GLM 3j (= LGAC fb); Brink variant has P2/r^2
    if b_variant == 'GLM':
        Sb = sum((-1)**(3 - l) * (x + y) * P[l] / R**(l + 1) for l in range(3))
    else:
        Sb = -(x + y) * (P[0] / R - P[1] / R**2 + P[2] / R**2)
    b = al * sp.exp(2 * be * (1 + Sb))
    psi = be * P[2] / R**3                          # GLM 3g
    A = (x**2 - 1) * (1 + a * b)**2 - (1 - y**2) * (b - a)**2                       # GLM 3d
    B = (x + 1 + (x - 1) * a * b)**2 + ((1 + y) * a + (1 - y) * b)**2               # GLM 3e
    C = ((x**2 - 1) * (1 + a * b) * (b - a - y * (a + b))
         + (1 - y**2) * (b - a) * (1 + a * b + x * (1 - a * b)))                    # GLM 3f
    lin_terms = [(x - y + (-1)**(2 - l) * (x + y)) / R**(l + 1) * P[l] for l in range(3)]
    if gamma_reading == 'inside':
        lin = be * sum(t - 2 for t in lin_terms)                                    # GLM 3h literal
    else:
        lin = be * (-2 + sum(lin_terms))                                            # LGAC (fgampr)
    gp = (sp.Rational(1, 2) * sp.log((x**2 - 1) / (x**2 - y**2))
          + beta2_coeff * be**2 / R**6 * (P[3]**2 - P[2]**2) + lin)                 # GLM 3h
    f = sp.exp(2 * psi) * A / B                                                     # GLM 3a
    omega = 2 * kk * sp.exp(-2 * psi) * C / A - 4 * kk * al / (1 - al**2)           # GLM 3b
    e2g = sp.exp(2 * gp) * A / ((x**2 - 1) * (1 - al**2)**2)                        # GLM 3c
    # Brink App. C Ernst potential (independent route to f and the twist potential)
    Ap = x * (1 + a * b) + sp.I * y * (b - a) + (1 - sp.I * a) * (1 - sp.I * b)
    Am = x * (1 + a * b) + sp.I * y * (b - a) - (1 - sp.I * a) * (1 - sp.I * b)
    ernst = sp.exp(2 * psi) * Am / Ap
    # Brink's gamma~ written his way: 1/2 ln((x^2-1)/(x^2-y^2)) - 1/2(ln(a/-alpha)+ln(b/alpha)) + c beta^2 (...)
    gB = (sp.Rational(1, 2) * sp.log((x**2 - 1) / (x**2 - y**2))
          - sp.Rational(1, 2) * (2 * be * (1 - Sa) + 2 * be * (1 + Sb))
          + beta2_coeff * be**2 / R**6 * (P[3]**2 - P[2]**2))
    # GLM Eq. 2 line element in (t, x, y, phi)
    rho2 = kk**2 * (x**2 - 1) * (1 - y**2)
    g = {
        'tt': -f,
        'tp': f * omega,
        'pp': rho2 / f - f * omega**2,
        'xx': kk**2 * e2g * (x**2 - y**2) / (f * (x**2 - 1)),
        'yy': kk**2 * e2g * (x**2 - y**2) / (f * (1 - y**2)),
    }
    return dict(a=a, b=b, psi=psi, A=A, B=B, C=C, gp=gp, gB=gB, f=f, omega=omega, e2g=e2g,
                ernst=ernst, Ap=Ap, Am=Am, g=g, rho2=rho2)


def params(M, a_spin, q, mp):
    """GLM Eq. 5: alpha = (-1+sqrt(1-chi^2))/chi, k = M(1-alpha^2)/(1+alpha^2), beta = q M^3/k^3."""
    M = mp.mpf(M); a_spin = mp.mpf(a_spin); q = mp.mpf(q)
    chi = a_spin / M
    alpha = (-1 + mp.sqrt(1 - chi**2)) / chi
    k = M * (1 - alpha**2) / (1 + alpha**2)
    beta = q * M**3 / k**3
    return alpha, beta, k


def C0(mp):
    alpha, beta, k = params(1, mp.mpf('0.9'), mp.mpf('0.95'), mp)
    return dict(M=mp.mpf(1), alpha=alpha, beta=beta, k=k)


def F5R(mp):
    alpha, beta, k = params(89, 80, mp.mpf(19) / 20, mp)
    return dict(M=mp.mpf(89), alpha=alpha, beta=beta, k=k)
