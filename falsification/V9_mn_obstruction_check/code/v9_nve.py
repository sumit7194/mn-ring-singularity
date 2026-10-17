"""V9 -- the bridge's own axial NVE for Manko-Novikov (q-anomaly subclass), from ansatz's manifest formulas.
Imports nothing from quantum. Branch: R = sqrt(x^2+y^2-1) implemented as x*sqrt(1+(y^2-1)/x^2) (=x on the axis).
Axis Gamma: y = cos(theta), theta = 0, L = 0.  A = g^{theta theta}|0, B = -d_y H|_{y=1}, xdot^2 = F1 (E^2/f - mu^2).
xi1 (= d theta) equation: xi1'' - (a'/a) xi1' + a b xi1 = 0, a = A/xdot, b = B/xdot.  Normal form z'' = r z,
p = -a'/a = -(A'/A - (xdot2)'/(2 xdot2)),  q = A B / xdot2,  r = p^2/4 + p'/2 - q."""
import sympy as sp
x, y = sp.symbols('x y')
def mn_functions(M, a_, beta):
    M, a_, beta = map(sp.nsimplify, (M, a_, beta))
    chi = a_/M; al = sp.nsimplify((-1 + sp.sqrt(1 - chi**2))/chi); k = sp.sqrt(M**2 - a_**2)
    R = x*sp.sqrt(1 + (y**2 - 1)/x**2); u = x*y/R
    P = [sp.Integer(1), u, (3*u**2 - 1)/2, (5*u**3 - 3*u)/2]
    psi = beta*P[2]/R**3
    A_ = -al*sp.exp(-2*beta*(-1 + sum((x - y)*P[l]/R**(l + 1) for l in range(3))))
    B_ = al*sp.exp(2*beta*(1 + sum((-1)**(3 - l)*(x + y)*P[l]/R**(l + 1) for l in range(3))))
    Am = (x**2 - 1)*(1 + A_*B_)**2 - (1 - y**2)*(B_ - A_)**2
    Bm = (x + 1 + (x - 1)*A_*B_)**2 + ((1 + y)*A_ + (1 - y)*B_)**2
    Cm = (x**2 - 1)*(1 + A_*B_)*(B_ - A_ - y*(A_ + B_)) + (1 - y**2)*(B_ - A_)*(1 + A_*B_ + x*(1 - A_*B_))
    f = sp.exp(2*psi)*Am/Bm
    om = 2*k*sp.exp(-2*psi)*Cm/Am - 4*k*al/(1 - al**2)
    e2g = (x**2 - 1)/(x**2 - y**2)*sp.exp(3*beta**2*(P[3]**2 - P[2]**2)/R**6
            + 2*beta*sum((x - y + (-1)**l*(x + y))*P[l]/R**(l + 1) - 2 for l in range(3)) + 8*beta)*Am/((x**2 - 1)*(1 - al**2)**2)
    return dict(f=f, om=om, e2g=e2g, k=k, al=al)

def axial_r(M, a_, beta, E, mu2):
    F = mn_functions(M, a_, beta); f, om, e2g, k = F['f'], F['om'], F['e2g'], F['k']
    F1 = f*(x**2 - 1)/(k**2*e2g*(x**2 - y**2))
    at1 = lambda e: e.subs(y, 1)
    f1, e2g1, F11 = at1(f), at1(e2g), at1(F1)
    fy1, F1y1, om1 = at1(sp.diff(f, y)), at1(sp.diff(F1, y)), at1(sp.diff(om, y))
    omega_on_axis = at1(om)
    px2 = (E**2/f1 - mu2)/F11
    A = f1/(k**2*e2g1*(x**2 - 1))
    B = -sp.Rational(1, 2)*(F1y1*px2 - E**2*f1*om1**2/(2*k**2*(x**2 - 1)) + E**2*fy1/f1**2)
    xdot2 = F11*(E**2/f1 - mu2)
    p = -(sp.diff(A, x)/A - sp.diff(xdot2, x)/(2*xdot2))
    q = A*B/xdot2
    r = p**2/4 + sp.diff(p, x)/2 - q
    return dict(r=r, A=A, B=B, xdot2=xdot2, omega_on_axis=omega_on_axis, F=F)
