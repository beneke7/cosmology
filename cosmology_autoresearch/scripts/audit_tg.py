"""Equation-level audit of Pszota & Van (2024), not an observational fit.
Dependencies: numpy, scipy. Run: python audit.py
Conventions: S(r) = 4*pi*G*rho(r), g = -d(phi)/dr, constant K > 0.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp


def g_coth(r, K=1.0, S=1.0):
    q = np.sqrt(K*S)
    return 1/(K*r) - q/(K*np.tanh(q*r))


def g_cot(r, K=1.0, S=1.0):
    q = np.sqrt(K*S)
    return -1/(K*r) + q/(K*np.tan(q*r))


def derivative_5point(f, r, h=2e-4):
    return (f(r-2*h)-8*f(r-h)+8*f(r+h)-f(r+2*h))/(12*h)


def equation_residual(f, r, K=1.0, S=1.0):
    g = f(r)
    return derivative_5point(f,r) + 2*g/r + S + K*g*g


def run_audit():
    grid = np.linspace(0.2,1.2,201)
    paper_residual = equation_residual(g_coth,grid)
    cot_residual = equation_residual(g_cot,grid)
    algebra_residual = paper_residual - 2*g_coth(grid)**2

    # Independent integrations of nonlinear g and transformed linear u.
    K = 0.4
    r0, rmax = 1e-4, 5.0
    source = lambda r: np.exp(-r*r)
    g0 = -r0/3 + (1/5-K/45)*r0**3
    u0 = 1-K*r0*r0/6+(K/20+K*K/120)*r0**4
    p0 = -K*r0/3+(K/5+K*K/30)*r0**3
    nonlinear = solve_ivp(
        lambda r,y: [-2*y[0]/r-source(r)-K*y[0]*y[0]],
        (r0,rmax),[g0],rtol=2e-11,atol=2e-13,dense_output=True)
    linear = solve_ivp(
        lambda r,y: [y[1],-2*y[1]/r-K*source(r)*y[0]],
        (r0,rmax),[u0,p0],rtol=2e-11,atol=2e-13,dense_output=True)
    if not (nonlinear.success and linear.success):
        raise RuntimeError('Integration failed')
    radii = np.linspace(0.01,rmax,1000)
    u,p = linear.sol(radii)
    direct = nonlinear.sol(radii)[0]
    transformed = p/(K*u)

    # Isolated constant-density sphere, qR=x, u -> 1 at infinity.
    # Interior u = sinc(qr)/cos(x); requires 0<x<pi/2 for
    # the everywhere-positive, regular branch connected to Newtonian gravity.
    # Exterior mass amplification is analytic, not a galaxy fit.
    xs = np.array([0.2,0.5,1.0,1.3,1.5])
    amp = 3*(np.tan(xs)/xs-1)/xs**2
    result = {
        'scope':'Equation checks and synthetic density only; no SPARC/cosmology fit; no Lean proof',
        'uniform_density_test': {
            'K':1.0,'S':1.0,'radius_range':[0.2,1.2],
            'paper_Eq23_max_absolute_residual':float(np.max(np.abs(paper_residual))),
            'corrected_cot_max_absolute_residual':float(np.max(np.abs(cot_residual))),
            'max_error_in_identity_R_paper_equals_2K_g_squared':float(np.max(np.abs(algebra_residual)))
        },
        'gaussian_density_crosscheck': {
            'K':K,'S':'exp(-r^2)','radius_range':[r0,rmax],
            'minimum_sampled_u':float(np.min(u)),
            'max_absolute_g_difference':float(np.max(np.abs(direct-transformed))),
            'max_relative_g_difference':float(np.max(np.abs((direct-transformed)/direct)))
        },
        'isolated_uniform_sphere': {
            'assumptions':'constant K>0, S>0; density truncates at R; regular centre; phi(infinity)=0; Newtonian-connected positive branch',
            'x_critical':float(np.pi/2),
            'x_qR':xs.tolist(),'M_app_over_M_b':amp.tolist()
        },
        'paper_K_flat_regime_speed_km_s':float((3.4e-5)**(-0.5))
    }
    assert np.max(np.abs(cot_residual)) < 1e-7
    assert np.max(np.abs(algebra_residual)) < 1e-7
    assert np.max(np.abs(direct-transformed)) < 1e-8
    return result


if __name__ == '__main__':
    print(json.dumps(run_audit(),indent=2))
