"""
model_simulation.py
Canonical Simulation Script for Zenodo v1.0.1
Title: Treatment Timing, Metabolic Modulation, and Bifurcation-Proximity Control in a Coupled Tumor-Immune-Drug Model
Author: Yuji Marutani (O'Valley Salon of Meta Design)
DOI: 10.5281/zenodo.21864631
"""

import numpy as np
from scipy.integrate import solve_ivp

# ==========================================
# 1. Canonical Parameter Dictionary (v1.0.1)
# ==========================================
base_params = {
    'r': 0.8,         'K': 1000.0,
    'k_T': 0.03,      'k_NK': 0.025,
    'MHC_I_0': 1.0,   'MHC_I_tum': 0.2,
    'NKG2DL': 0.8,    'K_NK': 0.5,
    'K_E': 30.0,      
    'V_max': 1.0,     'K_m': 0.5,       'P_gp': 2.0,    'K_ATP': 20.0,
    'gamma_D': 0.2,   'alpha': 1.5,     'beta': 2.5,
    'n': 2.0,         'h': 0.8,         'epsilon': 0.1,
    'g_m': 0.5,       'eta': 1.5,       'mu_m': 0.3,
    'chi': 0.05,      'p_T': 2.0,       'd_T': 0.05,    'S_T': 0.02,
    'p_NK': 1.5,      'd_NK': 0.04,     'S_NK': 0.02,
    'I_amp': 0.45     # Fixed Canonical Amplitude
}

# ==========================================
# 2. ODE System Definition
# ==========================================
def model_ode(t, y, p, f, delta_t, B_input, drug_on=True):
    C, D, M, E_T, E_NK = y
    C, D, M = max(0.0, C), max(0.0, D), max(0.0, M)
    E_T, E_NK = max(0.0, E_T), max(0.0, E_NK)
    
    theta_T = p['k_T'] * (p['MHC_I_tum'] / (p['MHC_I_0'] + p['MHC_I_tum']))
    theta_NK = p['k_NK'] * (1.0 - p['MHC_I_tum'] / p['MHC_I_0']) * (p['NKG2DL'] / (p['K_NK'] + p['NKG2DL']))
    kill_rate = p['alpha'] * (D**p['n'] / (p['h']**p['n'] + D**p['n'])) * (1.0 + p['beta'] / (M + p['epsilon']))
    
    I_t = p['I_amp'] * (1.0 + np.cos(2.0 * np.pi * f * (t - delta_t))) if (drug_on and t >= delta_t) else 0.0
    
    dC_dt = p['r'] * C * (1.0 - C / p['K']) - (theta_T * E_T + theta_NK * E_NK) * C - kill_rate * C
    
    atp_factor = M / (p['K_ATP'] + M + 1e-6)
    efflux = (p['V_max'] * atp_factor / (p['K_m'] + D + 1e-6)) * p['P_gp'] * D
    
    dD_dt = I_t - efflux - p['gamma_D'] * D
    dM_dt = p['g_m'] * C - p['eta'] * B_input * M - p['mu_m'] * M
    
    icd_signal = p['chi'] * kill_rate * C
    dE_T_dt = p['p_T'] * icd_signal * (E_T / (p['K_E'] + E_T)) - p['d_T'] * E_T + p['S_T']
    dE_NK_dt = p['p_NK'] * icd_signal * (E_NK / (p['K_E'] + E_NK)) - p['d_NK'] * E_NK + p['S_NK']
    
    return [dC_dt, dD_dt, dM_dt, dE_T_dt, dE_NK_dt]

def run_sim(f_val, dt_val, B_val, p_dict, drug_on=True, t_max=40):
    y0 = [800.0, 0.0, 100.0, 2.0, 2.0]
    sol = solve_ivp(
        lambda t, y: model_ode(t, y, p_dict, f_val, dt_val, B_val, drug_on=drug_on),
        (0, t_max), y0, method='BDF', t_eval=np.linspace(0, t_max, 1000)
    )
    return sol

if __name__ == "__main__":
    sol_untreated = run_sim(0.0, 0.0, 0.0, base_params, drug_on=False)
    print(f"Untreated C(40): {sol_untreated.y[0][-1]:.2f}")
