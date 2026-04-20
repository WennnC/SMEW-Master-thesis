# validation.py
"Created on 2026-03-01 by Wen Chai"

import numpy as np

# -----------------------------------------------------------------------
# base stats functions: RMSE, NRMSE, NSE
# -----------------------------------------------------------------------

def rmse(sim, obs):
    """
    Root Mean Square Error
    sim, obs: array-like
    """
    sim = np.array(sim)
    obs = np.array(obs)
    mask = ~np.isnan(obs) & ~np.isnan(sim)
    return np.sqrt(np.mean((sim[mask] - obs[mask])**2))


def nrmse(sim, obs, method='range'):
    """
    Normalized RMSE (Section 2.4.4)
    method='range' : RMSE / (max(obs) - min(obs))  
    method='mean'  : RMSE / mean(obs)
    returned as percentage (%)
    threshold: <30% could be taken into consideration (Moriasi et al., 2007)
    """
    sim = np.array(sim)
    obs = np.array(obs)
    mask = ~np.isnan(obs) & ~np.isnan(sim)
    obs_clean = obs[mask]
    
    r = rmse(sim, obs)
    
    if method == 'range':
        denom = np.max(obs_clean) - np.min(obs_clean)
    elif method == 'mean':
        denom = np.mean(obs_clean)
    
    if denom == 0:
        return np.nan
    
    return (r / denom) * 100   # [%]


def nse(sim, obs):
    """
    Nash-Sutcliffe Efficiency (Section 2.3.1,for adjusting model parameters Fd)
    NSE = 1: perfect; NSE < 0: model is worse than mean of obs; NSE > 0: model is better than mean of obs
    """
    sim = np.array(sim)
    obs = np.array(obs)
    mask = ~np.isnan(obs) & ~np.isnan(sim)
    sim, obs = sim[mask], obs[mask]
    
    numerator   = np.sum((obs - sim)**2)
    denominator = np.sum((obs - np.mean(obs))**2)
    
    if denominator == 0:
        return np.nan
    
    return 1 - numerator / denominator


# -----------------------------------------------------------------------
# Section 2.4.4: Objective validation function for multiple variables (NRMSE <30%)
# -----------------------------------------------------------------------

def validate_base(sim_dict, obs_dict, variables=None, nrmse_threshold=30.0):
    """
    
    sim_dict: {'pH': sim_pH_at_obs_times, 'Mg': sim_Mg, ...}
    obs_dict: {'pH': obs_pH,              'Mg': obs_Mg, ...}
    variables: Variable names to validate (default: all keys in sim_dict)
    nrmse_threshold: NRMSE < threshold is considered a pass (default:
    
    return dict: {'pH': {'RMSE': ..., 'NRMSE': ..., 'pass': True/False}, ...}
    """
    if variables is None:
        variables = list(sim_dict.keys())
    
    results = {}
    print(f"\n{'Variable':<15} {'RMSE':>10} {'NRMSE (%)':>12} {'Pass (<30%)':>12}")
    print("-" * 52)
    
    for var in variables:
        if var not in obs_dict:
            continue
        r   = rmse(sim_dict[var],  obs_dict[var])
        nr  = nrmse(sim_dict[var], obs_dict[var])
        passed = nr < nrmse_threshold if not np.isnan(nr) else False
        results[var] = {'RMSE': r, 'NRMSE': nr, 'pass': passed}
        flag = '✓' if passed else '✗'
        print(f"{var:<15} {r:>10.4f} {nr:>11.1f}% {flag:>12}")
    
    print()
    return results


# -----------------------------------------------------------------------
# Section 2.5.1: Si precipitation module evaluation (RMSE improvement >20%)
# -----------------------------------------------------------------------

def validate_Si_precip(sim_Si_base, sim_Si_precip, obs_Si,
                       improvement_threshold=20.0):
    """
    Compare RMSE of base model (no Si precipitation) vs. model with Si precipitation.
    
    sim_Si_base  : Silicate concentration without Si precipitation
    sim_Si_precip: Silicate concentration with Si precipitation
    obs_Si       : Observed silicate concentration
    
    return dict with RMSE_base, RMSE_precip, improvement(%), substantial(bool)
    """
    rmse_base   = rmse(sim_Si_base,   obs_Si)
    rmse_precip = rmse(sim_Si_precip, obs_Si)
    
    if rmse_base == 0:
        improvement = np.nan
    else:
        improvement = (rmse_base - rmse_precip) / rmse_base * 100  # [%]
    
    substantial = improvement > improvement_threshold if not np.isnan(improvement) else False
    
    print("\n--- Si Precipitation Module Evaluation (Section 2.5.1) ---")
    print(f"  RMSE (base model)      : {rmse_base:.4f}")
    print(f"  RMSE (with Si precip)  : {rmse_precip:.4f}")
    print(f"  Improvement            : {improvement:.1f}%")
    print(f"  Substantial (>20%)     : {'Yes ✓' if substantial else 'No ✗'}")
    print()
    
    return {
        'RMSE_base':    rmse_base,
        'RMSE_precip':  rmse_precip,
        'improvement':  improvement,
        'substantial':  substantial
    }


# -----------------------------------------------------------------------
# Helper function: Interpolate simulation results to observation times
# -----------------------------------------------------------------------

def interp_to_obs(t_sim, sim_values, t_obs):
    """
    Interpolate daily simulation results to irregular observation times.
    t_sim, t_obs: days (float or int), relative to start of experiment
    """
    return np.interp(t_obs, t_sim, sim_values)