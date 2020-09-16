
import numpy as np
import arviz as az


def proc(mcmc_mode, burn_frac, nsteps, trace):
    """
    """
    func_dict = {
        "16p": lambda x: np.percentile(x, 16),
        "median": lambda x: np.percentile(x, 50),
        "84p": lambda x: np.percentile(x, 84),
    }

    nburn = int(burn_frac * nsteps)
    if mcmc_mode == 'pymc3':
        az_data = az.from_pymc3(trace).sel(
            draw=slice(nburn, None))
        summ = az.summary(
            az_data, var_names=["x", "Intercept"], stat_funcs=func_dict)
        flat_samples = None

    elif mcmc_mode == 'emcee':
        var_names = (
            'alpha', 'b', 'in_scatt', 'out_mean', 'out_scatt', 'out_fr')
        nburn = int(burn_frac * nsteps)
        az_data = az.from_emcee(trace, var_names=var_names).sel(
            draw=slice(nburn, None))
        flat_samples = trace.get_chain(discard=nburn, flat=True)

        summ = az.summary(
            az_data, var_names=["alpha", "b"], stat_funcs=func_dict)

    return az_data, summ, flat_samples
