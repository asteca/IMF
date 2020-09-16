
import warnings
import numpy as np
from scipy.optimize import differential_evolution as DE
from astropy.stats import histogram as ashist
import pymc3 as pm
from . import linear_bayes_fit


def massResample(Nruns, min_mass, max_mass, mass_mean, mass_std):
    """
    """
    y_vals, x_edges = [], []
    for _ in range(Nruns):
        # Re-sample mass values
        mass_sample = np.random.normal(mass_mean, mass_std)
        # Obtain histogram
        yy, xx = ashist(mass_sample, bins=15, density=True)
        xx = .5 * (xx[1:] + xx[:-1])
        y_vals += list(yy)
        x_edges += list(xx)

    x_edges, y_vals = np.array(x_edges), np.array(y_vals)

    # Filter stars outside of given mass range(min_mass, max_mass)
    msk1 = (x_edges >= min_mass) & (x_edges <= max_mass)
    x_edges0, y_vals0 = x_edges[~msk1], y_vals[~msk1]
    x_edges, y_vals = x_edges[msk1], y_vals[msk1]

    return x_edges0, y_vals0, x_edges, y_vals


def maxLkl(mass_mean, mass_std, Nruns):
    """
    Method defined in Khalaj & Baumgardt (2013):
    https://academic.oup.com/mnras/article/434/4/3236/960889
    and used in Sheikhi et al. (2016):
    https://academic.oup.com/mnras/article/457/1/1028/989829
    """
    min_mass, max_mass = .01, 100

    def minfunc(alpha, x, xmin, xminmax, N):
        y = abs(alpha - (1 + N / (
                np.sum(np.log(x / xmin)) -
                N * (np.log(xminmax) / (1 - xminmax**(alpha - 1))))))
        return y

    # Estimate alpha using DE algorithm.
    bounds = [(1.5, 3.5)]
    alpha_lst = []
    for _ in range(Nruns):
        if _ % 10 == 0:
            print(_)
        # Re-sample mass values
        mass_sample = np.random.normal(mass_mean, mass_std)
        args = (mass_sample, mass_sample.min(),
                mass_sample.max() / mass_sample.min(), mass_sample.size)
        result = DE(minfunc, bounds, maxiter=5000, args=args)
        alpha_lst.append(result.x[0])

    alpha = np.mean(alpha_lst)
    _16p, _84p = np.percentile(alpha_lst, (16, 84))
    print(alpha, _16p, _84p)
    intercept = (1 - alpha) / (max_mass**(1 - alpha) - min_mass**(1 - alpha))

    import matplotlib.pyplot as plt

    # bins = np.linspace(min_mass, max_mass, 5000)
    yy, xx = ashist(mass_mean, bins='knuth', density=True)
    xx = .5 * (xx[1:] + xx[:-1])
    y_vals_log = 10**intercept * xx**(-alpha)
    plt.plot(xx, y_vals_log)

    vv = np.linspace(.01, 10., 1000)
    v_mean = []
    for h in vv:
        v_mean.append(np.mean((h * yy - y_vals_log)**2))
    idx = np.argmin(v_mean)
    print(vv[idx], np.mean((vv[idx] * yy - y_vals_log)**2))
    plt.scatter(xx, vv[idx] * yy, alpha=.5)

    plt.loglog()
    plt.show()

    import pdb; pdb.set_trace()  # breakpoint afdf1aaa //

    return alpha, intercept


def emceeRun(x_edges, y_vals, nsteps, nwalkers):
    """
    """
    # Obtain IMF factor.
    # Mask -inf values that can appear when a bin contains 0 elements.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_log = np.log10(y_vals)
    msk = y_log == -np.inf
    x_log = np.log10(x_edges[~msk])
    y_log = y_log[~msk]
    # No uncertainties
    xe, ye = np.zeros(x_log.size), np.zeros(x_log.size)
    data = np.array([x_log, y_log, xe, ye]).T

    # Bayesian linear fit
    trace = linear_bayes_fit.fit_data(data, nsteps, nwalkers)

    return trace


def pyMC3Run(x_edges, y_vals, nsteps):
    """
    """
    # Mask -inf values that can appear when a bin contains 0 elements.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_log = np.log10(y_vals)
    msk = y_log == -np.inf
    x_log = np.log10(x_edges[~msk])
    y_log = y_log[~msk]

    data = dict(x=x_log, y=y_log)

    with pm.Model() as model_robust:
        family = pm.glm.families.StudentT()
        pm.glm.GLM.from_formula('y ~ x', data, family=family)
        trace_robust = pm.sample(nsteps)

    return trace_robust


def LSFRun(min_mass, max_mass, mass_mean):
    """
    Least Square Fitting
    Sigma interval for polyfit:
    https://stackoverflow.com/a/28528966/1391441
    """
    yy, xx = ashist(mass_mean, bins=15, density=True)  # 'knuth'
    x_edges = .5 * (xx[1:] + xx[:-1])

    # Split stars range
    msk = (x_edges >= min_mass) & (x_edges <= max_mass)
    if (~msk).sum() > 0:
        x_edges0, y_vals0 = x_edges[~msk], yy[~msk]
    else:
        x_edges0, y_vals0 = [], []
    x_edges, y_vals = x_edges[msk], yy[msk]

    # Mask -inf values that can appear when a bin contains 0 elements.
    y_log = np.log10(y_vals)
    msk = y_log == -np.inf
    x_log = np.log10(x_edges[~msk])
    y_log = y_log[~msk]
    # Least squares fit
    pars, cov_matrix = np.polyfit(x_log, y_log, 1, cov=True)
    alpha, b = pars
    alpha_std = cov_matrix[0][0]

    return x_edges0, y_vals0, x_edges, y_vals, (alpha, b, alpha_std)
