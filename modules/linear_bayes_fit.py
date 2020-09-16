
import warnings
import emcee
from emcee import moves
import numpy as np
from scipy.optimize import differential_evolution as DE
from . import update_progress


def fit_data(
    data, nsteps, nwalkers,
        priorlimits=[(-2.8, -1.8), (-1000., 1000.), (0.001, 100.), (-10., 10.),
                     (0.001, 1000.), (0., 1.)]):
    """
    This code will fit a straight line with intrinsic dispersion to data with
    (optionally) covariant errors on both the independent and dependent
    variables, and takes care of outliers using a mixture model approach.

    The free parameters in the model are:

    * slope: slope of the fitted line.

    * intercept: intercept of the fitted line.

    * intrinsic scatter ('sigma_intrinsic'): Hogg, Bovy, Lang (2010):
      "intrinsic Gaussian variance V, orthogonal to the line."

    * outlier mean ('y_outlier'): mean of outliers.

    * outlier scatter ('sigma_outlier'): standard deviation of outliers.

    * outlier fraction ('outlier_fraction'): fraction of outliers in data.
      Hogg, Bovy, Lang (2010): "the probability that a data point is bad (or,
      more properly, the amplitude of the bad-data distribution function in the
      mixture)."


    Parameters
    ----------

    data : np.ndarray
        Should have shape (N,4) (if no covariances on errors) or (N,5) (if
        covariant errors). Should be in the order (x, y, dx, dy) or
        (x, y, dx, dy, dxy).

    priorlimits : np.ndarray
        Upper and lower values for each of the model parameters. The
        outlier fraction which has a flat prior between 0 and 1. The limits
        should be provided in the order:
        [slope, intercept, intrinsic scatter, outlier mean, outlier scatter,
        outlier fraction] (so that the array has 6x2 elements).

    nwalkers : int
        The number of emcee walkers to use in the fit.

    nsteps : int
        The number of steps each walker should take in the MCMC.

    Returns
    -------

    samples : np.array
        All the chains (flattened) minus the burn-in period.

    point_estim : np.array
        (16th, 50th, 84th) percentile for each of the 6 fitted parameters in
        the order: (slope, intercept, intrinsic scatter, outlier mean,
        outlier deviation, outlier fraction).

    """

    mv = [
        (moves.DESnookerMove(), 0.1),
        (moves.DEMove(), 0.9 * 0.9),
        (moves.DEMove(gamma0=1.0), 0.9 * 0.1),
    ]
    # mv = moves.StretchMove()

    # Unpack and check data.
    if data.shape[1] == 4:
        # No correlations on the errors
        x, y, dx, dy = data.T
        dxy = np.zeros_like(x)
    elif data.shape[1] == 5:
        # Data with dxy correlations
        x, y, dx, dy, dxy = data.T
    else:
        raise ValueError("'data' must have 4 or 5 columns, not {}. \
                Try transposing your data.".format(data.shape[1]))

    # Run the sampler hiding some annoying warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # The number of dimensions is fixed.
        ndim = 6

        # # Estimate initial values using DE algorithm.
        # def minfunc(model):
        #     return -full_posterior(model, x, y, dx, dy, dxy, priorlimits)
        # # bmin, bmax = priorlimits[0::2], priorlimits[1::2]
        # bmin, bmax = np.array(priorlimits).T
        # bounds = list(zip(*[bmin, bmax]))
        # pstart = DE(minfunc, bounds, maxiter=5000).x
        # # Sample ball around the max posterior point.
        # p0 = emcee.utils.sample_ball(
        #     pstart, 0.01 * np.ones_like(pstart), size=nwalkers)
        # # Make sure there are no negative outlier fractions
        # p0[:, -1] = np.abs(p0[:, -1])

        slope0 = np.random.uniform(*priorlimits[0], nwalkers)
        intercept0 = np.random.uniform(*priorlimits[1], nwalkers)
        i_scatter0 = np.random.uniform(*priorlimits[2], nwalkers)
        out_mean0 = np.random.uniform(*priorlimits[3], nwalkers)
        out_scatter0 = np.random.uniform(*priorlimits[4], nwalkers)
        out_fraction0 = np.random.uniform(*priorlimits[5], nwalkers)
        p0 = np.array([
            slope0, intercept0, i_scatter0, out_mean0, out_scatter0,
            out_fraction0]).T

        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, full_posterior,
            args=[x, y, dx, dy, dxy, priorlimits], moves=mv)

        for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
            if i % 100 and i < (nsteps - 1):
                continue
            update_progress.updt(nsteps, i + 1)

    # tau = sampler.get_autocorr_time(tol=0)
    # tau_mean = np.nanmean(tau, 0)

    # print("\nTaus: ", tau)
    # print("Mean Tau: {:.0f}".format(tau_mean))
    # print("Neff: {:.0f}".format(flat_samples.shape[0] / tau_mean))

    return sampler


def full_posterior(params, x, y, dx, dy, dxy, priorlimits):
    """
    The log-posterior of the data given the full mixture model of the linear
    function and the outlier distribution.

    Parameters
    ----------

    params : np.ndarray or list
        [slope,intercept, intrinsic scatter, outlier mean,
            outlier standard deviation, outlier fraction]

    Returns
    -------

    float
        The posterior of the parameters given the data.
    """
    if log_priors(params, priorlimits) == -np.inf:
        return -np.inf
    else:
        return log_priors(params, priorlimits) +\
            full_log_likelihood(params, x, y, dx, dy, dxy)


def log_priors(params, priorlimits):
    """
    Prior probabilities on the parameters, given upper and lower limits on each
    parameter. Jeffreys priors are used for the intrinsic and outlier standard
    deviations, and a prior that is flat in Arctan(slope) is used for the
    slope. For everything else, priors are uniform within the given limits.

    Parameters
    ----------

    params : np.ndarray or list
        [slope, intercept, intrinsic scatter, outlier mean, outlier standard
         deviation, outlier fraction]

    Returns
    -------

    float
        The prior density of these parameters.
    """
    m, b, sigma_intrinsic, y_outlier, sigma_outlier, outlier_fraction = params
    mlo, mhi, blo, bhi, silo, sihi, yolo, yohi, solo, sohi, oflo, ofhi =\
        [item for sublist in priorlimits for item in sublist]

    if m < mlo or m > mhi or b < blo or b > bhi or sigma_intrinsic < silo or\
        sigma_intrinsic > sihi or sigma_outlier < solo or sigma_outlier > sohi\
        or y_outlier < yolo or y_outlier > yohi or outlier_fraction < oflo or\
            outlier_fraction > ofhi:
        return -np.inf
    else:
        return -np.log(1. + m * m) - np.log(sigma_intrinsic) -\
            np.log(sigma_outlier)


def full_log_likelihood(params, x, y, dx, dy, dxy):
    """
    The log-likelihood of the data given the full mixture model of the linear
    function and the outlier distribution.

    This is basically E1. (17) in Hogg, Bovy, Lang (2010), accounting for the
    intrinsic scatter term in 'likelihood_line()'.

    Returns
    -------

    float
        The likelihood of the data given this set of model parameters.
    """
    m, b, sigma_intrinsic, y_outlier, sigma_outlier, outlier_fraction = params

    lkl_line = likelihood_line([m, b, sigma_intrinsic], x, y, dx, dy, dxy)
    out_dist = outlier_distribution(
        [y_outlier, sigma_outlier], x, y, dx, dy, dxy)

    return np.sum(np.log(
        (1. - outlier_fraction) * lkl_line + outlier_fraction * out_dist))


def likelihood_line(params, x, y, dx, dy, dxy):
    """
    Likelihood for the linear function.

    Returns
    -------

    float
        The likelihood of the data given this set of model parameters.
    """
    m, b, sigma_intrinsic = params
    theta = np.arctan(m)

    sint, cost = np.sin(theta), np.cos(theta)

    # Perpendicular distance to the line
    delta = -sint * x + cost * y - cost * b

    # Projection of covariance matrix along line
    Sigma_dd = sint**2. * dx**2. - np.sin(2. * theta) * dxy + cost**2. * dy**2.

    lkl_line = (2. * np.pi * (Sigma_dd + sigma_intrinsic**2.))**-.5 *\
        np.exp(-delta**2. / (2. * (Sigma_dd + sigma_intrinsic**2.)))

    return lkl_line


def outlier_distribution(params, x, y, dx, dy, dxy):
    """
    The likelihood for the outlier distribution, which is modeled as a uniform
    distribution in x and a Gaussian distribution in y with some mean y0 and
    standard deviation sigma0.

    Returns
    -------

    float
        The likelihood of the data given this set of model parameters.
    """
    y_outlier, sigma_outlier = params
    sigma_total2 = sigma_outlier**2. + dy**2.

    out_dist = (2. * np.pi * sigma_total2)**-0.5 *\
        np.exp(-.5 * (y - y_outlier)**2. / sigma_total2)

    return out_dist
