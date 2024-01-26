
import numpy as np
from astropy.stats import histogram as ashist


def singleBinarRatio(binar_cut, mass_mean, binar_probs, phot):
    """
    """
    def probsRange(xdata):
        """
        """
        # Define 10 bins in x data
        x_range = np.linspace(xdata.min(), xdata.max(), 10)
        Nratios = []
        min_m = x_range[0]
        for mag_m in x_range[1:]:
            msk = (xdata >= min_m) & (xdata <= mag_m)
            min_m = mag_m
            binar_probs_msk = binar_probs[msk]
            msk_b = binar_probs_msk <= binar_cut
            Nsingle, Nbinar = binar_probs_msk[msk_b], binar_probs_msk[~msk_b]
            if Nsingle.size + Nbinar.size > 0:
                Nratios.append(Nbinar.size / (Nsingle.size + Nbinar.size))
            else:
                Nratios.append(0)

        return x_range, np.array(Nratios)

    all_Nratios = []
    for xdata in (phot.T[1], mass_mean):
        all_Nratios.append(probsRange(xdata))

    return all_Nratios


def maxLkl(mass, alpha_min=-1, alpha_max=5):
    """
    Method defined in Khalaj & Baumgardt (2013):
    https://academic.oup.com/mnras/article/434/4/3236/960889
    and used in Sheikhi et al. (2016):
    https://academic.oup.com/mnras/article/457/1/1028/989829
    """

    def minfunc(alpha, x, xmin, xminmax, N):
        y = abs(
            alpha - (
                1 + N / (
                    np.sum(np.log(x / xmin))
                    - N * (np.log(xminmax) / (1 - xminmax**(alpha - 1)))
                )
            )
        )
        idx = np.argmin(y)
        return alpha[idx]

    # Slope on the original sample
    N = mass.size
    xmin, xmax = mass.min(), mass.max()
    xminmax = xmax / xmin
    if alpha_min <= 1:
        # Protect from -1 divergence
        alpha_r1 = list(np.linspace(alpha_min, .95, 2500))
        alpha_r2 = list(np.linspace(1.05, alpha_max, 2500))
        alpha_vals = np.array(alpha_r1 + alpha_r2)
    else:
        alpha_vals = np.linspace(alpha_min, alpha_max, 5000)
    alpha_lkl = minfunc(alpha_vals, mass, xmin, xminmax, N)

    N = mass.size
    alpha_std = (1/np.sqrt(N)) * (
        (alpha_lkl-1)**(-2) - np.log(xminmax)**2*(
            (xminmax**(alpha_lkl-1))/(1-xminmax**(alpha_lkl-1))**2))**(-.5)

    return alpha_lkl, alpha_std


def binnedIMF(mass, bins):
    """
    """
    # Histogram of all the masses
    yy, xx = ashist(mass, bins=bins, density=True)

    # Remove possible nans
    yy[np.isnan(yy)] = 0.
    # Align x values
    xx = .5 * (xx[1:] + xx[:-1])
    # Remove empty bins. Not sure if this has any impact
    msk = yy != 0
    yy, xx = yy[msk], xx[msk]

    # Perform LSF fit in logarithmic values
    alpha, intercept, alpha_std = logLSF(xx, yy)

    return bins, xx, yy, -alpha, intercept, alpha_std


def logLSF(xx, yy):
    """
    """
    # Least squares fit
    y_log = np.log10(yy)
    # Mask -inf values that can appear when a bin contains 0 elements.
    msk = y_log == -np.inf
    x_log = np.log10(xx[~msk])
    y_log = y_log[~msk]
    pars, cov_matrix = np.polyfit(x_log, y_log, 1, cov=True)
    alpha, intercept = pars
    alpha_std = np.sqrt(cov_matrix[0][0])

    return alpha, intercept, alpha_std
