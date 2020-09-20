
import numpy as np
# from scipy.optimize import differential_evolution as DE


def maxLkl(
    Nruns, min_mass, max_mass, full_mr_mean, mass_mean, mass_std, alpha_bounds,
        maxiter=10000):
    """
    Method defined in Khalaj & Baumgardt (2013):
    https://academic.oup.com/mnras/article/434/4/3236/960889
    and used in Sheikhi et al. (2016):
    https://academic.oup.com/mnras/article/457/1/1028/989829
    """
    # min_mass, max_mass = .01, 100

    def minfunc(alpha, x, xmin, xminmax, N):
        y = abs(alpha - (1 + N / (
                np.sum(np.log(x / xmin)) -
                N * (np.log(xminmax) / (1 - xminmax**(alpha - 1))))))
        idx = np.argmin(y)
        return alpha[idx]

    # Range for alpha
    # bounds = [alpha_bounds]
    N = mass_mean.size
    xmin = mass_mean.min()
    xminmax = mass_mean.max() / xmin
    alpha_vals = np.linspace(alpha_bounds[0], alpha_bounds[1], 1000)
    # y = abs(alpha_vals - (1 + N / (
    #         np.sum(np.log(mass_mean / xmin)) -
    #         N * (np.log(xminmax) / (1 - xminmax**(alpha_vals - 1))))))

    # Values obtained using the mean masses
    # args = (mass_mean, mass_mean.min(),
    #         mass_mean.max() / mass_mean.min(), mass_mean.size)
    # alpha_lkl = DE(minfunc, bounds, maxiter=maxiter, args=args).x[0]
    alpha_lkl = minfunc(alpha_vals, mass_mean, xmin, xminmax, N)

    print("Bootstraping slope's distribution")
    # Estimate bootstrap for alpha using DE algorithm.
    alpha_lst = []
    for _ in range(Nruns):
        # Re-sample mass values
        mass_sample = np.random.normal(mass_mean, mass_std)

        # Apply mass range
        msk = (mass_sample >= min_mass) & (mass_sample <= max_mass)
        mass_sample = mass_sample[msk]

        # Masses can not be smaller than this value
        mass_sample = np.clip(mass_sample, a_min=0.01, a_max=None)
        # args = (mass_sample, mass_sample.min(),
        #         mass_sample.max() / mass_sample.min(), mass_sample.size)
        # result = DE(minfunc, bounds, maxiter=maxiter, args=args)
        # alpha_lst.append(result.x[0])

        N = mass_sample.size
        xmin = mass_sample.min()
        xminmax = mass_sample.max() / xmin
        alpha_lst.append(minfunc(alpha_vals, mass_sample, xmin, xminmax, N))

    # intercept = (1 - alpha) / (max_mass**(1 - alpha) - min_mass**(1 - alpha))
    alpha_bootstrp = np.array(alpha_lst)

    print("Creating dictionary of slopes for different mass ranges")
    # Create dictionary of slopes. Store the Likelihood alpha obtained
    # with the filtered mass range
    alpha_ranges = {"[{:.2f}, {:.2f}]".format(
        mass_mean.min(), mass_mean.max()): alpha_lkl}
    # Estimate slopes using the full magnitude range
    mr = np.linspace(full_mr_mean.min(), full_mr_mean.max(), 5)
    mr1, mr2, mr3 = (mr[0], mr[1]), (mr[0], mr[2]), (mr[0], mr[3])
    mr4, mr5, mr6 = (mr[1], mr[2]), (mr[1], mr[3]), (mr[1], mr[4])
    mr7, mr8, mr9 = (mr[2], mr[3]), (mr[2], mr[4]), (mr[3], mr[4])
    mr10 = (mr[0], mr[4])
    for mrx in (mr1, mr2, mr3, mr4, mr5, mr6, mr7, mr8, mr9, mr10):
        mi, mh = mrx
        msk = (full_mr_mean >= mi) & (full_mr_mean < mh)

        mass_msk = full_mr_mean[msk]
        # args = (mass_msk, mass_msk.min(),
        #         mass_msk.max() / mass_msk.min(), mass_msk.size)
        # alpha_ranges["[{:.2f}, {:.2f})".format(mi, mh)] = DE(
        #     minfunc, bounds, maxiter=maxiter, args=args).x[0]

        N = mass_msk.size
        xmin = mass_msk.min()
        xminmax = mass_msk.max() / xmin
        alpha_ranges["[{:.2f}, {:.2f})".format(mi, mh)] = minfunc(
            alpha_vals, mass_msk, xmin, xminmax, N)

    return alpha_lkl, alpha_bootstrp, alpha_ranges
