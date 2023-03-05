
import numpy as np


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


def maxLkl(mass, alpha_bounds, bootsrp_args, mass_full_range):
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
    alpha_vals = np.linspace(alpha_bounds[0], alpha_bounds[1], 5000)
    alpha_lkl = minfunc(alpha_vals, mass, xmin, xminmax, N)

    Nruns, mass_min, mass_max, mass_std = bootsrp_args

    # if mass_std is not None:
    #     print("Bootstraping slope's distribution")
    # Estimate bootstrap for alpha using DE algorithm.
    alpha_lst = []
    for _ in range(Nruns):
        # Re-sample mass values
        if mass_std is not None:
            mass_sample = np.random.normal(mass, mass_std)
            # mass_sample = np.random.choice(mass, mass.size)
        else:
            mass_sample = np.random.choice(mass, mass.size)

        # Apply mass range
        msk = (mass_sample >= mass_min) & (mass_sample <= mass_max)
        mass_sample = mass_sample[msk]

        # Masses can not be smaller than this value
        mass_sample = np.clip(mass_sample, a_min=0.08, a_max=None)
        N, xmin = mass_sample.size, mass_sample.min()
        xminmax = mass_sample.max() / xmin
        alpha_lst.append(minfunc(
            alpha_vals, mass_sample, xmin, xminmax, N))

    # intercept = (1 - alpha)/(mass_max**(1 - alpha)-mass_min**(1 - alpha))
    alpha_bootstrp = np.array(alpha_lst)

    if mass_std is None:
        return alpha_lkl, alpha_bootstrp

    # print("Creating dictionary of slopes for different mass ranges")
    # Create dictionary of slopes. Store the Likelihood alpha obtained
    # with the filtered mass range
    alpha_ranges = {"[{:.2f}, {:.2f}]".format(
        mass.min(), mass.max()): alpha_lkl}
    # Estimate slopes using the full magnitude range
    mr = np.linspace(mass_full_range.min(), mass_full_range.max(), 5)
    mr1, mr2, mr3 = (mr[0], mr[1]), (mr[0], mr[2]), (mr[0], mr[3])
    mr4, mr5, mr6 = (mr[1], mr[2]), (mr[1], mr[3]), (mr[1], mr[4])
    mr7, mr8, mr9 = (mr[2], mr[3]), (mr[2], mr[4]), (mr[3], mr[4])
    mr10 = (mr[0], mr[4])
    for mrx in (mr1, mr2, mr3, mr4, mr5, mr6, mr7, mr8, mr9, mr10):
        mi, mh = mrx
        msk = (mass_full_range >= mi) & (mass_full_range < mh)

        mass_msk = mass_full_range[msk]
        N, xmin = mass_msk.size, mass_msk.min()
        xminmax = mass_msk.max() / xmin
        if abs(xminmax - 1.) < 0.001:
            continue

        alpha_ranges["[{:.2f}, {:.2f})".format(mi, mh)] = minfunc(
            alpha_vals, mass_msk, xmin, xminmax, N)

    return alpha_lkl, alpha_bootstrp, alpha_ranges
