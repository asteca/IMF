
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from astropy.stats import histogram as ashist
import seaborn as sn
# from . import IMF


def main(
    Nruns, min_mag, max_mag, min_mass, max_mass, binar_min, binar_max,
    binar_probs, phot_bin_used, phot_bin_unused, mass_mean_phot_msk,
    mass_mean_mass_msk, alpha_lkl, alpha_bootstrp, alpha_ranges, alpha_min,
        alpha_max):
    """
    """
    print("Plotting..")
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(8, 8)

    ax = plt.subplot(gs[0:2, 0:2])
    ax.minorticks_on()
    plt.title("P_min, P_max = [{:.2f}, {:.2f}]".format(binar_min, binar_max))
    _, edges, _ = plt.hist(binar_probs, 20)
    plt.xlabel("Binary system probability")
    plt.axvline(x=binar_min, c='r')
    plt.axvline(x=binar_max, c='r')
    ax.axvspan(binar_min, binar_max, alpha=0.2, color='red')
    ax.set_yticks([])
    plt.xlim(0., 1.)

    ax = plt.subplot(gs[0:2, 2:4])
    ax.minorticks_on()
    plt.scatter(*phot_bin_unused, s=5, c='k', lw=0, alpha=.5)
    plt.scatter(*phot_bin_used, s=10, c='r', lw=0, alpha=.5,
                label="N={}".format(phot_bin_used.shape[1]))
    ymin, ymax = ax.get_ylim()
    plt.axhline(min_mag, c='grey', ls=':')
    plt.axhline(max_mag, c='grey', ls=':')
    plt.legend(fontsize=8)
    plt.ylim(ymin, ymax)
    plt.gca().invert_yaxis()

    alpha_16p, alpha_50p, alpha_84p = np.percentile(
        alpha_bootstrp, (16, 50, 84))
    alpha_mean, alpha_std = alpha_bootstrp.mean(), alpha_bootstrp.std()

    ax = plt.subplot(gs[0:2, 4:6])
    ax.minorticks_on()
    plt.title(r"Bootstrap distribution (N={})".format(Nruns))
    sn.kdeplot(alpha_bootstrp)
    # plt.hist(alpha_bootstrp, 20, color='grey', density=True)
    plt.axvline(alpha_16p, c='orange', ls=':', label="16p")
    plt.axvline(alpha_50p, c='red', ls=':', label="median")
    plt.axvline(alpha_84p, c='orange', ls=':', label="84p")
    plt.axvline(alpha_mean, c='green', ls='--', label="mean")
    txt = "Original sample\n(bias={:.3f})".format(alpha_mean - alpha_lkl)
    plt.axvline(alpha_lkl, c='k', ls='-', label=txt)
    plt.legend(fontsize=8)
    plt.xlabel(r"$\alpha$")

    ax = plt.subplot(gs[2:4, 0:4])
    plt.title("Mass range: [{:.2f}, {:.2f}]".format(min_mass, max_mass))
    ax.grid(ls=':', lw=.5)
    for bins in (5, 10, 25, 50):
        xmin, xmax, ymin, ymax, intercept = binnedIMF(
            ax, mass_mean_mass_msk, bins)

    # Plot Likelihood results
    # Best fit line. Use the (last) LSF intercept for vertical alignment
    x0 = np.linspace(mass_mean_mass_msk.min(), mass_mean_mass_msk.max(), 100)
    y_vals_log = 10**intercept * x0**(-alpha_lkl)
    txt = r"$\alpha_{{Lkl}}={:.3f}\pm{:.3f}$".format(
        alpha_lkl, alpha_std)
    plt.plot(x0, y_vals_log, c='k', lw=2, ls='--', label=txt, zorder=5)

    plt.axvline(x=min_mass, c='grey', ls=':', lw=1.5)
    plt.axvline(x=max_mass, c='grey', ls=':', lw=1.5)
    sn.kdeplot(mass_mean_phot_msk, label="KDE of all masses")
    # region where the mass break in the IMF is usually placed
    ax.axvspan(.5, 1., alpha=0.1, color='grey')
    plt.xlabel(r"$m\,[M_{\odot}]$")
    plt.ylabel(r"$\xi(m) \Delta m$")
    plt.xlim(max(0.001, xmin), xmax + 1.)
    plt.ylim(max(0.001, ymin), ymax + 10)
    plt.minorticks_on()
    plt.loglog()
    ax.set_xticks([0.5, 1, 2, 5])
    ax.set_yticks([])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    plt.legend()

    ax = plt.subplot(gs[2:4, 4:6])
    ax.minorticks_on()
    ax.tick_params(labelbottom=False)
    plt.title("Slope values for several mass ranges")
    ax.grid(ls=':', lw=.5)
    for i, (k, v) in enumerate(alpha_ranges.items()):
        if i == 0:
            plt.scatter(i, v, s=50, c='k', marker='o', label="{}".format(k))
        else:
            plt.scatter(i, v, marker='x', label="{}".format(k))
    plt.legend(ncol=3, fontsize=8, framealpha=.3)
    plt.axhline(alpha_min, c='r', ls='--')
    plt.axhline(alpha_max, c='r', ls='--')

    fig.tight_layout()
    plt.savefig("output/IMF.png", dpi=150, bbox_inches='tight')

    print("Finished")


def binnedIMF(ax, mass_mean_mass_msk, bins):
    """
    """
    # Histogram of all the masses
    yy, xx = ashist(mass_mean_mass_msk, bins=bins, density=True)
    # Align x values
    xx = .5 * (xx[1:] + xx[:-1])
    # Remove empty bins. Not sure if this has any impact
    msk = yy != 0
    yy, xx = yy[msk], xx[msk]

    # Perform LSF fit in logarithmic values
    alpha, intercept, alpha_std = logLSF(xx, yy)

    # Scatter plot
    plt.scatter(
        xx, yy, s=70 - bins, marker='x', alpha=.8,
        label=r"$\alpha_{{LSF}}={:.3f}\pm{:.3f}$ (b={})".format(
            -alpha, alpha_std, bins))
    # Extract limits from this plot
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Line fit
    y_vals_log = 10**intercept * xx**alpha
    plt.plot(xx, y_vals_log, ls=":")

    return xmin, xmax, ymin, ymax, intercept


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
