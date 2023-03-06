
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from astropy.stats import histogram as ashist
import seaborn as sn
from .IMF import salpeter55, millerscalo79, kroupa01, chabrier03individual,\
    chabrier03system


def main(
    fname, binar_cut, mag_min, mag_max, mass_min, mass_max, Nruns,
    all_Nratios, binar_probs, phot_bin_used, phot_bin_unused,
    mass_mean_phot_msk, mass_mean_mass_msk, alpha_lkl, alpha_bootstrp,
        alpha_ranges, alpha_min, alpha_max, sampled_IMFs):
    """
    """
    # print("Plotting..")
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(8, 8)

    ax = plt.subplot(gs[0:2, 0:2])
    ax.minorticks_on()
    plt.title("P_cut = {:.2f}".format(binar_cut))
    hist, edges = np.histogram(binar_probs, bins=20)
    # Histogram normalized to max value equal to 1
    plt.bar(
        (edges[1:] + edges[:-1]) * .5, hist / hist.max(),
        width=(edges[1] - edges[0]))

    plt.xlabel("Binary system probability")
    binar_min, binar_max = 0., binar_cut
    plt.axvline(x=binar_min, c='r')
    plt.axvline(x=binar_max, c='r')
    ax.axvspan(binar_min, binar_max, alpha=0.1, color='red')

    xx, yy = np.linspace(.01, .99, 100), []
    Nt = len(binar_probs)
    for bl in xx:
        yy.append((binar_probs >= bl).sum() / Nt)
    idx = np.argmin(abs(xx - binar_cut))
    plt.plot(xx, yy, c='k', label="Binary fraction ({:.2f})".format(yy[idx]),
             zorder=6)
    plt.legend(fontsize=8)
    # ax.set_yticks([])
    plt.xlim(0., 1.)

    ax = plt.subplot(gs[0:1, 2:4])
    x_range, Nratios = all_Nratios[0]
    width = .8 * (x_range[1] - x_range[0])
    # Center bins
    x_range = .5 * (x_range[1:] + x_range[:-1])
    single_probs = 1. - Nratios
    plt.bar(x_range, single_probs, label="Prob single", width=width)
    plt.bar(x_range, Nratios, bottom=single_probs, label="Prob binar",
            width=width)
    plt.legend(fontsize=8)
    plt.xlabel("Magnitude")

    ax = plt.subplot(gs[1:2, 2:4])
    x_range, Nratios = all_Nratios[1]
    width = .8 * (x_range[1] - x_range[0])
    # Center bins
    x_range = .5 * (x_range[1:] + x_range[:-1])
    single_probs = 1. - Nratios
    plt.bar(x_range, single_probs, width=width)
    plt.bar(x_range, Nratios, bottom=single_probs, width=width)
    plt.xlabel("Mass")
    plt.gca().invert_xaxis()

    ax = plt.subplot(gs[0:2, 4:6])
    ax.minorticks_on()
    plt.scatter(*phot_bin_unused, c='r', lw=0, alpha=.5)
    plt.scatter(*phot_bin_used, c='k', lw=0, alpha=.5,
                label="N={}".format(phot_bin_used.shape[1]))
    ymin, ymax = ax.get_ylim()
    plt.axhline(mag_min, c='grey', ls=':')
    plt.axhline(mag_max, c='grey', ls=':')
    plt.legend(fontsize=8)
    plt.ylim(ymin, ymax)
    plt.gca().invert_yaxis()
    plt.xlabel("Color")
    plt.ylabel("Mag")

    alpha_16p, alpha_50p, alpha_84p = np.percentile(
        alpha_bootstrp, (16, 50, 84))
    alpha_mean, alpha_std = alpha_bootstrp.mean(), alpha_bootstrp.std()

    ax = plt.subplot(gs[2:4, 0:2])
    ax.minorticks_on()
    plt.title(r"Bootstrap distribution (N={})".format(Nruns))
    sn.kdeplot(alpha_bootstrp)
    # plt.hist(alpha_bootstrp, 20, color='grey', density=True)
    plt.axvline(alpha_16p, c='orange', ls=':', label="16p")
    plt.axvline(alpha_50p, c='red', ls=':', label="median")
    plt.axvline(alpha_84p, c='orange', ls=':', label="84p")
    plt.axvline(alpha_mean, c='green', ls='--', label="mean")
    boot_bias = alpha_mean - alpha_lkl
    txt = "Original sample\n(bias={:.3f})".format(boot_bias)
    plt.axvline(alpha_lkl, c='k', ls='-', label=txt)
    plt.legend(fontsize=8)
    plt.xlabel(r"$\alpha$")
    plt.xlim(alpha_mean - 3 * alpha_std, alpha_mean + 3 * alpha_std)

    ax = plt.subplot(gs[2:4, 2:6])
    plt.title("Mass range: [{:.2f}, {:.2f}]".format(mass_min, mass_max))
    ax.grid(ls=':', lw=.5)

    # LSF histogram fits
    for bins in (5, 10, 25):
        xmin, xmax, ymin, ymax, intercept = binnedIMF(
            ax, mass_mean_mass_msk, bins)

    # Bias corrected value
    # On bootstrap bias:  https://stats.stackexchange.com/a/488223/10416
    alpha_lkl -= boot_bias

    # Plot Likelihood results
    # Best fit line. Use the (last) LSF intercept for vertical alignment
    x0 = np.linspace(mass_mean_mass_msk.min(), mass_mean_mass_msk.max(), 100)
    y_vals_log = 10**intercept * x0**(-alpha_lkl)
    txt = r"$\alpha_{{Lkl}}={:.3f}\pm{:.3f}$".format(
        alpha_lkl, alpha_std)
    plt.plot(x0, y_vals_log, c='k', lw=2, ls='--', label=txt, zorder=5)

    m = np.linspace(.1, 100, 1000)
    for label, imf in zip(
            ['Salpeter (55)', 'Miller-Scalo (79)', 'Kroupa (01)',
             'Chabrier (03, indiv)', 'Chabrier (03, system)'],
            [salpeter55, millerscalo79, kroupa01, chabrier03individual,
             chabrier03system]):
        plt.plot(m, imf(m) / imf(1), ls='--', label=label)

    plt.axvline(x=mass_min, c='grey', ls=':', lw=1.5)
    plt.axvline(x=mass_max, c='grey', ls=':', lw=1.5)
    sn.kdeplot(mass_mean_phot_msk, label="KDE of all masses", ls=':')
    # region where the mass break in the IMF is usually placed
    ax.axvspan(.5, 1., alpha=0.1, color='grey')
    plt.xlabel(r"$m\,[M_{\odot}]$")
    plt.ylabel(r"$\xi(m) \Delta m$")
    plt.minorticks_on()
    plt.loglog()
    ax.set_xticks([0.1, 0.5, 1, 2, 5])
    # ax.set_yticks([])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    # The x,y limits need to be better estimated
    plt.xlim(max(0.01, xmin), xmax + .5)
    plt.ylim(max(0.001, ymin), ymax + 10)
    plt.legend(fontsize=9)

    ax = plt.subplot(gs[4:6, 0:2])
    ax.minorticks_on()
    ax.tick_params(labelbottom=False)
    plt.title("Slope values for several mass ranges")
    ax.grid(ls=':', lw=.5)
    for i, (k, v) in enumerate(alpha_ranges.items()):
        v -= boot_bias
        if i == 0:
            plt.scatter(i, v, s=50, c='k', marker='o', label="{}".format(k))
        else:
            plt.scatter(i, v, marker='x', label="{}".format(k))
    plt.legend(ncol=3, fontsize=8, framealpha=.3)
    plt.axhline(alpha_min, c='r', ls='--')
    plt.axhline(alpha_max, c='r', ls='--')

    ax = plt.subplot(gs[4:6, 2:6])
    plt.title(r"Slope for sampled IMFs, $M_T=${:.0f}".format(
        mass_mean_mass_msk.sum()))
    ax.grid(ls=':', lw=.5)
    for k, (mass_IMF, Lkl_IMF, bootstrp_IMF) in sampled_IMFs.items():
        # LSF histogram fits
        for bins in (25,):
            xmin, xmax, ymin, ymax, intercept = binnedIMF(
                ax, mass_IMF, bins)
        # Plot Likelihood results
        x0 = np.linspace(mass_IMF.min(), mass_IMF.max(), 100)
        y_vals_log = 10**intercept * x0**(-Lkl_IMF)
        txt = r"$\alpha_{{Lkl}}={:.3f}\pm{:.3f}$, {}".format(
            Lkl_IMF, bootstrp_IMF.std(), k)
        plt.plot(x0, y_vals_log, label=txt, zorder=5)

    ax.axvspan(.5, 1., alpha=0.1, color='grey')
    plt.xlabel(r"$m\,[M_{\odot}]$")
    plt.ylabel(r"$\xi(m) \Delta m$")
    plt.minorticks_on()
    plt.loglog()
    ax.set_xticks([0.1, 0.5, 1, 2, 5])
    # ax.set_yticks([])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    plt.xlim(max(0.01, xmin), xmax + .5)
    plt.ylim(max(0.001, ymin), ymax + 10)
    plt.legend(ncol=2, fontsize=9)

    fig.tight_layout()
    plt.savefig(
        "output/{}.png".format(fname), dpi=150, bbox_inches='tight')

    print("Finished")


def binnedIMF(ax, mass_mean_mass_msk, bins):
    """
    """
    # Histogram of all the masses
    yy, xx = ashist(mass_mean_mass_msk, bins=bins, density=True)
    # Remove possible nans
    yy[np.isnan(yy)] = 0.
    # Align x values
    xx = .5 * (xx[1:] + xx[:-1])
    # Remove empty bins. Not sure if this has any impact
    msk = yy != 0
    yy, xx = yy[msk], xx[msk]

    # Perform LSF fit in logarithmic values
    alpha, intercept, alpha_std = logLSF(xx, yy)

    # Scatter plot
    plt.scatter(
        xx, yy, s=100 - 2 * bins, marker='x', alpha=.8,
        label=r"$\alpha_{{LSF}}={:.3f}\pm{:.3f}$ (b={})".format(
            -alpha, alpha_std, bins))
    # Extract limits from this plot
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # # Line fit
    # y_vals_log = 10**intercept * xx**alpha
    # plt.plot(xx, y_vals_log, ls=":", lw=1.5)

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
