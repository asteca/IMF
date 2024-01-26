
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy.stats import gaussian_kde
from .IMF import salpeter55


def main(
    fname, binar_cut, mag_min, mag_max,
    all_Nratios, mass_smsk, binar_probs, phot_bin_used, phot_bin_unused,
    mass_phot_msk, alpha_lkl, alpha_std, LSF_fits
):
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
    #
    plt.xlabel("Binary system probability")
    binar_min, binar_max = 0., binar_cut
    plt.axvline(x=binar_min, c='r')
    plt.axvline(x=binar_max, c='r')
    ax.axvspan(binar_min, binar_max, alpha=0.1, color='red')
    #
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
    #
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
    plt.scatter(*phot_bin_used, c='k', lw=0, alpha=.5,
                label=f"Single systems ({phot_bin_used.shape[1]})")
    plt.scatter(*phot_bin_unused, c='r', lw=0, alpha=.5,
                label=f"Binary systems ({phot_bin_unused.shape[1]})")
    ymin, ymax = ax.get_ylim()
    plt.axhline(mag_min, c='grey', ls=':')
    plt.axhline(mag_max, c='grey', ls=':')
    plt.legend(fontsize=8)
    plt.ylim(ymin, ymax)
    plt.gca().invert_yaxis()
    plt.xlabel("Color")
    plt.ylabel("Mag")

    ax = plt.subplot(gs[2:4, 0:2])
    plt.hist(mass_smsk, 20, alpha=.5, label=f"Single systems ({len(mass_smsk)})")
    plt.hist(mass_phot_msk, 20, alpha=.5,
             label=f"Single systems in mag range ({len(mass_phot_msk)})")
    plt.legend()
    plt.xlabel("Mass")

    ax = plt.subplot(gs[2:4, 2:6])
    # plt.title("Mass range: [{:.2f}, {:.2f}]".format(mass_min, mass_max))
    ax.grid(ls=':', lw=.5)

    # Plot Likelihood results
    # Best fit line. Use the (last) LSF intercept for vertical alignment
    intercept = LSF_fits[-1][4]
    x0 = np.linspace(mass_phot_msk.min(), mass_phot_msk.max(), 100)
    y_vals_log = 10**intercept * x0**(-alpha_lkl)
    txt = r"$\alpha_{{Lkl}}={:.3f}\pm{:.3f}$".format(
        alpha_lkl, alpha_std)
    plt.plot(x0, y_vals_log, c='k', lw=2, ls='--', label=txt, zorder=5)

    for bins, xx, yy, alpha, intercept, alpha_std in LSF_fits:
        # Scatter plot
        plt.scatter(
            xx, yy, s=100 - 2 * bins, marker='x', alpha=.8,
            label=r"$\alpha_{{LSF}}={:.3f}\pm{:.3f}$ (N_bins={})".format(
                alpha, alpha_std, bins))
        # Extract limits from this plot
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

    m = np.linspace(.1, 100, 1000)
    label, imf = 'Salpeter 1955', salpeter55
    x_norm, y_norm = .5*(xmin + xmax), np.median(y_vals_log)
    plt.plot(m, y_norm * imf(m) / imf(x_norm), ls='-.', c='purple', label=label)
    # plt.axvline(x=mass_min, c='grey', ls=':', lw=1.5)
    # plt.axvline(x=mass_max, c='grey', ls=':', lw=1.5)
    density = gaussian_kde(mass_phot_msk)
    xs = np.linspace(mass_phot_msk.min(), mass_phot_msk.max(), 100)
    plt.plot(xs, density(xs), c='k', label="KDE of all masses", ls=':')
    # region where the mass break in the IMF is usually placed
    ax.axvspan(.5, 1., alpha=0.1, color='grey')
    plt.xlabel(r"$m\,[M_{\odot}]$")
    plt.ylabel(r"$\xi(m) \Delta m$")
    plt.minorticks_on()
    # plt.loglog()
    plt.yscale("log")
    # ax.set_xticks([0.1, 0.5, 1, 2, 5])
    # ax.set_yticks([])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    # The x,y limits need to be better estimated
    plt.xlim(max(0.01, xmin), xmax + .1)
    plt.ylim(max(0.001, ymin), ymax + 5)
    plt.legend(fontsize=9)

    fig.tight_layout()
    plt.savefig(
        "output/{}.png".format(fname), dpi=150, bbox_inches='tight')

    print("Finished")
