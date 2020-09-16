
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import arviz as az
import pymc3 as pm
from . import IMF


def main(
    lkl_alpha, lkl_intercep, x_edges0_lsf, y_vals0_lsf, x_edges_lsf,
        y_vals_lsf, pars_lsf):
    # mcmc_mode, x_edges0, y_vals0, x_edges, y_vals, sampler, flat_samples,
    # az_data, summ,):
    """
    """
    print("Plotting..")

    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)
    #
    IMFPlot(ax, lkl_alpha, lkl_intercep)
    LSFPlot(x_edges0_lsf, y_vals0_lsf, x_edges_lsf, y_vals_lsf, pars_lsf)
    #
    ax.loglog()
    # Auto scale:
    # https://stackoverflow.com/q/21920233/1391441
    ax.set_xticks([0.5, 1, 2, 5])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

    plt.minorticks_on()
    plt.legend()
    plt.show()
    fig.tight_layout()
    plt.savefig("output/{}_IMF.png", dpi=150, bbox_inches='tight')

    # tracePlot(mcmc_mode, az_data)
    # if mcmc_mode == 'emcee':
    #     slope_name = 'alpha'
    # elif mcmc_mode == 'pymc3':
    #     slope_name = 'x'
    # convergPlots(mcmc_mode, az_data, slope_name)

    print("Finished")


def tracePlot(mcmc_mode, az_data):
    """
    """
    fig = plt.figure(figsize=(20, 24))
    az.plot_trace(az_data)
    fig.tight_layout()
    plt.savefig(
        "output/{}_trace.png".format(mcmc_mode), dpi=150, bbox_inches='tight')


def convergPlots(mcmc_mode, az_data, slope_name):
    """
    """
    fig = plt.figure(figsize=(20, 24))
    gs = gridspec.GridSpec(8, 6)

    ax = plt.subplot(gs[0:2, 0:2])
    az.plot_posterior(az_data, var_names=[slope_name], ax=ax)

    ax = plt.subplot(gs[2:4, 0:2])
    az.plot_autocorr(
        az_data, var_names=[slope_name], combined=True, max_lag=200, ax=ax)

    ax = plt.subplot(gs[2:4, 2:4])
    az.plot_ess(az_data, var_names=[slope_name], kind="evolution", ax=ax)

    fig.tight_layout()
    plt.savefig(
        "output/{}_convergence.png".format(mcmc_mode), dpi=150,
        bbox_inches='tight')


def IMFPlot(
    ax, lkl_alpha, lkl_intercep, mcmc_mode, summ, x_edges0, y_vals0, x_edges, y_vals, sampler,
        flat_samples):
    """
    """
    alpha, b = summ['median'].values
    _16p, _84p = summ['16p'].values[0], summ['84p'].values[0]

    plt.title("Mode: {}".format(mcmc_mode))

    # Sampled histograms (max 1000 points)
    Nmax = 1000
    plt.scatter(x_edges[:Nmax], y_vals[:Nmax], c='grey', alpha=.25, zorder=2)
    plt.scatter(x_edges0[:Nmax], y_vals0[:Nmax], c='r', alpha=.25)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Best fit line (medians)
    y_vals_log = 10**b * x_edges**alpha
    txt = r"$\alpha={:.3f}_{{{:.3f}}}^{{{:.3f}}}$".format(-alpha, -_84p, -_16p)
    plt.plot(x_edges, y_vals_log, c='g', label=txt, zorder=5)

    # # Plot (min, max) percentiles region
    # if flat_samples is not None:
    #     x0 = np.linspace(min(x_edges), max(x_edges), 100)
    #     inds = np.random.randint(len(flat_samples), size=100)
    #     y2 = []
    #     for ind in inds:
    #         sample = flat_samples[ind]
    #         alpha, b = sample[0], sample[1]
    #         y2.append(10**b * x0**alpha)
    #     _min, _max = np.percentile(y2, [16, 84], 0)
    #     plt.fill_between(x0, _min, _max, color="C1", zorder=0, alpha=.1)
    # else:
    #     pm.plot_posterior_predictive_glm(sampler, samples=100)

    # Maximum likelihood best fit line
    y_vals_log = 10**lkl_intercep * x_edges**(-lkl_alpha)
    txt = r"$\alpha$={:.3f}".format(lkl_alpha)
    plt.plot(x_edges, y_vals_log, c='r', label=txt, zorder=5)
    import pdb; pdb.set_trace()  # breakpoint d1959e4a //


    # x0 = np.linspace(0.01, 10, 100)
    # yIMF = IMF.chabrier2001_log(x0)
    # vert_dist = np.linspace(-1., 1., 1000)
    # v_mean = []
    # for h in vert_dist:
    #     v_mean.append(np.mean(cdist(
    #         np.array([x_edges, y_vals_log]).T, np.array([x0, h + yIMF]).T)))
    # idx = np.argmin(v_mean)
    # print(idx, vert_dist[idx])
    # plt.plot(
    #     x0, yIMF + vert_dist[idx], ls=':', c='k', label="Chabrier (2001), Log")

    # plt.plot(
    #     x0, yIMF, ls=':', c='r', label="Chabrier (2001), Log")

    # plt.axvline(x=1, ls='--', c='k')
    plt.xlabel(r"$m\,[M_{\odot}]$")
    plt.ylabel(r"$\xi(m) \Delta m$")
    plt.xlim(max(0.001, xmin), xmax + 1.)
    plt.ylim(max(0.001, ymin), ymax + 10)


def LSFPlot(xx0, yy0, xx, yy, pars_lsf):
    """
    """
    alpha, b, alpha_std = pars_lsf
    plt.plot(
        xx, 10**(b) * xx**(alpha),
        c='#7BC5D4', label=r"$\alpha_{{LSF}}={:.3f}\pm{:.3f}$".format(
            -alpha, alpha_std))
    plt.scatter(xx, yy, c='cyan', s=25)  # label=r'$M_{i}^{obs}$')
    if len(xx0) > 0:
        plt.scatter(xx0, yy0, c='orange')

    # plt.axvline(x=1, ls='--', c='k')
    # plt.xlabel(r"$m\,[M_{\odot}]$")
    # plt.ylabel(r"$\xi(m) \Delta m$")
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.minorticks_on()
    # plt.legend()
