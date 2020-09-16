
import time as t
from modules import data_IO, mass_analysis, arvizData, makePlot


def main():
    """
    TODO in place for #96
    1. use the LF maximum magnitude cut, otherwise the incompleteness of low
    mass stars impacts on the slope estimation.
    2. perhaps fit two distinct lines: m<1 and m>1? This is where most IMFs
    seem to make a break in the slope.
    4. Apply completeness corrections?
    """

    Nruns, min_mass, max_mass, mcmc_mode, burn_frac, nsteps, nwalkers =\
        data_IO.readINI()

    inputfiles = data_IO.readFiles()

    for file in inputfiles:
        print(file)
        mass_data = data_IO.dataRead(file)
        mass_mean, mass_std = mass_data['Mass_mu'], mass_data['Mass_std']

        lkl_alpha, lkl_intercep = mass_analysis.maxLkl(
            mass_mean, mass_std, Nruns)

        # # Generate Nruns histograms by randomly re-sampling the masses
        # x_edges0, y_vals0, x_edges, y_vals = mass_analysis.massResample(
        #     Nruns, min_mass, max_mass, mass_mean, mass_std)

        # st = t.time()
        # if mcmc_mode == 'pymc3':
        #     print("pyMC3 analysis")
        #     trace = mass_analysis.pyMC3Run(x_edges, y_vals, nsteps)
        # elif mcmc_mode == 'emcee':
        #     print("emcee analysis")
        #     trace = mass_analysis.emceeRun(x_edges, y_vals, nsteps, nwalkers)
        # print("Time elapsed: {:.0f} s".format(t.time() - st))

        # # Prepare MCMC data in Arviz format
        # az_data, summ, flat_samples = arvizData.proc(
        #     mcmc_mode, burn_frac, nsteps, trace)

        x_edges0_lsf, y_vals0_lsf, x_edges_lsf, y_vals_lsf, pars_lsf =\
            mass_analysis.LSFRun(min_mass, max_mass, mass_mean)

        makePlot.main(
            lkl_alpha, lkl_intercep,
            # mcmc_mode, x_edges0, y_vals0, x_edges, y_vals, trace,
            # flat_samples, az_data, summ,
            x_edges0_lsf, y_vals0_lsf, x_edges_lsf, y_vals_lsf, pars_lsf)


if __name__ == '__main__':
    main()
