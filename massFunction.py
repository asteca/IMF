
import warnings
import numpy as np
from numpy.random import MT19937, RandomState, SeedSequence
from modules import data_IO, mass_analysis, makePlot


def main():
    """
    TODO: Apply completeness corrections?
    """
    seed = np.random.randint(100000000)
    # print("Random seed: {}".format(seed))
    RandomState(MT19937(SeedSequence(seed)))

    mag_min, mag_max, binar_cut, make_plot_flag = data_IO.readINI()

    inputfiles = data_IO.readFiles()

    for file in inputfiles:
        fname = str(file).split('/')[-1].split('.')[0]

        mass, mass_std, binar_probs, phot = data_IO.dataRead(file)

        # Mask photometry and masses given the binary probability cut
        smsk = binar_probs <= binar_cut
        mass_smsk, mass_std_smsk, phot_smsk = mass[smsk], \
            mass_std[smsk], phot[smsk]

        # Apply magnitude range.
        # mass_std_phot_msk not used because uncertainties in masses are not used
        msk = (phot_smsk[:, 1] >= mag_min) & (phot_smsk[:, 1] <= mag_max)
        mass_phot_msk, mass_std_phot_msk = mass_smsk[msk], \
            mass_std_smsk[msk]
        if msk.sum() < 10:
            warnings.warn("Less than 10 stars identified as single systems. "
                          + "Can not process", UserWarning)
            continue

        # # Apply mass range. Used by the likelihood analysis
        # msk = (mass_phot_msk >= mass_min) & (mass_phot_msk <= mass_max)
        # mass_mass_msk, mass_std_mass_msk = mass_phot_msk[msk],\
        #     mass_std_phot_msk[msk]
        # if msk.sum() < 10:
        #     warnings.warn("Less than 10 stars identified as single systems. "
        #                   + "Can not process", UserWarning)
        #     continue

        # Maximum likelihood analysis
        alpha_lkl, alpha_std = mass_analysis.maxLkl(mass_phot_msk)
        print("Cluster: {} ; alpha = {:.3f}+/-{:.3f} ({})".format(
            fname, alpha_lkl, alpha_std, msk.sum()))

        # LSF histogram fits
        LSF_fits = []
        for bins in (5, 10, 20):
            LSF_fits.append(mass_analysis.binnedIMF(mass_phot_msk, bins))

        if make_plot_flag:
            make_plot(
                mag_min, mag_max,
                binar_cut, fname, mass, binar_probs, phot, smsk, mass_smsk, phot_smsk,
                mass_phot_msk, alpha_lkl, alpha_std, LSF_fits)


def make_plot(
    mag_min, mag_max, binar_cut,
    fname, mass, binar_probs, phot, smsk, mass_smsk, phot_smsk, mass_phot_msk,
    alpha_lkl, alpha_std, LSF_fits
):
    """ """
    phot_bin_used, phot_bin_unused = phot_smsk.T, phot[~smsk].T

    all_Nratios = mass_analysis.singleBinarRatio(
        binar_cut, mass, binar_probs, phot)

    makePlot.main(
        fname, binar_cut, mag_min, mag_max, all_Nratios, mass_smsk,
        binar_probs, phot_bin_used, phot_bin_unused, mass_phot_msk,
        alpha_lkl, alpha_std, LSF_fits)


if __name__ == '__main__':
    main()
