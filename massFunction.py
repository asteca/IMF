
import numpy as np
from numpy.random import MT19937, RandomState, SeedSequence
from modules import data_IO, mass_analysis, makePlot


def main():
    """
    TODO: Apply completeness corrections?
    """
    seed = np.random.randint(100000000)
    print("Random seed: {}".format(seed))
    RandomState(MT19937(SeedSequence(seed)))

    masses_type, Nruns, alpha_min, alpha_max, mag_min, mag_max, mass_min,\
        mass_max, binar_min, binar_max = data_IO.readINI()

    inputfiles = data_IO.readFiles()

    for file in inputfiles:
        fname = str(file).split('/')[1].split('.')[0]
        print("Processing: {}".format(fname))
        print("Reading data from file")
        mass_mean, mass_std, binar_probs, phot = data_IO.dataRead(
            masses_type, file)

        # Mask photometry and masses given the binary fraction probabilities
        # range
        bmsk = (binar_probs >= binar_min) & (binar_probs <= binar_max)
        mass_mean_bmsk, mass_std_bmsk, phot_bmsk = mass_mean[bmsk],\
            mass_std[bmsk], phot[bmsk]
        # For plotting
        phot_bin_used, phot_bin_unused = phot_bmsk.T, phot[~bmsk].T

        # Apply magnitude range, This keeps the full mass range
        msk = (phot_bmsk[:, 1] >= mag_min) & (phot_bmsk[:, 1] <= mag_max)
        mass_mean_phot_msk, mass_std_phot_msk = mass_mean_bmsk[msk],\
            mass_std_bmsk[msk]

        # Apply mass range. Used by the likelihood analysis
        msk = (mass_mean_phot_msk >= mass_min) &\
            (mass_mean_phot_msk <= mass_max)
        mass_mean_mass_msk, mass_std_mass_msk = mass_mean_phot_msk[msk],\
            mass_std_phot_msk[msk]

        # Maximum likelihood analysis
        alpha_lkl, alpha_bootstrp, alpha_ranges = mass_analysis.maxLkl(
            Nruns, mass_min, mass_max, mass_mean_phot_msk, mass_mean_mass_msk,
            mass_std_mass_msk, (alpha_min, alpha_max))

        makePlot.main(
            Nruns, mag_min, mag_max, mass_min, mass_max, binar_min, binar_max,
            binar_probs, phot_bin_used, phot_bin_unused, mass_mean_phot_msk,
            mass_mean_mass_msk, alpha_lkl, alpha_bootstrp, alpha_ranges,
            alpha_min, alpha_max)


if __name__ == '__main__':
    main()
