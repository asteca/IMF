
from modules import data_IO, mass_analysis, makePlot


def main():
    """
    TODO: Apply completeness corrections?
    """

    masses_type, Nruns, alpha_min, alpha_max, min_mag, max_mag, min_mass,\
        max_mass, binar_min, binar_max = data_IO.readINI()

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
        msk = (phot_bmsk[:, 1] >= min_mag) & (phot_bmsk[:, 1] <= max_mag)
        mass_mean_phot_msk, mass_std_phot_msk = mass_mean_bmsk[msk],\
            mass_std_bmsk[msk]

        # Apply mass range. Used by the likelihood analysis
        msk = (mass_mean_phot_msk >= min_mass) &\
            (mass_mean_phot_msk <= max_mass)
        mass_mean_mass_msk, mass_std_mass_msk = mass_mean_phot_msk[msk],\
            mass_std_phot_msk[msk]

        # Maximum likelihood analysis
        alpha_lkl, alpha_bootstrp, alpha_ranges = mass_analysis.maxLkl(
            Nruns, min_mass, max_mass, mass_mean_phot_msk, mass_mean_mass_msk,
            mass_std_mass_msk, (alpha_min, alpha_max))

        makePlot.main(
            Nruns, min_mag, max_mag, min_mass, max_mass, binar_min, binar_max,
            binar_probs, phot_bin_used, phot_bin_unused, mass_mean_phot_msk,
            mass_mean_mass_msk, alpha_lkl, alpha_bootstrp, alpha_ranges,
            alpha_min, alpha_max)


if __name__ == '__main__':
    main()
