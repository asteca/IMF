
import pickle
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d


"""
Taken from Wikipedia's IMF page

The package https://github.com/keflavich/imf has some more (I think,
24-09-2019).
"""


def salpeter55(m):
    alpha = 2.35
    return m**-alpha


def millerscalo79(m):
    return np.where(m > 1, salpeter55(m), salpeter55(1))


def chabrier03individual(m):
    k = 0.158 * np.exp(-(-np.log10(0.08))**2 / (2 * 0.69**2))
    return np.where(m <= 1, 0.158 * (1. / m) * np.exp(-(
        np.log10(m) - np.log10(0.08))**2 / (2 * 0.69**2)), k * m**-2.3)


def chabrier03system(m):
    k = 0.086 * np.exp(-(-np.log10(0.22))**2 / (2 * 0.57**2))
    return np.where(m <= 1, 0.086 * (1. / m) * np.exp(
        -(np.log10(m) - np.log10(0.22))**2 / (2 * 0.57**2)), k * m**-2.3)


def kroupa01(m):
    """
    Kroupa (2001), 'On the variation of the initial mass function', Eq (2)
    """
    return np.where(
        m < 0.08, m**-0.3,
        np.where(m < 0.5, 0.08**-0.3 * (m / 0.08)**-1.3,
                 0.08**-0.3 * (0.5 / 0.08)**-1.3 * (m / 0.5)**-2.3))


# def kroupa2002(x):
#     """
#     Kroupa (2002) Salpeter (1995) piecewise IMF taken from MASSCLEAN
#     article, Eq. (2) & (3), p. 1725
#     """
#     alpha = [-0.3, -1.3, -2.3]
#     m0, m1, m2 = [0.01, 0.08, 0.5]
#     factor = [(1. / m1) ** alpha[0], (1. / m1) ** alpha[1],
#               ((m2 / m1) ** alpha[1]) * ((1. / m2) ** alpha[2])]
#     if m0 < x <= m1:
#         i = 0
#     elif m1 < x <= m2:
#         i = 1
#     elif m2 < x:
#         i = 2
#     return factor[i] * (x ** alpha[i])

# def kroupa1993(x):
#     """
#     Kroupa, Tout & Gilmore. (1993) piecewise IMF.
#     http://adsabs.harvard.edu/abs/1993MNRAS.262..545K
#     Eq. (13), p. 572 (28)
#     """
#     alpha = [-1.3, -2.2, -2.7]
#     m0, m1, m2 = [0.08, 0.5, 1.]
#     factor = [0.035, 0.019, 0.019]
#     if m0 < x <= m1:
#         i = 0
#     elif m1 < x <= m2:
#         i = 1
#     elif m2 < x:
#         i = 2
#     return factor[i] * (x ** alpha[i])


# def chabrier2001_log(x):
#     """
#     Chabrier (2001) lognormal form of the IMF.
#     http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
#     Eq (7)
#     """
#     imf_vals = (1. / (np.log(10) * x)) * 0.141 * \
#         np.exp(-((np.log10(x) - np.log10(0.1)) ** 2) / (2 * 0.627 ** 2))

#     # Normalize PDF
#     int_imf = np.trapz(imf_vals, x)
#     imf_vals /= int_imf
#     return imf_vals


# def chabrier2001_exp(x):
#     """
#     Chabrier (2001) exponential form of the IMF.
#     http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
#     Eq (8)
#     """
#     return 3. * x ** (-3.3) * np.exp(-(716.4 / x) ** 0.25)


# def salpeter1955(x):
#     """
#     Salpeter (1955)  IMF.
#     https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
#     """
#     return x ** -2.35

def invTrnsfSmpl(masses, IMF_name, m_low=0.08, m_high=150):
    """
    IMF inverse transform sampling.

    Asked here: https://stackoverflow.com/q/21100716/1391441
    """
    imf_id = IMF_name.replace(' ', '_')

    # IMF mass interpolation step and grid values.
    mass_step = 0.05
    mass_values = np.arange(m_low, m_high, mass_step)

    #
    # *** USE THIS BLOCK TO RE_GENERATE THE CDFs ***
    #
    # # Generate and write CDFs
    # for IMF_name in ('Salpeter (55)', 'Miller-Scalo (79)', 'Kroupa (01)',
    #                  'Chabrier (03, indiv)', 'Chabrier (03, system)'):
    #     print(IMF_name)
    #     if IMF_name == 'Salpeter (55)':
    #         IMF_func = salpeter55
    #     elif IMF_name == 'Miller-Scalo (79)':
    #         IMF_func = millerscalo79
    #     elif IMF_name == 'Kroupa (01)':
    #         IMF_func = kroupa01
    #     elif IMF_name == 'Chabrier (03, indiv)':
    #         IMF_func = chabrier03individual
    #     elif IMF_name == 'Chabrier (03, system)':
    #         IMF_func = chabrier03system
    #     imf_id = IMF_name.replace(' ', '_')
    #     # Sample the CDF
    #     CDF_samples = []
    #     for m in mass_values:
    #         CDF_samples.append(quad(IMF_func, m_low, m)[0])
    #     # Normalize values
    #     CDF_samples = np.array(CDF_samples) / max(CDF_samples)
    #     with open('modules/{}.pickle'.format(imf_id), 'wb') as ff:
    #         pickle.dump(CDF_samples, ff, protocol=pickle.HIGHEST_PROTOCOL)
    # breakpoint()

    # Read CDFs from file
    with open('modules/{}.pickle'.format(imf_id), 'rb') as ff:
        CDF_samples = pickle.load(ff)

    CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()

    # Inverse CDF
    inv_cdf = interp1d(CDF_samples, mass_values)

    def sampled_inv_cdf(N):
        mr = np.random.rand(N)
        mr = mr[(mr >= CDF_min) & (mr <= CDF_max)]
        return inv_cdf(mr)

    # Sample in chunks of 100 stars until the maximum defined mass is reached.
    M_total, M_min, M_max = masses.sum(), masses.min(), masses.max()
    mass_IMF = []
    while np.sum(mass_IMF) < M_total:
        mass_samples = sampled_inv_cdf(100)
        msk = (mass_samples >= M_min) & (mass_samples <= M_max)
        if msk.sum() > 0:
            mass_IMF += mass_samples[msk].tolist()
    mass_IMF = np.array(mass_IMF)
    # print(mass_IMF.sum())

    i = np.argmin(abs(np.cumsum(mass_IMF) - M_total))
    mass_IMF = mass_IMF[:i + 1] * 1

    # Remove stars until the M_total value is approximated
    # i = 1
    # while mass_IMF.sum() > M_total:
    #     mass_IMF = mass_IMF[:-i]
    #     i += 5
    # print(mass_IMF.sum())
    # print(mass_IMF2.sum())

    return mass_IMF
