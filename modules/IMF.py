
import numpy as np


"""
Define any number of IMFs.

The package https://github.com/keflavich/imf has some more (I think,
24-09-2019).
"""


def kroupa1993(x):
    """
    Kroupa, Tout & Gilmore. (1993) piecewise IMF.
    http://adsabs.harvard.edu/abs/1993MNRAS.262..545K
    Eq. (13), p. 572 (28)
    """
    alpha = [-1.3, -2.2, -2.7]
    m0, m1, m2 = [0.08, 0.5, 1.]
    factor = [0.035, 0.019, 0.019]
    if m0 < x <= m1:
        i = 0
    elif m1 < x <= m2:
        i = 1
    elif m2 < x:
        i = 2
    return factor[i] * (x ** alpha[i])


def kroupa2002(x):
    """
    Kroupa (2002) Salpeter (1995) piecewise IMF taken from MASSCLEAN
    article, Eq. (2) & (3), p. 1725
    """
    alpha = [-0.3, -1.3, -2.3]
    m0, m1, m2 = [0.01, 0.08, 0.5]
    factor = [(1. / m1) ** alpha[0], (1. / m1) ** alpha[1],
              ((m2 / m1) ** alpha[1]) * ((1. / m2) ** alpha[2])]
    if m0 < x <= m1:
        i = 0
    elif m1 < x <= m2:
        i = 1
    elif m2 < x:
        i = 2
    return factor[i] * (x ** alpha[i])


def chabrier2001_log(x):
    """
    Chabrier (2001) lognormal form of the IMF.
    http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
    Eq (7)
    """
    imf_vals = (1. / (np.log(10) * x)) * 0.141 * \
        np.exp(-((np.log10(x) - np.log10(0.1)) ** 2) / (2 * 0.627 ** 2))

    # Normalize PDF
    int_imf = np.trapz(imf_vals, x)
    imf_vals /= int_imf
    return imf_vals


def chabrier2001_exp(x):
    """
    Chabrier (2001) exponential form of the IMF.
    http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
    Eq (8)
    """
    return 3. * x ** (-3.3) * np.exp(-(716.4 / x) ** 0.25)


def salpeter1955(x):
    """
    Salpeter (1955)  IMF.
    https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
    """
    return x ** -2.35
