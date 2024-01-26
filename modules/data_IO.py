
import numpy as np
import configparser
from pathlib import Path
from astropy.io import ascii


def readINI():
    """
    Read .ini config file
    """
    in_params = configparser.ConfigParser()

    ini_file = Path("params.not_tracked.ini")
    if not ini_file.is_file():
        ini_file = Path("params.ini")
    in_params.read(ini_file)

    # Data columns
    ipars = in_params['Input parameters']
    make_plot, mag_min, mag_max, binar_cut = ipars.getboolean('make_plot'), \
        ipars.getfloat('mag_min'), ipars.getfloat('mag_max'), \
        ipars.getfloat('binar_cut')

    return mag_min, mag_max, binar_cut, make_plot


def dataRead(file_in):
    """
    Read data file
    """
    data = ascii.read(file_in)
    mass, mass_std = data['M1'], data['M1_std']
    binar_probs = data['P_binar']
    phot = np.array([data['Col'], data['Mag']]).T

    return mass, mass_std, binar_probs, phot


def readFiles():
    """
    Read files from the input folder
    """
    path = "input"
    files = []
    for pp in Path(path).iterdir():
        if pp.is_file():
            files += [pp]

    return files
