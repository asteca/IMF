
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
    masses_type, Nruns, alpha_min, alpha_max, mag_min, mag_max, mass_min,\
        mass_max, binar_cut = ipars.get('masses_type'),\
        ipars.getint('Nruns'),\
        ipars.getfloat('alpha_min'), ipars.getfloat('alpha_max'),\
        ipars.getfloat('mag_min'), ipars.getfloat('mag_max'),\
        ipars.getfloat('mass_min'), ipars.getfloat('mass_max'),\
        ipars.getfloat('binar_cut')

    if masses_type not in ('single', 'binar'):
        raise ValueError("The 'masses_type' value is not valid")

    return masses_type, Nruns, alpha_min, alpha_max, mag_min, mag_max,\
        mass_min, mass_max, binar_cut


def dataSave(data, file_out):
    """
    Create output data file
    """
    ascii.write(data, file_out, overwrite=True)


def dataRead(masses_type, file_in):
    """
    Read data file
    """
    data = ascii.read(file_in)

    if masses_type == 'single':
        mass, mass_std = data['M1'], data['M1_std']
    elif masses_type == 'binar':
        mass, mass_std = data['M2'], data['M2_std']

    binar_probs = data['P_binar']
    phot = np.array([data['Col'], data['Mag']]).T

    return mass, mass_std, binar_probs, phot


def readFiles():
    """
    Read files from the input folder
    """
    files = []
    for pp in Path('input').iterdir():
        if pp.is_file():
            files += [pp]
        # else:
        #     files += [arch for arch in pp.iterdir()]

    return files
