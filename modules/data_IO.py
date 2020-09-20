
import numpy as np
import configparser
from pathlib import Path
from astropy.io import ascii


def readINI():
    """
    Read .ini config file
    """

    # def vtype(var):
    #     tp, v = var.split('_')
    #     if tp == 'int':
    #         return int(v)
    #     elif tp == 'float':
    #         return float(v)
    #     elif tp == 'bool':
    #         return bool(strtobool(v))
    #     elif tp == 'str':
    #         return v

    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    # Data columns
    ipars = in_params['Input parameters']
    masses_type, Nruns, alpha_min, alpha_max, min_mag, max_mag, min_mass,\
        max_mass, binar_min, binar_max = ipars.get('masses_type'),\
        ipars.getint('Nruns'),\
        ipars.getfloat('alpha_min'), ipars.getfloat('alpha_max'),\
        ipars.getfloat('min_mag'), ipars.getfloat('max_mag'),\
        ipars.getfloat('min_mass'), ipars.getfloat('max_mass'),\
        ipars.getfloat('binar_min'), ipars.getfloat('binar_max')

    if masses_type not in ('single', 'binar'):
        raise ValueError("The 'masses_type' value is not valid")

    if binar_max <= binar_min:
        raise ValueError("The maximum binary fraction value must be larger "
                         "than the minimum value")

    return masses_type, Nruns, alpha_min, alpha_max, min_mag, max_mag,\
        min_mass, max_mass, binar_min, binar_max


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
        mass_mean, mass_std = data['Mass_mu'], data['Mass_std']
    elif masses_type == 'binar':
        mass_mean, mass_std = data['Mass_binar_mu'], data['Mass_binar_std']

    binar_probs = data['P_binar']
    phot = np.array([data['Col'], data['Mag']]).T

    return mass_mean, mass_std, binar_probs, phot


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
