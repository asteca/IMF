
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
    M_hist = in_params['Mass histogram parameters']
    Nruns, min_mass, max_mass = M_hist.getint('Nruns'),\
        M_hist.getfloat('min_mass'), M_hist.getfloat('max_mass')

    mcmc_p = in_params['MCMC parameters']
    mcmc_mode, burn_frac, nsteps, nwalkers = mcmc_p.get('mode'),\
        mcmc_p.getfloat('burn_frac'), mcmc_p.getint('nsteps'),\
        mcmc_p.getint('nwalkers')

    if mcmc_mode not in ("emcee", "pymc3"):
        raise ValueError("Mode not recognized: {}".format(mcmc_mode))

    return Nruns, min_mass, max_mass, mcmc_mode, burn_frac, nsteps, nwalkers


def dataSave(data, file_out):
    """
    Create output data file
    """
    ascii.write(data, file_out, overwrite=True)


def dataRead(file_in):
    """
    Read data file
    """
    return ascii.read(file_in)


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
