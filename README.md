# IMF

Estimate the IMF for an observed cluster, given its individual masses (and uncertainties). This is a rather complicated processes to generalize, so I moved it outside of ASteCA and into its own package.

Related literature:

- [CCD photometric and mass function study of nine young Large Magellanic Cloud star clusters, Kumar et al. (2008)](http://mnras.oxfordjournals.org/content/386/3/1380.full)
- [On the function describing the stellar initial mass function, Maschberger  (2012)](http://arxiv.org/abs/1212.0939)
- [Mass distribution and structural parameters of Small Magellanic Cloud star clusters, Maia et al. (2013)](http://adsabs.harvard.edu/abs/2013arXiv1310.5934M).
- [The Range of Variation of the Mass of the Most Massive Star in Stellar Clusters Derived from 35 Million Monte Carlo Simulations, Popescu & Hanson (2014)](http://adsabs.harvard.edu/abs/2014ApJ...780...27P)


## Installation

Install the requirements in a `conda` environment with:

```
$ conda create --name imfenv python=3.8 numpy matplotlib astropy
$ conda activate imfenv
(imfenv) $ conda install -c conda-forge emcee pymc3 arviz
```

## Analysis

1. The per-star masses are re-sampled using a Normal distribution, and a new histogram is generated each time. This histograms are what will be fed to the MCMC samplers.
2. The selected MCMC sampler is applied. The package performs two different analysis employing `emcee` and `PyMC3`:
  1. The `emcee` mode uses the method developed in [LinearBayes3](https://github.com/Gabriel-p/LinearBayes3) which accounts for outliers in the model. This method requires substantially more steps to properly converge, compared to the `pymc3` mode.
  2. The `pymc3` mode applies the [Robust Linear Regression](https://docs.pymc.io/notebooks/GLM-robust.html) method. The basis of the method is that "*minimizing the squared distance of the regression line is the same as maximizing the likelihood of a Normal distribution with the mean coming from the regression line*". The method employs a Students T distribution