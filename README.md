# IMF

Estimate the IMF for an observed cluster, given its individual masses (and uncertainties). This is a rather complicated processes to generalize, so I moved it outside of ASteCA and into its own package.

## Installation

Install the requirements in a `conda` environment with:

```
$ conda create --name imfenv python=3.8 numpy matplotlib seaborn astropy
$ conda activate imfenv
```

## Analysis

The method is originally developed in [Khalaj & Baumgardt (2013)](https://academic.oup.com/mnras/article/434/4/3236/960889) and used in [Sheikhi et al. (2016)](https://academic.oup.com/mnras/article/457/1/1028/989829)

The analysis requires the user to select first which stars will be processed. The `binar_min, binar_max` parameters filter stars by their probability of being a binary system. The `masses_type` parameter decides which masses will be used in the analysis: single or binary system masses.

After this the range for the slope is given by the `alpha_min, alpha_max` parameters, and the ranges for magnitude and mass by the `mag_min, mag_max` and `mass_min, mass_max` respectively.

Finally, the number of bootstrap runs is controlled by the `Nruns` parameter.


## Output

![Alt text](output.png?raw=true)

The plots in the output figure represent the following analysis:

**A**.
The histogram shows the distribution of binary probabilities (BPs) estimated by `ASteCA`. This is the probability that a given star in the CMD is a binary system. Hence: a low BP value means that the star is most likely a single system, and a large BP means it is most likely a binary system.

The region shaded light red shows the BPs of the stars that were selected for the subsequent analysis. This selection is performed by:

1. fixing the `binar_cut` parameter, which represent the probability value that separates single stars from binary systems. Stars with BPs below this value will be considered single stars, and binary systems for stars above this value.
2. fixing the `masses_type` parameter to either `single` or `binar`. Selecting the former means that the code will select stars with BPs *below* the `binar_cut` value, and consider all of them to be single systems. The latter on the other hand selects stars *above* the `binar_cut` value and considers all of them to be binary systems.

For example, if you want to analyze the IMF of the sequence of single systems, you'd select something like:

    binar_cut      = .5
    masses_type    = single

If you wanted to analyze the IMF of the most probable binary systems, you could select:

    binar_cut      = .75
    masses_type    = binar

The black curve shows the *binary fraction* of the cluster for different values of the `binar_cut` parameter. It has a value of 1 for `binar_cut=0` (since all stars are considered to be binaries) and 0 for `binar_cut=1` (since now all stars are considered to be single systems)

**B**.
Shows the probability of a star in a given magnitude range (top plot) or mass range (bottom plot) of being either a binary system (orange bar), or a single system (blue bar). Notice that both bars add up to 1 since the probability of being a binary is 1 minus the probability of being a single system (no systems of higher number are considered, i.e. ternary, etc.)

**C**.
CMD of the input data. The stars colored in red are those that are used in the analysis, selected as described in plot **A**. The rest (black stars) are ignored.

The user can also restrict the analysis to a magnitude and/or mass range using the parameters:

    mag_min        = 0
    mag_max        = 25
    mass_min       = 1.
    mass_max       = 100

It is recommended to use a `mass_min` value of 1, since below this mass value the IMF will most likely be a poor fit.

**D**.
Bootstrap distribution for the fitted slope (alpha parameter in the IMF)

**E**.
The dashed black line shows the maximum likelihood fit (as described in *Khalaj & Baumgardt 2013*) on the selected set of stars shown in plot **C** as red stars. The remaining colored dashed lines are the slopes of several theoretical IMFs, for comparison.

The colored crosses show the slope values that one would obtain applying the simpler histogram fit on the selected set of stars, for three different number of bins used (5, 10, 25).

The gray shaded region below `m=1` shows the region where most of the theoretical IMFs split their functions. The brown dotted line is the KDE of all the observed stars, for comparison.

**F**.
The black point to the left of the plot is the same slope value shown in plot ***E** for the likelihood analysis performed on the selected mass range. The rest of the colored crosses are the slopes obtained with the same analysis, but using different mass ranges (shown in the legend)

**G**.
This is the same likelihood analysis performed for our data in plot **E**, but applied instead to synthetic data sampled for the theoretical IMFs listed. All the IMFs are sampled in the same mass range as the observed data, such that the total mass is the same. The idea is to see what does the likelihood method return for data sampled from theoretical IMFs in this mass range.

The crosses represent the simple histogram fit for each theoretical data sample, using 25 bins.



--------------------------------------------------------------------
Related literature:

- [CCD photometric and mass function study of nine young Large Magellanic Cloud star clusters, Kumar et al. (2008)](http://mnras.oxfordjournals.org/content/386/3/1380.full)
- [On the function describing the stellar initial mass function, Maschberger  (2012)](http://arxiv.org/abs/1212.0939)
- [Mass distribution and structural parameters of Small Magellanic Cloud star clusters, Maia et al. (2013)](http://adsabs.harvard.edu/abs/2013arXiv1310.5934M).
- [The Range of Variation of the Mass of the Most Massive Star in Stellar Clusters Derived from 35 Million Monte Carlo Simulations, Popescu & Hanson (2014)](http://adsabs.harvard.edu/abs/2014ApJ...780...27P)