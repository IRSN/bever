R Package bever
================
Yves Deville
2023-05-22

# Description

**bever** is a R package funded by IRSN/Behrig.

This package is devoted to some Bayesian approaches in Extreme Value
(EV) modelling, with focus on Markov Chain Monte-Carlo inference. The
package does not implement by itself the Bayesian estimation or the
inference of any EV model. These tasks are left to the **revdbayes**
package or to some other MCMC sampler as can be implemented by using
**JAGS** or **Stan**. The **bever** package aims at producing some
classical results or plots such as Return Level tables or plots. Its
main goal is to allow a comparison of the frequentist and Bayesian
inference methods when applied to EV models.

The name **bever** relates to **B**ayes **E**xtreme-**V**alue. In some
modern or ancient European languages, the word relates to the *river*
or/and to one of its famous inhabitant. Hydrology may be one field of
application.

## Important notes

This package is still in its early stage of development

- The use of block maxima with block duration differing from 1 is not
  well tested.

- For now it is possible to draw the predictive curve on the same
  graphics as the RL plot only for GEV models, not for Poisson-GP
  models. The reason is that some analysis is still required for the
  transformation of plotting positions in “Renext style” to use them in
  “block style”.

# Installation

## Using the *remotes* package

In an R session use

``` r
library(remotes)
install_github("IRSN/bever", dependencies = TRUE)
```

This should install the package and make it ready to use.

Mind that by default this does not build the vignette shipped with the
package (long-form documentation). To build the vignette, use instead

``` r
install_github("IRSN/bever", dependencies = TRUE, build_vignettes = TRUE)
```

The installation will then take a longer time but the vignette will be
accessible from the help of the package (link above the “Help Pages”
section).

You can also select a specific branch or a specific commit by using the
suitable syntax for `install_github`. For instance to install the branch
`develop` use

``` r
install_github("IRSN/bever@develop", dependencies = TRUE)
```

See the **remotes** package documentation for more details.

## Clone, build and install

### Cloning the repository

You can also clone the repository to install the package. If you do not
have yet a local `bever` repository, use `git clone` to clone the
`bever` repository

``` bash
git clone https://github.com/IRSN/bever
```

This will create a `bever` sub-directory of the current directory,
i.e. the directory from which the git command was issued. Of course this
can work only if you have the authorisation to clone.

### Installation

Move to the parent directory of your cloned repository and use the
following command from a terminal to create a tarball source file

``` bash
R CMD build bever
```

This will produce a source tarball `bever_x.y.z` where `x`, `y` and `z`
stand for the major, minor and patch version numbers. Then you can
install from a command line

``` bash
R CMD INSTALL bever_x.y.z
```

Note that you must also have all the packages required by **bever**
installed.

You can also use the **RStudio** IDE to install the package.
