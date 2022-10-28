# ixpeobssim
Simulation and analysis framework for the Imaging X-ray Polarimetry Explorer

This is the public repository for the package described at [SoftwareX, Volume 19, July 2022, 101194](https://www.sciencedirect.com/science/article/pii/S2352711022001169)

![Tag](https://img.shields.io/github/v/tag/lucabaldini/ixpeobssim)
[![PyPI version](https://badge.fury.io/py/ixpeobssim.svg)](https://badge.fury.io/py/ixpeobssim)
![PyPI - Downloads](https://img.shields.io/pypi/dw/ixpeobssim)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/ixpeobssim/badge/?version=latest)](https://ixpeobssim.readthedocs.io/en/latest/?badge=latest)

![GitHub repo size](https://img.shields.io/github/repo-size/lucabaldini/ixpeobssim)
![GitHub contributors](https://img.shields.io/github/contributors/lucabaldini/ixpeobssim)
![GitHub last commit](https://img.shields.io/github/last-commit/lucabaldini/ixpeobssim)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/lucabaldini/ixpeobssim)
![Libraries.io dependency status for GitHub repo](https://img.shields.io/librariesio/github/lucabaldini/ixpeobssim)

![GitHub pull requests](https://img.shields.io/github/issues-pr/lucabaldini/ixpeobssim)
![GitHub issues](https://img.shields.io/github/issues/lucabaldini/ixpeobssim)
![GitHub closed issues](https://img.shields.io/github/issues-closed/lucabaldini/ixpeobssim)

Please refer to the
[documentation](https://ixpeobssim.readthedocs.io/en/latest/?badge=latest)
and, should you encounter any problem, feel free to use our
[issue tracker](https://github.com/lucabaldini/ixpeobssim/issues).


# Synopsis

ixpeobssim is a simulation and analysis framework specifically developed for
the Imaging X-ray Polarimetry Explorer (IXPE). Given a source model and the
response functions of the telescopes, it is designed to produce realistic
simulated observations, in the form of event lists in FITS format, containing a
strict superset of the information included in the publicly released IXPE data
products. The core simulation capabilities are complemented by a full suite of
post-processing applications which support the spatial, spectral, and temporal
models needed for analysis of typical polarized X-ray sources, allowing for the
implementation of complex, polarization-aware analysis pipelines, and
facilitating the interoperation with the standard visualization and analysis
tools traditionally in use by the X-ray community.

Although much of the framework is specific to IXPE, the modular nature of the
underlying implementation makes it potentially straightforward to adapt it to
different missions with polarization capabilities.


# Installation

In a nutshell, you should be able to install `ixpeobssim` via pip

```
pip install ixpeobssim
```

Refer to the [installation instructions](https://ixpeobssim.readthedocs.io/en/latest/installation.html) for more in-depth information.


# Documentation

The ixpeobssim documentation is hosted on [readthedocs](https://ixpeobssim.readthedocs.io/en/latest/index.html), the most
useful pointers being:

* [overview](https://ixpeobssim.readthedocs.io/en/latest/overview.html)
* [quick start](https://ixpeobssim.readthedocs.io/en/latest/quickstart.html)
* [showcase](https://ixpeobssim.readthedocs.io/en/latest/showcase.html)
* [application reference](https://ixpeobssim.readthedocs.io/en/latest/reference.html)

For convenience, a static version of the latest documentation in pdf format is available at
[this link](https://ixpeobssim.readthedocs.io/_/downloads/en/latest/pdf/).
