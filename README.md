Ultraviolet Explorer (UVEX) Imager Exposure Time Calculator
-----------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

Overview
--------

This is the exposure time calculator for the UVEX imagers. Basic usage is demonstrated in notebooks/example.ipynb. It takes input sources, positions, and observation times, and performs various sensitivity calculations, such as required exposure times, limiting magnitudes, and signal-to-noise ratios.

Telescope configuration is loaded from uvex_response.

Dependencies
------------

We require at least python 3.12 for this installation. We'll assume for the installation 
purposes below that you have a miniconda installation of python.

Installation
------------

We recommend checking out this repo via git to ensure that you have the most recent
version of the code.

Make a local git folder if you don't have one, or navigate to where you want to check out
your git repo:

```
> mkdir git
> cd git
> git clone https://github.com/uvex-mission/uvex-imager-etc.git 
```

...which will download the files into the uvex-imager-etc directory.

```
> cd uvex-imager-etc
```

You'll want to create a standalone python environment for the UVEX ETC. We'll assume here
that you're using conda.

```
> conda create --name uvex-etc python=3.12
> conda activate uvex-etc
```

Once this is complete, install the UVEX ETC:

```
> pip install .
```

This should make uvex_imager_etc importable anywhere.

This UVEX ETC does not come pre-installed with the latest UVEX response curves as these are managed separately in uvex_response. To finish setting up the ETC, download the latest UVEX CALDB from [TBD] and place the unzipped, dated folder in uvex_imager_etc/response_files.

Example Usage
-------------

Exposure time calculator examples are found in uvex-imager-etc/notebooks/example.ipynb.
