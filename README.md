Ultraviolet Explorer (UVEX) Imager Exposure Time Calculator
-----------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge


Overview
--------

This is the exposure time calculator for the UVEX imagers.

Dependecies
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

Once this is complete, install a dev version of the UVEX ETC:

```
> pip install -e .
```

This should make uvex_imager_etc importable anywhere and you should be good to go.

ETCs and Notebooks
-------------------

Exposure time calculator examples are found in uvex-imager-etc/notebooks (coming soon)
