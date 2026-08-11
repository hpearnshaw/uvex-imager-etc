import os
import json
import numpy as np
from numpy import pi as PI

import astropy.units as u

from synphot import Empirical1D, SpectralElement
from synphot.models import ConstFlux1D

from . import backgrounds

response_files_dir = os.path.join(os.path.dirname(__file__), 'response_files')

class UVEX():
    '''
        Container class to load telescope CBE values and bandpasses
        from uvex_response generated config
        
        Parameters
        ----------
        caldb : string
            Directory name for a UVEX CALDB version
    '''
    def __init__(self, caldb=None):
    
        # Check for existence of CALDBs
        avail_caldb = [f for f in os.listdir(response_files_dir) if not f.startswith('.')]
        if len(avail_caldb) == 0:
            raise ValueError("No available CALDBs in response_files.")
        
        # Define CALDB we are using
        if caldb is None:
            # TODO: code to get latest version
            self.caldb = '20260723_v0.1a'
        else:
            if caldb not in avail_caldb:
                raise ValueError(f"CALDB version {caldb} not available in response_files.")
            self.caldb = caldb
        config_file = os.path.join(response_files_dir,self.caldb,'config','response_files.json')
        with open(config_file) as f:
            self.config = json.load(f)
        caldb_dir = os.path.join(response_files_dir,self.caldb,'etc','imager')
    
        # Characteristics
        self.EPD = float(self.config['uvex']['EPD']['value']) * u.cm
        self.PIX_UM = float(self.config['uvex']['PIXEL_UM']['value']) * u.micron
        self.PLATE_SCALE = float(self.config['uvex']['PLATE_SCALE']['value']) * u.arcsec / self.PIX_UM
        self.PIXEL = (self.PIX_UM * self.PLATE_SCALE)**2 # Pixel area

        self.READ_NOISE = {
            'nuv': self.config['read_noise']['nuv']['value'] * u.electron,
            'fuv': self.config['read_noise']['fuv']['value'] * u.electron,
        }
        self.DARK_CURRENT = {
            'nuv': self.config['dark_current']['nuv']['value'] * u.electron / u.s,
            'fuv': self.config['dark_current']['fuv']['value'] * u.electron / u.s,
        }

        self.NPIX = 10.15 # Not currently in uvex_response
        
        # Ly-alpha background level in kR
        self.lya_kr = 2
        
        # Load in filter bandpasses from uvex_response generated config as SpectralElement objects
        nuv_data = np.genfromtxt(os.path.join(caldb_dir, 'nuv_bandpass.txt'))
        self.nuv_bandpass = SpectralElement(Empirical1D,
                                            points = nuv_data[:,0], lookup_table = nuv_data[:,1])
        fuv_data = np.genfromtxt(os.path.join(caldb_dir, 'fuv_bandpass.txt'))
        self.fuv_bandpass = SpectralElement(Empirical1D,
                                            points = fuv_data[:,0], lookup_table = fuv_data[:,1])
        
        # Placeholder values
        # TODO: get actual Cherenkov bandpasses from uvex-response
        self.nuv_cherenkov_bandpass = SpectralElement(Empirical1D,
                                            points = nuv_data[:,0], lookup_table = nuv_data[:,1]/100.)
        self.fuv_cherenkov_bandpass = SpectralElement(Empirical1D,
                                            points = fuv_data[:,0], lookup_table = fuv_data[:,1]/100.)
        # Other background-related bandpasses as needed

    @property
    def AREA(self):
        return PI * (self.EPD*0.5)**2
