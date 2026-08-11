import numpy as np

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.stats import signal_to_noise_oir_ccd

from synphot import SourceSpectrum, Empirical1D, SpectralElement, Observation
from synphot.models import ConstFlux1D

from . import uvex
from . import backgrounds

class ETC():
    '''
        Class to hold ETC-related properties and operations
        
        Parameters
        ----------
        coordinate : SkyCoord
            Source coordinates as SkyCoord object
            Defaults to an 'average' location 15-deg out of Galactic Plane
        
        obstime : Time
            Time of observation for each source
        
        source : Quantity or SourceSpectrum
            Source magnitude/flux as Quantity or a synphot SourceSpectrum object
        
        TODO: add a source_type option to generate a specific spectral shape
        if source is provided as a magnitude (for now just a constant)
        
        telescope : UVEX
            UVEX object containing a particular telescope configuration
        
        Attributes
        ----------
        nuv_exposure : default exposure duration
        fuv_exposure : default exposure duration
        n_nuv : number of NUV frames in a dwell
        n_fuv : number of FUV frames in a dwell
        coord : SkyCoord object with coordinate list
        n_coord : number of coordinates
        obstime : Time object of observation time list
        source : List of SourceSpectrum objects
    '''
    def __init__(self, source=None, coordinate=None, obstime=None, telescope=None):
        # Standard observing dwell definition
        self.nuv_exposure = 300*u.s
        self.fuv_exposure = 900*u.s
        self.n_nuv = 3
        self.n_fuv = 1
        
        # Ingest sources
        if source is not None:
            if isinstance(source, u.quantity.Quantity):
                # TODO: Add capacity to generate a range of spectrum types for given magnitudes
                if len(source) == 1:
                    self.source = [SourceSpectrum(ConstFlux1D, amplitude=source)]
                    self.n_source = 1
                else:
                    self.source = [SourceSpectrum(ConstFlux1D, amplitude=s) for s in source]
                    self.n_source = len(source)
            elif isinstance(source, SourceSpectrum):
                # Directly assign the spectrum
                self.source = [source]
                self.n_source = 1
            elif isinstance(source, list) | isinstance(source, np.ndarray):
                if isinstance(source[0], SourceSpectrum):
                    # Directly assign list/array of spectra
                    self.source = source
                    self.n_source = len(source)
            else:
                raise ValueError("Source must be a flux Quantity or synphot SourceSpectrum (or list thereof)")
        else:
            # If a source isn't provided then do we want some default source?
            print('No source provided')
            self.source = None
        
        # Set source locations (used for calculating background)
        if coordinate is not None:
            if not isinstance(coordinate, SkyCoord):
                raise ValueError("Coordinate must be a `SkyCoord` object.")
            if (len(coordinate) > 1) & (self.n_source > 1):
                if len(coordinate) != self.n_source:
                    raise ValueError("Length of coordinate must match number of sources.")
            self.coord = coordinate
            self.n_coord = len(coordinate)
        else:
            # Default 'average' location 15-deg out of Galactic Plane
            self.coord = SkyCoord(120., 15., unit=u.deg, frame='galactic')
            self.n_coord = 1
        
        # Set the observation times (used for calculating background)
        if obstime is not None:
            if not isinstance(obstime, Time):
                raise ValueError("Obstime must be a `Time` object.")
            elif (len(obstime) > 1) and (len(obstime) != self.n_coord):
                raise ValueError("Length of obstime must be 1 or equal to length of coordinate.")
            elif len(obstime) == self.n_coord:
                self.obstime = obstime
            elif len(obstime) == 1:
                time = obstime * self.n_coord
                self.obstime = Time(time, scale='utc', format='iso')
        else:
            # Generate default obstimes based on number of coordinates
            time = ['2030-06-01 09:00:00'] * self.n_coord
            self.obstime = Time(time, scale='utc', format='iso')
        
        if telescope is not None:
            if not isinstance(telescope, uvex.UVEX):
                raise ValueError("Telescope must be a `UVEX` object.")
            self.telescope = telescope
        else:
            # Initialize a UVEX object with default parameters
            self.telescope = uvex.UVEX()
        
        # Initialize source and background count rates
        self.source_count_rate = {}
        self.background_count_rate = {}
        
        # TODO: Add functionality to switch certain background effects on and off
        # Dark current, sky components, Cherenkov, scattered light

    # Functions
    def _calc_source_count_rate(self):
        '''
            Calculate and set the count rate for all sources
        '''
        nuv_rate, fuv_rate = np.array([]), np.array([])
        for s in self.source:
            # TODO: Make this more efficient
            nuv_obs = Observation(s, self.telescope.nuv_bandpass)
            nuv_rate = np.append(nuv_rate, nuv_obs.countrate(area=self.telescope.AREA).value)
            fuv_obs = Observation(s, self.telescope.fuv_bandpass)
            fuv_rate = np.append(fuv_rate, fuv_obs.countrate(area=self.telescope.AREA).value)
        self.source_count_rate['nuv'] = nuv_rate * u.electron / u.s
        self.source_count_rate['fuv'] = fuv_rate * u.electron / u.s
    
    def _calc_background_count_rate(self):
        '''
            Calculate and set the background rate for all background locations and times
        '''
        # Calculate backgrounds
        self.background_count_rate['nuv'] = backgrounds.make_nuv_background(self.telescope, self.coord, self.obstime)
        self.background_count_rate['fuv'] = backgrounds.make_fuv_background(self.telescope, self.coord, self.obstime)

    def get_snr(self, exptime=None, n_frames=None, n_dwells=None, band='nuv'):
        """Calculate the SNR of an observation of a point source with UVEX.

        Parameters
        ----------
        exptime : Quantity
            Exposure time
        
        n_frames : int
            Number of exposures added together
            
        n_dwells : int
            Sets exptime and n_frames to a specific number of dwells
            exptime and n_frame inputs are ignored in this case
            Dwells are defined as 3 x 300s NUV exposures or 1 x 900s FUV exposure
        
        band : 'nuv' or 'fuv'
            The UVEX band in which to calculate SNR
        
        Returns
        -------
        float
            The signal to noise ratio
        """
        # Determine inputs
        band = band.lower()
        if not ((band == 'nuv') | (band == 'fuv')):
            raise ValueError(f"band must be 'nuv' or 'fuv'; got {band}")
        
        if n_dwells is not None:
            if band == 'nuv':
                exptime = 300 * u.s
                n_frames = 3
            elif band == 'fuv':
                exptime = 900 * u.s
                n_frames = 1
        else:
            if not isinstance(exptime, u.quantity.Quantity):
                raise ValueError("Exptime must be a Quantity.")
            if not isinstance(n_frames, int):
                raise ValueError("n_frames must be a positive integer.")
        
        # Load appropriate read noise and dark current from telescope
        dark_current = self.telescope.DARK_CURRENT[band]
        read_noise = self.telescope.READ_NOISE[band]
        npix = self.telescope.NPIX
        
        # Trigger generation of count rates if necessary
        if band not in self.source_count_rate: self._calc_source_count_rate()
        if band not in self.background_count_rate: self._calc_background_count_rate()
        
        source = self.source_count_rate[band]
        sky = self.background_count_rate[band]
    
        snr = signal_to_noise_oir_ccd(exptime.to(u.s).value,
            source.value,
            sky.value,
            dark_eps=dark_current.value,
            rd=read_noise.value,
            npix=npix,
            gain=1.
        )
        snr *= np.sqrt(n_frames)
    
        return snr
        
    
    # TODO: more ETC calculations
    # Get exposure time and/or number of dwells (for given SNR)
    # Get limiting mag (either exposure time * frames or n_dwells)
    
    # TODO: setters for ETC setup (reset source/bg rates as required)
    # Set source (reset source and background count rates)
    # Set obstime (reset background count rates)
    # Set coords (reset background count rates)
    # Set telescope (reset source and background count rates)
