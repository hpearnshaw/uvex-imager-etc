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
                # Source is provided as a quantity - treat as a flat spectrum
                # TODO: Add capacity to generate a range of spectrum types for given magnitudes
                if len(source) == 1:
                    self.source = [SourceSpectrum(ConstFlux1D, amplitude=source)]
                    self.n_source = 1
                else:
                    self.source = [SourceSpectrum(ConstFlux1D, amplitude=s) for s in source]
                    self.n_source = len(source)
                self.source_info = f'Constant spectrum at {source:.2f}'
            elif isinstance(source, SourceSpectrum):
                # Directly assign the spectrum
                self.source = [source]
                self.n_source = 1
                self.source_info = 'User-defined spectrum'
            elif isinstance(source, list) | isinstance(source, np.ndarray):
                if isinstance(source[0], SourceSpectrum):
                    # Directly assign list/array of spectra
                    self.source = source
                    self.n_source = len(source)
                self.source_info = 'User-defined spectra'
            else:
                raise ValueError("Source must be a flux Quantity or synphot SourceSpectrum (or list thereof)")
        else:
            # No need to define a source for limiting magnitude calculations
            self.source = None
            self.source_info = 'None'
        
        # Set source locations (used for calculating background)
        if coordinate is not None:
            if not isinstance(coordinate, SkyCoord):
                raise ValueError("Coordinate must be a `SkyCoord` object.")
            if (coordinate.size > 1) & (self.n_source > 1):
                if coordinate.size != self.n_source:
                    raise ValueError("Length of coordinate must match number of sources.")
            self.coord = coordinate
            self.n_coord = coordinate.size
        else:
            # Default 'average' location 15-deg out of Galactic Plane
            self.coord = SkyCoord(120., 15., unit=u.deg, frame='galactic')
            self.n_coord = 1
        
        # Set the observation times (used for calculating background)
        if obstime is not None:
            if not isinstance(obstime, Time):
                raise ValueError("Obstime must be a `Time` object.")
            elif obstime.size > 1 and (obstime.size != self.n_coord):
                raise ValueError("Length of obstime must be 1 or equal to length of coordinate.")
            elif obstime.size == self.n_coord:
                self.obstime = obstime
            elif obstime.size == 1:
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
        # Dark current, sky components, Cherenkov, scattered light?

    # Functions
    def get_info(self):
        '''
            Returns current information about ETC setup
        '''
        print(f'UVEX version: {self.telescope.get_caldb()}')
        if self.n_source > 1: n_s = f' x {self.n_source}'
        else: n_s = ''
        print(f'Source: {self.source_info}{n_s}')
        print(f'Source position: {self.coord}')
        print(f'Observation time: {self.obstime}')
    
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
        
    def _req_source(self, k, exposure, bgd_rate, read_noise, neff):
        """
        Isolate source flux to get at least SNR of k in exposure seconds

        Parameters
        -----------
        k : float
            Desired SNR
        exposure: float
            Exposure in seconds
        bgd_rate : float
            Combined sky and dark current
        read_noise : float
            Read noise per pixel
        neff : float
            Effective number of pixels
        """
        c = neff * k**2 * (read_noise**2 + exposure*(bgd_rate))
        source =  (k**2 + np.sqrt(k**4 + 4*c))/ (2*exposure)
        return source * u.ct / u.s
        
        
    def _calc_exposure(self, k, src_rate, bgd_rate, read_noise, neff):
        """
        Compute the time to get to a given significance (k) given the source rate,
        the background rate, the read noise, and the number
        of effective background pixels. Inversion of the standard CCD SNR equation.
        
        Parameters
        -----------
        k : float
            Desired SNR
        src_rate : float
            Source count rate
        bgd_rate : float
            Combined sky and dark current
        read_noise : float
            Read noise per pixel
        neff : float
            Effective number of pixels
        """
        denom = 2 * src_rate**2

        nom1 = (k**2) * (src_rate + neff*bgd_rate)
        nom2 = ( k**4 *(src_rate + neff*bgd_rate)**2 +
                        4 * k**2 * src_rate**2 * neff * read_noise**2)**(0.5)
        exposure = (nom1 + nom2) / denom
        return exposure * u.s

    def get_snr(self, exptime=None, n_frames=None, n_dwells=None, band='nuv'):
        """
        Calculate the SNR of an observation of a source with UVEX.

        Parameters
        ----------
        exptime : Quantity
            Exposure time
        
        n_frames : int
            Number of exposures added together
            
        n_dwells : int
            Sets exptime and n_frames to a specific number of dwells
            exptime and n_frame inputs are ignored in this case
            Dwells are defined using ETC properties
        
        band : 'nuv' or 'fuv'
            The UVEX band in which to calculate SNR
        
        Returns
        -------
        float array
            The signal to noise ratio
        """
        # Determine inputs
        band = band.lower()
        if not ((band == 'nuv') | (band == 'fuv')):
            raise ValueError(f"band must be 'nuv' or 'fuv'; got {band}")
        
        if n_dwells is not None:
            if band == 'nuv':
                exptime = self.nuv_exposure
                n_frames = self.n_nuv * n_dwells
            elif band == 'fuv':
                exptime = self.fuv_exposure
                n_frames = self.n_fuv * n_dwells
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
            gain=1.)
        snr *= np.sqrt(n_frames)
    
        return snr
    
    
    def get_limiting_mag(self, snr=5., exptime=None, n_frames=None, n_dwells=None, band='nuv'):
        """
        Get the limiting magnitude at a certain location and time for given SNR and exposure
        
        Does not require any source information to be loaded

        Parameters
        ----------
        snr : float
            Desired signal-to-noise ratio
        
        exptime : Quantity
            Exposure time
        
        n_frames : int
            Number of exposures added together
            
        n_dwells : int
            Sets exptime and n_frames to a specific number of dwells
            exptime and n_frame inputs are ignored in this case
            Dwells are defined using ETC properties
        
        band : 'nuv' or 'fuv'
            The UVEX band in which to calculate limiting magnitude
        
        Returns
        -------
        float array
            The limiting magnitude for each position
        """
        # Determine inputs
        band = band.lower()
        if not ((band == 'nuv') | (band == 'fuv')):
            raise ValueError(f"band must be 'nuv' or 'fuv'; got {band}")
        
        if n_dwells is not None:
            if band == 'nuv':
                exptime = self.nuv_exposure
                n_frames = self.n_nuv * n_dwells
                bandpass = self.telescope.nuv_bandpass
            elif band == 'fuv':
                exptime = self.fuv_exposure
                n_frames = self.n_fuv * n_dwells
                bandpass = self.telescope.fuv_bandpass
        else:
            if not isinstance(exptime, u.quantity.Quantity):
                raise ValueError("Exptime must be a Quantity.")
            if not isinstance(n_frames, int):
                raise ValueError("n_frames must be a positive integer.")
        
        # Load appropriate read noise and dark current from telescope
        dark_current = self.telescope.DARK_CURRENT[band].value
        read_noise = self.telescope.READ_NOISE[band].value
        npix = self.telescope.NPIX
        
        # Trigger generation of count rates if necessary
        if band not in self.background_count_rate: self._calc_background_count_rate()
        
        # Get reference count rate
        m_ref = 22*u.ABmag
        sp = SourceSpectrum(ConstFlux1D, amplitude=m_ref)
        obs_band = Observation(sp, bandpass)
        ref_rate = obs_band.countrate(area=self.telescope.AREA)
        
        # Get the required source rate per exposure
        per_exp_snr = snr/np.sqrt(n_frames)
        req_rate = self._req_source(per_exp_snr, exptime.to(u.s).value,
                                    self.background_count_rate[band].value + dark_current,
                                    read_noise, npix)
        ratio = req_rate / ref_rate
        m_limit = m_ref - (2.5*np.log10(ratio))*u.mag
        
        return m_limit
    
    
    def get_exposure(self, snr=5., band='nuv'):
        """
        Get the required exposure time to detect point source to a given SNR
        
        Returns both an exposure time for a single observation and the number of
        standard dwells required to stack to reach the required SNR

        Parameters
        ----------
        snr : float
            Desired signal-to-noise ratio
        
        band : 'nuv' or 'fuv'
            The UVEX band in which to calculate exposure time/dwells
        
        Returns
        -------
        exptime : float array
            The required exposure time in seconds for each source
        
        n_dwells : int
            The required number of standard dwells for each source
        """
        # Determine inputs
        band = band.lower()
        if not ((band == 'nuv') | (band == 'fuv')):
            raise ValueError(f"band must be 'nuv' or 'fuv'; got {band}")
        if band == 'nuv':
            exposure = self.nuv_exposure
        elif band == 'fuv':
            exposure = self.fuv_exposure
        
        # Load appropriate read noise and dark current from telescope
        dark_current = self.telescope.DARK_CURRENT[band].value
        read_noise = self.telescope.READ_NOISE[band].value
        npix = self.telescope.NPIX
        
        # Trigger generation of count rates if necessary
        if band not in self.source_count_rate: self._calc_source_count_rate()
        if band not in self.background_count_rate: self._calc_background_count_rate()
        
        # Get the required single-exposure time
        exptime = self._calc_exposure(snr, self.source_count_rate[band].value,
                                      self.background_count_rate[band].value + dark_current,
                                      read_noise, npix)
        
        # Get the required number of standard dwells
        snr_per_frame = self.get_snr(exposure, n_frames=1, band=band)
        n_frames = np.ceil((snr / snr_per_frame)**2)
        
        if band == 'FUV':
            n_dwells = n_frames
        else:
            n_dwells = np.ceil(n_frames / 3)
        
        return exptime, n_dwells



    
    # TODO: more ETC calculations
    # Get exposure time and/or number of dwells (for given SNR)
    
    # TODO: setters for ETC setup (reset source/bg rates as required)
    # Set source (reset source and background count rates)
    # Set obstime (reset background count rates)
    # Set coords (reset background count rates)
    # Set telescope (reset source and background count rates)
