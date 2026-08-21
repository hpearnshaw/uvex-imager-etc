"""
Tests for uvex_imager_etc.backgrounds

These exercise the individual background components (Lyman-alpha, Galactic
diffuse, Zodiacal, Cherenkov) as well as the combined make_nuv_background /
make_fuv_background entry points.

Because these functions perform real synphot/astropy radiometric
calculations, most assertions here check units, shapes, finiteness, and
*relative* behavior (e.g. "doubling kr doubles the peak flux") rather than
hard-coded reference numbers, so the suite doesn't need pre-computed
"golden" values to stay meaningful.
"""
import numpy as np
import pytest
import astropy.units as u
from synphot import SourceSpectrum

from uvex_imager_etc import backgrounds


class TestLymanAlpha:
    def test_returns_source_spectrum(self, telescope):
        spec = backgrounds.make_lyman_spec(telescope)
        assert isinstance(spec, SourceSpectrum)

    def test_default_kr_matches_telescope_lya_kr(self, telescope):
        default = backgrounds.make_lyman_spec(telescope)
        explicit = backgrounds.make_lyman_spec(telescope, kr=telescope.lya_kr)
        assert default(1216 * u.AA) == explicit(1216 * u.AA)

    def test_flux_scales_linearly_with_kr(self, telescope):
        low = backgrounds.make_lyman_spec(telescope, kr=2)
        high = backgrounds.make_lyman_spec(telescope, kr=4)
        ratio = (high(1216 * u.AA) / low(1216 * u.AA)).decompose().value
        assert ratio == pytest.approx(2.0, rel=1e-2)

    def test_peaked_at_lyman_alpha(self, telescope):
        spec = backgrounds.make_lyman_spec(telescope)
        on_line = spec(1216 * u.AA)
        off_line = spec(1400 * u.AA)
        assert on_line > off_line


class TestGalacticSpec:
    def test_warns_for_in_plane_latitude(self, telescope):
        lat = [5.0, 30.0] * u.deg
        with pytest.warns(UserWarning, match="Galactic background invalid"):
            backgrounds.make_galactic_spec(telescope, lat, "nuv")

    def test_no_warning_outside_the_plane(self, telescope, recwarn):
        lat = [20.0, -20.0] * u.deg
        backgrounds.make_galactic_spec(telescope, lat, "nuv")
        assert not any("Galactic background invalid" in str(w.message) for w in recwarn.list)

    def test_boundary_latitude_of_15deg_does_not_warn(self, telescope, recwarn):
        # abs(lat) < 15 deg triggers the warning, so exactly 15 deg should not.
        lat = [15.0] * u.deg
        backgrounds.make_galactic_spec(telescope, lat, "nuv")
        assert not any("Galactic background invalid" in str(w.message) for w in recwarn.list)

    @pytest.mark.parametrize("band", ["nuv", "fuv"])
    def test_returns_one_spectrum_per_latitude(self, telescope, band):
        lat = [20.0, -20.0, 40.0] * u.deg
        spectra = backgrounds.make_galactic_spec(telescope, lat, band)
        assert len(spectra) == len(lat)
        assert all(isinstance(s, SourceSpectrum) for s in spectra)

    def test_single_latitude_still_returns_a_list(self, telescope):
        spectra = backgrounds.make_galactic_spec(telescope, [30.0] * u.deg, "nuv")
        assert isinstance(spectra, list)
        assert len(spectra) == 1


class TestZodiacalSpec:
    def test_load_zodi_spatial_returns_finite_callable(self, telescope):
        model = backgrounds.load_zodi_spatial()
        result = model(np.array([10.0, 50.0]), np.array([30.0, 90.0]))
        assert np.all(np.isfinite(result))

    def test_zodi_spec_returns_one_spectrum_per_scale_value(self, telescope):
        scale = np.array([77.0, 150.0, 300.0])
        spectra = backgrounds.zodi_spec(telescope, scale=scale)
        assert len(spectra) == len(scale)
        assert all(isinstance(s, SourceSpectrum) for s in spectra)

    def test_zodi_spec_scales_with_input_scale(self, telescope):
        low = backgrounds.zodi_spec(telescope, scale=np.array([77.0]))[0]
        high = backgrounds.zodi_spec(telescope, scale=np.array([154.0]))[0]
        wave = 5000 * u.AA
        ratio = (high(wave) / low(wave)).decompose().value
        assert ratio == pytest.approx(2.0, rel=1e-2)

class TestCherenkovSpectrum:
    def test_output_shape_matches_wavelength_grid(self, telescope):
        wave = np.arange(1000, 2000) * u.AA
        spec = backgrounds.gen_cherenkov_spectrum(telescope, wave)
        assert len(spec) == len(wave)

    def test_output_is_non_negative(self, telescope):
        wave = np.arange(1000, 2000) * u.AA
        spec = backgrounds.gen_cherenkov_spectrum(telescope, wave)
        assert np.all(spec.value >= 0)


class TestCombinedBackgrounds:
    def test_nuv_background_single_pointing(self, telescope, default_coord, default_obstime):
        rate = backgrounds.make_nuv_background(telescope, default_coord, default_obstime)
        assert rate.unit.is_equivalent(u.electron / u.s)
        assert rate.size == 1
        assert np.all(np.isfinite(rate.value))
        assert np.all(rate.value > 0)

    def test_fuv_background_single_pointing(self, telescope, default_coord, default_obstime):
        rate = backgrounds.make_fuv_background(telescope, default_coord, default_obstime)
        assert rate.unit.is_equivalent(u.electron / u.s)
        assert rate.size == 1
        assert np.all(np.isfinite(rate.value))
        assert np.all(rate.value > 0)

    def test_nuv_background_multi_pointing_matches_coord_length(
        self, telescope, multi_coord, default_obstime
    ):
        rate = backgrounds.make_nuv_background(telescope, multi_coord, default_obstime)
        assert rate.size == len(multi_coord)
        assert np.all(np.isfinite(rate.value))

    def test_nuv_background_multi_time_matches_obstime_length(
        self, telescope, default_coord, multi_obstime
    ):
        rate = backgrounds.make_nuv_background(telescope, default_coord, multi_obstime)
        assert rate.size == len(multi_obstime)
        assert np.all(np.isfinite(rate.value))
