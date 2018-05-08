"""Driver for command-line bosssky script
"""
from __future__ import division, print_function

import os.path

import numpy as np

import astropy.table
import astropy.io.fits as fits
import astropy.constants
import astropy.units as u

import wpca

import bossdata.bits

import bossdata

from six import binary_type


def get_sky(plate, mjd, output_path, verbose=False):
    """Extract individual sky exposures used for a single coadd.

    The results are saved to a single FITS file named sky-<plate>-<mjd>.fits
    with the following HDUs:

     - 0: header with cards PLATE, MJD, NFIBERS, NEXP
     - 1: table with per-exposure metadata.
     - 2: table with per-sky-fiber plugmap metadata.
     - remaining ImageHDUs named xWLEN, xWDISP, xRDNOISE, xFLAT, xFLUX, xIVAR,
       xMASK, where x=B,R identifies the band. Each HDU contains a 2D array
       with shape (NFIBERS * NEXP, NWLEN), with NWLEN=4112(B) or 4128(R).

    xFLUX is  in flat-field corrected electrons, with corresponding pipeline
    inverse variances xIVAR. xFLUX * xFLAT and xRDNOISE are in units of
    detected electrons.  xWLEN and xWDISP are in angstroms.

    Parameters
    ----------
    plate : int
        Plate number identifying the coadd.
    mjd : int
        MJD number identifying the coadd.
    output_path : str
        Path where output FITS file should be written.
    verbose : bool
        Print progress updates when True.

    Returns
    -------
    astropy.table.Table
        Table of metadata for the extracted exposures.
    """
    tag = f'PLATE {plate:05d} MJD {mjd:05d} PATH {output_path}'
    if verbose:
        print('Getting {}'.format(tag))
    # Initialize output data.
    last_nexp = None
    plugmaps = []
    wlens = {'b': [], 'r': []}
    wdisps = {'b': [], 'r': []}
    fluxes = {'b': [], 'r': []}
    ivars = {'b': [], 'r': []}
    flats = {'b': [], 'r': []}
    rdnoises = {'b': [], 'r': []}
    masks = {'b': [], 'r': []}
    obskeys = ('EXPOSURE', 'TAI-BEG', 'EXPTIME', 'AZ', 'ALT', 'AIRMASS',
               'PRESSURE', 'AIRTEMP',
               'RDNOISE0', 'RDNOISE1', 'RDNOISE2', 'RDNOISE3')
    obsvals = {key: [] for key in obskeys}
    # Size of each amplifier in raw image pixels along (wlen, tracex) axes.
    ampsize = {'b': (2056, 2048), 'r': (2064, 2057)}
    # ampx[band] tabulates whether each wavelength index is readout by
    # amplifier 0/2 (=0) or 1/3 (=1).
    ampx = {'b': 1 * (np.arange(4112) >= 2056),
            'r': 1 * (np.arange(4128) >= 2064)}
    # amplifer[band] is a function that takes a traceset as input an returns an
    # array that tabulates whether each wavelength index is readout by
    # amplifier 0-3.
    amplifier = {'b': lambda x: 2 * (x >= 2048) + ampx['b'],
                 'r': lambda x: 2 * (x >= 2057) + ampx['r']}
    # Scaling such that RMS = rdnoise_scale * RDNOISEn * neff.
    rdnoise_scale = (4 * np.pi) ** 0.25
    # Conversion from constant log-lambda pixels to wavelength ratio.
    wdisp_const = 1e-4 * np.log(10)
    # Allowed pixel mask bits.
    valid_mask = (1 << 32) - 1
    # Slices of valid data to save. These trim pixels at each end where
    # IVAR=0 or other serious pixel mask bits are often set.
    valid_slices = {'b': slice(767, 3299), 'r': slice(483, 3668)}
    # Initialize data access.
    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()
    # Loop over spectrographs.
    expected_fibers = []
    for specidx in 1, 2:
        # Load list of science exposures used for this spectrograph's coadd.
        fiber = 500 * (specidx - 1) + 1
        spec_name = finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
        exposures = bossdata.spec.SpecFile(mirror.get(spec_name)).exposures
        for band in 'b', 'r':
            camera = '{}{}'.format(band, specidx)
            use = valid_slices[band]
            # Loop over science exposures for this camera.
            nexp = exposures.num_by_camera[camera]
            if not (last_nexp is None or nexp == last_nexp):
                print(f'Different nexp for {camera} {tag}')
                return None
            last_nexp = nexp
            for expidx in range(nexp):
                # Load this camera's spFrame file.
                name = exposures.get_exposure_name(expidx, camera, 'spFrame')
                path = mirror.get(finder.get_plate_path(plate, name))
                spFrame = bossdata.plate.FrameFile(path, calibrated=False)
                # Lookup this spectrograph's sky fibers.
                sky_name = binary_type('SKY             ', 'ascii')
                fiberidx = np.where(
                    spFrame.plug_map['OBJTYPE'] == sky_name)[0]
                if expidx == 0 and band == 'b':
                    # Save plugmap metadata.
                    plugmaps.append(
                        spFrame.plug_map[
                            ['FIBERID', 'RA', 'DEC', 'XFOCAL', 'YFOCAL']
                            ][fiberidx])
                    if specidx == 2:
                        plugmap = astropy.table.vstack(plugmaps)
                if specidx == 1 and band == 'b':
                    # Record observation metadata.
                    for key in obskeys:
                        try:
                            value = spFrame.header[key]
                        except KeyError:
                            value = -999  # invalid value for int/float types
                        obsvals[key].append(value)
                # Load the sky fiber data.
                fibers = spFrame.plug_map['FIBERID'][fiberidx].data
                assert np.all(fiberidx == spFrame.get_fiber_offsets([fibers]))
                if expidx == 0 and band == 'b':
                    expected_fibers.append(fibers)
                    if verbose:
                        print('Found {} sky fibers on spec{}: {}.'.format(
                            len(fibers), specidx,
                            ','.join([str(f) for f in fibers])))
                else:
                    if not np.all(fibers == expected_fibers[specidx - 1]):
                        print('Did not get expected fibers for {} exp {}'
                              .format(camera, expidx))
                data = spFrame.get_valid_data(
                    fibers, include_sky=True, include_wdisp=True,
                    use_ivar=True, pixel_quality_mask=valid_mask)
                if verbose:
                    print('Reading {} for exposure {} / {}...'
                          .format(camera, expidx + 1, nexp))
                assert data.shape == (len(fibers), 2 * ampsize[band][0])
                mask = spFrame.get_pixel_masks(fibers)
                masks[band].append(mask[:, use])
                # Identify pixels with valid data.
                valid = ~data['ivar'].mask
                bad_fibers = ~np.any(valid, axis=1)
                if verbose and np.any(bad_fibers):
                    print('  bad fibers: {}'.format(fibers[bad_fibers]))
                ivar = data['ivar'].data
                assert np.all(ivar[valid] > 0)
                ivars[band].append(ivar[:, use])
                # Load the superflat and trace vectors for sky fibers.
                superflat = spFrame.get_superflat(fibers)
                tracex = spFrame.hdulist[7].read()[fiberidx]
                # Load fiberflat and neff vectors from this camera's spFlat.
                name = exposures.get_exposure_name(expidx, camera, 'spFlat')
                path = mirror.get(finder.get_plate_path(plate, name))
                with fits.open(path) as spFlat:
                    fiberflat = spFlat[0].data[fiberidx]
                    neff = bossdata.plate.TraceSet(spFlat[3]).get_y()[fiberidx]
                if np.any(neff[valid] <= 0):
                    print(f'WARNING: neff <= 0 for {camera} {expidx} {tag}')
                # Lookup the per-amplifier readnoise values.
                readnoises = np.array([
                    spFrame.header['RDNOISE{}'.format(amp)]
                    for amp in range(4)], dtype=np.float32)
                # Determine which amplifier (0-3) each pixel along the trace is
                # read out by and scale to RMS readnoise per wavelength pixel.
                amp = amplifier[band](tracex)
                rdnoise = rdnoise_scale * readnoises[amp] * neff
                rdnoises[band].append(rdnoise[:, use].astype(np.float32))
                # Combine the superflat and fiberflat.
                flat = superflat * fiberflat
                assert np.all(flat[valid] > 0)
                flats[band].append(flat[:, use])
                # Save wavelength solutions in angstroms.
                wlen = data['wavelength'].data
                wlens[band].append(wlen[:, use])
                # Save wavelength dispersions in angstroms.
                wdisp = data['wdisp'].data
                assert np.all(wdisp[valid] > 0)
                wdisp = wlen * np.expm1(wdisp_const * wdisp)
                wdisps[band].append(wdisp[:, use])
                # Save the combined flat-fielded sky models + residuals,
                # which might be negative due to readnoise.
                flux = data['flux'].data + data['sky'].data
                fluxes[band].append(flux[:, use])
    # Build observation metadata table.
    obslist = astropy.table.Table()
    for key in obskeys:
        obslist[key] = obsvals[key]
    # Build the output HDU list.
    hdus = fits.HDUList()
    cards = dict(PLATE=plate, MJD=mjd, NFIBERS=len(plugmap), NEXP=nexp)
    hdus.append(fits.PrimaryHDU(header=fits.Header(cards)))
    hdus.append(fits.table_to_hdu(obslist))
    hdus[-1].name = 'OBSLIST'
    hdus.append(fits.table_to_hdu(plugmap))
    hdus[-1].name = 'PLUGMAP'
    for band in 'b', 'r':
        Band = band.upper()
        # Combine arrays for each band and save an an image HDU.
        hdus.append(fits.ImageHDU(np.vstack(wlens[band]),
                                  name='{}WLEN'.format(Band)))
        hdus.append(fits.ImageHDU(np.vstack(wdisps[band]),
                                  name='{}WDISP'.format(Band)))
        hdus.append(fits.ImageHDU(np.vstack(rdnoises[band]),
                                  name='{}RDNOISE'.format(Band)))
        hdus.append(fits.ImageHDU(np.vstack(flats[band]),
                                  name='{}FLAT'.format(Band)))
        hdus.append(fits.ImageHDU(np.vstack(fluxes[band]),
                                  name='{}FLUX'.format(Band)))
        hdus.append(fits.ImageHDU(np.vstack(ivars[band]),
                                  name='{}IVAR'.format(Band)))
        hdus.append(fits.ImageHDU(np.vstack(masks[band]),
                                  name='{}MASK'.format(Band)))
    name = os.path.join(output_path, 'sky-{}-{}.fits'.format(plate, mjd))
    hdus.writeto(name, overwrite=True)
    print('Completed {}'.format(tag))
    return obslist


def smooth_variance(hdus, min_valid_frac=0.9, n_pca=20,
                    denoise_steps=3, verbose=False):
    """Smooth variance estimates using weighted PCA.
    """
    hdr = hdus[0].header
    plate = hdr['PLATE']
    mjd = hdr['MJD']
    nfibers = hdr['NFIBERS']
    nexp = hdr['NEXP']
    tag = f'PLATE {plate:05d} MJD {mjd:05d}'
    if verbose:
        print(f'Smoothing {tag} with {nexp} exposures for {nfibers} fibers.')

    pixel_mask = (1 << 32) - 1 - bossdata.bits.bitmask_from_text(
        bossdata.bits.SPPIXMASK, 'BRIGHTSKY|BADSKYCHI|REDMONSTER')

    for j, band in enumerate('BR'):

        # Look up arrays for this band.
        wlen = hdus[band + 'WLEN'].data
        rdnoise = hdus[band + 'RDNOISE'].data
        flat = hdus[band + 'FLAT'].data
        nelec = hdus[band + 'FLUX'].data * flat
        nspectra, npixels = wlen.shape

        # Only use pixels with IVAR>0 and only allowed mask bits set.
        ivar = hdus[band + 'IVAR'].data
        ivar[~np.isfinite(ivar)] = 0
        valid = (ivar > 0) & (
            np.bitwise_and(hdus[band + 'MASK'].data, pixel_mask) == 0)

        # Count the number of valid pixels in each spectrum.
        nvalid = np.count_nonzero(valid, axis=1)

        # Drop any spectra without enough valid pixels.
        good_spectra = nvalid >= 0.9 * npixels
        good_spectra_idx = np.where(good_spectra)[0]
        wlen = wlen[good_spectra]
        rdnoise = rdnoise[good_spectra]
        nelec = nelec[good_spectra]
        valid = valid[good_spectra]
        if verbose:
            print('Dropped {} / {} spectra with no good pixels.'
                  .format(np.count_nonzero(~good_spectra), nspectra))

        # Iteratively drop pixels until every element of the
        # npixels x npixels matrix VALID.T VALID is True.
        wgt = 1. * valid
        matrix = np.dot(wgt.T, wgt)
        good_pixels = np.ones(npixels, bool)
        ngood = npixels
        while True:
            good_pixels_idx = np.where(good_pixels)[0]
            submatrix = matrix[np.ix_(good_pixels_idx, good_pixels_idx)]
            if np.all(submatrix > 0):
                break
            # Find the pixel responsible for the row with the most zeros.
            nzeros = np.count_nonzero(submatrix == 0, axis=0)
            most = np.argmax(nzeros)
            worst = good_pixels_idx[most]
            good_pixels[worst] = False
            ngood -= 1
            if ngood < min_valid_frac * npixels:
                print(f'Not enough valid {band} pixels for wPCA in {tag}')
                return False

        if verbose:
            print('Dropped {} / {} pixels for non-zero weight matrix.'
                  .format(np.count_nonzero(~good_pixels), npixels))
        wlen = wlen[:, good_pixels]
        rdnoise = rdnoise[:, good_pixels]
        nelec = nelec[:, good_pixels]
        valid = valid[:, good_pixels]

        # Clip nelec to zero.
        nelec = np.clip(nelec, a_min=0, a_max=None)

        # Iteratively solve for a de-noised shot-noise variance.
        evar = nelec.copy()
        niter = 1
        median_nelec = np.maximum(1, np.median(nelec, axis=0))
        while niter <= denoise_steps:
            weights = (
                np.clip(evar, a_min=0, a_max=None) + rdnoise ** 2) ** -0.5
            weights[~valid] = 0.
            new_evar = wpca.WPCA(n_components=n_pca).fit_reconstruct(
                nelec, weights=weights)
            # Clip estimated shot noise to zero.
            new_evar = np.clip(new_evar, a_min=0, a_max=None)
            # Measure mean squared error between original and new estimates.
            MSE = np.mean(((new_evar - evar) / median_nelec) ** 2)
            if verbose:
                print('niter {} MSE {}'.format(niter, MSE))
            evar = new_evar
            niter += 1

        # Build a de-noised ivar in the original indexing.
        eivar = np.zeros_like(evar)
        eivar[valid] = (evar[valid] + rdnoise[valid] ** 2) ** -1
        eivar_full = np.zeros_like(hdus[band + 'IVAR'].data)
        eivar_full[np.ix_(good_spectra_idx, good_pixels_idx)] = eivar

        # Scale by the flat so that xEIVAR is directly comparable to xIVAR.
        eivar_full *= flat ** 2

        # Save or replace this array in the HDU list.
        hdu_name = band + 'EIVAR'
        if hdu_name in hdus:
            hdus[hdu_name].data = eivar_full
        else:
            hdus.append(fits.ImageHDU(eivar_full, name=hdu_name))

    return True


def downsample(hdus, flatfielded=True, verbose=True):
    """Downsample each spectrum to 100A bands in flux density units.

    Must run sky_smooth() before calling this function.
    """
    hdr = hdus[0].header
    plate = hdr['PLATE']
    mjd = hdr['MJD']
    nfibers = hdr['NFIBERS']
    nexp = hdr['NEXP']
    tag = f'PLATE {plate:05d} MJD {mjd:05d}'
    if verbose:
        print(f'Downsampling {tag}.')

    # Downsample to 100A bands from 3700-10400A.
    edges = 100. * np.arange(37, 104)
    nbands = len(edges) - 1
    dsflux = np.zeros((nexp * nfibers, nbands))
    dsivar = np.zeros((nexp * nfibers, nbands))

    # Work in units of (optionally flat-fielded) electrons per Angstrom per
    # second, converted to 1e-13 ergs assuming each electron corresponds
    # to one photon. The multiplier 1e-13 is chosen to give values ~1.
    hc = (astropy.constants.h * astropy.constants.c).to(
        1e-13 * u.erg * u.Angstrom).value
    exposures = astropy.table.Table(hdus['OBSLIST'].data)
    exptime = np.repeat(exposures['EXPTIME'].data, nfibers).reshape(-1, 1)
    for band in 'BR':
        wlen = hdus[band + 'WLEN'].data
        energy_per_photon = hc / wlen
        scale = energy_per_photon / (exptime * np.gradient(wlen, axis=1))
        eflux = hdus[band + 'FLUX'].data * scale
        eivar = hdus[band + 'EIVAR'].data / scale ** 2
        if not flatfielded:
            flat = hdus[band + 'FLAT'].data
            eflux *= flat
            nonzero = flat > 0
            eivar[nonzero] /= flat[nonzero] ** 2
        eflux_wgtd = eivar * eflux
        # Downsample each individual spectrum.
        for i in range(nexp * nfibers):
            idx = np.searchsorted(wlen[i], edges)
            for j in range(nbands):
                if idx[j] == idx[j + 1]:
                    continue
                sl = slice(idx[j], idx[j + 1])
                wsum = eivar[i, sl].sum()
                if wsum > 0:
                    dsflux[i, j] += eflux_wgtd[i, sl].sum()
                    dsivar[i, j] += wsum
    nonzero = (dsivar > 0)
    dsflux[nonzero] /= dsivar[nonzero]

    # Save or replace the downsampled arrays in the HDU list.
    hdu_name = 'DSFLUX' if flatfielded else 'DSFLUXE'
    if hdu_name in hdus:
        hdus[hdu_name].data = dsflux
    else:
        hdus.append(fits.ImageHDU(dsflux, name=hdu_name))
    hdu_name = 'DSIVAR' if flatfielded else 'DSIVARE'
    if hdu_name in hdus:
        hdus[hdu_name].data = dsivar
    else:
        hdus.append(fits.ImageHDU(dsivar, name=hdu_name))

    return True
