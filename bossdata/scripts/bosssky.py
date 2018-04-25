"""Driver for command-line bosssky script
"""
from __future__ import division, print_function

import os.path

import numpy as np

import astropy.table
import astropy.io.fits as fits

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
    inverses variance xIVAR. xFLUX * xFLAT and xRDNOISE are in units of
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
    valid_slices = {'b': slice(767, 3299), 'r': slice(483, 3668) }
    # Initialize data access.
    finder = bossdata.path.Finder()
    mirror = bossdata.remote.Manager()
    # Loop over spectrographs.
    expected_fibers = []
    for specidx in 1, 2:
        # Load the list of science exposures used for this spectrograph's coadd.
        fiber = 500 * (specidx - 1) + 1
        spec_name = finder.get_spec_path(plate, mjd, fiber=fiber, lite=True)
        exposures = bossdata.spec.SpecFile(mirror.get(spec_name)).exposures
        for band in 'b', 'r':
            camera = '{}{}'.format(band, specidx)
            use = valid_slices[band]
            # Loop over science exposures for this camera.
            nexp = exposures.num_by_camera[camera]
            assert last_nexp is None or nexp == last_nexp, \
                'Cameras do not have matching exposures.'
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
                    plugmaps.append(spFrame.plug_map[
                        ['FIBERID','RA','DEC','XFOCAL','YFOCAL']][fiberidx])
                    if specidx == 2:
                        plugmap = astropy.table.vstack(plugmaps)
                if specidx == 1 and band == 'b':
                    # Record observation metadata.
                    for key in obskeys:
                        try:
                            value = spFrame.header[key]
                        except KeyError:
                            value = -999 # invalid value for int/float types
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
                    fibers, include_sky=True, include_wdisp=True, use_ivar=True,
                    pixel_quality_mask=valid_mask)
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
                assert np.all(neff[valid] > 0)
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
    if verbose:
        print('Wrote {}'.format(name))
    return obslist
