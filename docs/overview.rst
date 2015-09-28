===================================
Overview of SDSS Spectroscopic Data
===================================

This package is primarily intended for working with data from the `SDSS-III BOSS survey <https://www.sdss3.org/surveys/boss.php>`_, but can also be used to access older data from SDSS-I/II and newer data from the `SEQUELS ancillary program <http://www.sdss.org/dr12/algorithms/ancillary/boss/sequels/>`_ and the `SDSS-IV eBOSS survey <http://www.sdss.org/surveys/eboss/>`_ (see :doc:`/envvar` for details).

BOSS data consists of `spectroscopic observations <http://www.sdss.org/dr12/spectro/spectro_basics/>`_ of astrophysical `targets <http://www.sdss.org/dr12/algorithms/boss_target_selection/>`_. An observation is identified by a triplet of numbers (PLATE,MJD,FIBER). Most BOSS targets only have a single observation. Each observation consists of several 15-minute exposures using red and blue cameras with overlapping wavelength coverage that are combined to give a single co-added spectrum.

The table below summarizes the different files produced by the `spectroscopic pipeline <http://www.sdss.org/dr12/spectro/pipeline/>`_ containing the individual and combined exposures contributing to each observation. Files contain from 1 to 1000 spectra, with some duplication between files.  Each file provides wavelength, flux, inverse variance, `mask bits <https://www.sdss3.org/dr9/algorithms/bitmask_sppixmask.php>`_ and subtracted sky for each of its spectra.

====== ====== ===== ==== ====== ====== ============================================ =================================
Type   Size   #Tgts #Exp Coadd? Calib? Datamodel                                    Bossdata Class
====== ====== ===== ==== ====== ====== ============================================ =================================
lite   0.2Mb      1    0      Y      Y :datamodel:`lite <spectra/lite/PLATE4/spec>` :class:`bossdata.spec.SpecFile`
spec   1.7Mb      1  ALL      Y      Y :datamodel:`spec <spectra/PLATE4/spec>`      :class:`bossdata.spec.SpecFile`
plate  110Mb   1000    0      Y      Y :datamodel:`plate <PLATE4/spPlate>`          :class:`bossdata.plate.PlateFile`
cframe 75Mb     500    1      N      Y :datamodel:`cframe <PLATE4/spCFrame>`        :class:`bossdata.plate.FrameFile`
frame  30Mb     500    1      N      N :datamodel:`frame <PLATE4/spFrame>`          :class:`bossdata.plate.FrameFile`
====== ====== ===== ==== ====== ====== ============================================ =================================

The following examples show how the same combined spectrum can be :ref:`plotted <bossplot>` from lite files and plate files::

    bossplot --plate 6641 --mjd 56383 --fiber 30
    bossplot --plate 6641 --mjd 56383 --fiber 30 --platefile

Individual exposures can also be plotted using either spec files, cframe files or frame files::

    bossplot --plate 6641 --mjd 56383 --fiber 30 --exposure 0
    bossplot --plate 6641 --mjd 56383 --fiber 30 --exposure 2 --cframe
    bossplot --plate 6641 --mjd 56383 --fiber 30 --exposure 2 --frame

Note that the indexing of exposures is different for spec files, which only index exposures used in the final coadd, and (c)frame files which index all available exposures. The indices used in the example all refer to exposure number 00158842, which can be verified by adding the ``--verbose`` option to these commands. The difference between the cframe and frame files is that the frame gives fluxes in units of flat-fielded detected electrons, before the step of calibrating fluxes using standard stars.
