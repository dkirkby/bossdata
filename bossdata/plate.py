# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Access BOSS plate data products.
"""

from __future__ import division, print_function

class Plan(object):
    """The plan file for configuring the BOSS pipeline to combine exposures of a single plate.

    The datamodel for plan files is at
    http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/spPlan.html.
    Combined plan files are small text files that list the per-spectrograph (b1,b2,r1,r2)
    spFrame files used for a single coadd.

    Args:
        path(str): The local path to a plan file.
    """
    def __init__(self, path):
        self.plate = None
        self.exposures = {}
        with open(path,'r') as f:
            for line in f:
                if not line.startswith('SPEXP'):
                    continue
                tokens = line.split()
                # token[1] is the plate number and must be identical for all lines.
                if not self.plate:
                    self.plate = tokens[1]
                elif self.plate != tokens[1]:
                    raise RuntimeError('Internal format error: unexpected plate {0} in {1}.'.format(
                        tokens[1],path))
                # Build an exposure record from the remaining tokens.
                exposure = dict(MJD=tokens[2], EXPTIME=float(tokens[5]),
                    B1=tokens[7], B2=tokens[8], R1=tokens[9], R2=tokens[10])
                # Record this exposure under the appropriate category.
                flavor = tokens[4]
                if flavor in self.exposures:
                    self.exposures[flavor].append(exposure)
                else:
                    self.exposures[flavor] = [exposure]
        print(self.exposures)
