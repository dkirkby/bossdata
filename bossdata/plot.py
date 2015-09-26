# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Support for plotting BOSS spectscopic data in different formats.

These functions use the optional matplotlib dependency so will raise an
ImportError if this is not installed.  Functions do not create figures or
call :func:`matplotlib.pyplot.show` before exiting, to provide the maximum
flexibility.  To display a single plot, you can use the following wrapper::

    plt.figure(figsize=(11,8.5))
    # ... call one of the plot functions ...
    plt.tight_layout()
    plt.show()

See the :doc:`/examples` for details.
"""

from __future__ import division, print_function

import numpy as np


def by_fiber(data, mask=None, subsets=dict(), percentile_cut=0.0,
             plot_label=None, data_label=None,
             default_options=dict(marker='o', lw=0.5, s=60)):
    """Plot per-fiber data values in fiber order.

    This is a useful plot to show any dependence of the data value on a
    fiber's position on the CCD and slithead.  Both spectrographs are
    superimposed on the same plot.  The points for each fiber are color-coded
    according to their associated data value using the same scheme as
    :func:`focal_plane`.

    Args:
        data(numpy.ndarray): A 1D array of data values to plot, where the
            array index matches the fiber number and all fibers are included.
        mask(numpy.ndarray): An optional 1D array of boolean values with True
            values used to mask out values in the data array.  Masked values
            will not be plotted and will not be used to calculate the plot
            data range.
        subsets(dict): A dictionary of fiber subsets that will be separately
            identified in the plot.  Each dictionary must define values for
            two keys: 'options' and 'fibers'.  The options are a dictionary
            of arguments passed to :func:`matplotlib.pyplot.scatter` and used
            to style the subset.  The fibers value is used to index the data
            array to pick out the subset's data values.
        percentile_cut(float): Data will be clipped to this percentile value on
            both sides of its distribution.  Use a value of zero (the default)
            for no clipping.
        plot_label(str): A label identifying this plot that will be displayed
            in the top-left corner.
        data_label(str): A label identifying the data values that will be
            used to label the y axis.
        default_options(dict): A dictionary of options passed to
            :func:`matplotlib.pyplot.scatter` that is used to draw data points.
            Options in a subset dictionary override any values here. Fibers not
            in any subset are drawn using these default options.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri

    if mask is not None:
        valid = ~mask
    else:
        valid = np.ones_like(data, dtype=bool)
    vmin, vmax = np.percentile(
        data[valid],(percentile_cut, 100. - percentile_cut))

    # Draw lines connecting the fibers of each spectrograph.
    num_fibers = len(data) // 2
    fiber_num = np.arange(1, num_fibers + 1)

    data1 = data[:num_fibers]
    data2 = data[num_fibers:]
    plt.plot(fiber_num, data1, '-', label='1-{}'.format(num_fibers))
    plt.plot(fiber_num, data2, ':', label='{}-{}'.format(
        num_fibers + 1, 2 * num_fibers))

    # Draw each fiber using the specified marker.
    remaining = valid
    for subset_name, subset_config in subsets.iteritems():
        try:
            subset = remaining & subset_config['fibers']
        except Exception as e:
            raise ValueError(
                'Invalid subset configuration for {0}.'.format(subset_name))
        mask1 = subset[:num_fibers]
        mask2 = subset[num_fibers:]
        options = dict(default_options)
        if 'options' in subset_config:
            options.update(subset_config.get('options'))
        plt.scatter(fiber_num[mask1], data1[mask1], vmin=vmin, vmax=vmax,
            c=data1[mask1], label=subset_name, **options)
        plt.scatter(fiber_num[mask2], data2[mask2], vmin=vmin, vmax=vmax,
            c=data2[mask2], **options)
        remaining = remaining & ~subset
    plt.legend(loc='lower left', scatterpoints=1, markerscale=1.0)
    # Plot any remaining fibers using the default options.
    mask1 = remaining[:num_fibers]
    mask2 = remaining[num_fibers:]
    plt.scatter(fiber_num[mask1], data1[mask1],
                vmin=vmin, vmax=vmax, c=data1[mask1], **default_options)
    plt.scatter(fiber_num[mask2], data2[mask2],
                vmin=vmin, vmax=vmax, c=data2[mask2], **default_options)

    plt.xlim(0, num_fibers + 1)
    if vmin < vmax:
        pad = 0.02 * (vmax - vmin)
        plt.ylim(vmin - pad, vmax + pad)

    plt.xlabel('Fiber Position')
    if data_label is not None:
        plt.ylabel(data_label)
    if plot_label is not None:
        xy = 0.01, 0.99
        plt.annotate(plot_label, xy=xy, xytext=xy, xycoords='axes fraction',
                     textcoords='axes fraction',
                     horizontalalignment='left', verticalalignment='top',
                     fontsize=20, color='black')


def focal_plane(xfocal, yfocal, data, mask=None, subsets=dict(),
                background=None, numbered=None, percentile_cut=0.0,
                mesh_refinement=0, plot_label=None, data_label=None,
                show_background_mesh=False, number_color='red',
                default_options=dict(marker='o', lw=0.5, s=60), rmax=350.0):
    """Plot per-fiber data values using focal-plane positions.

    This is a useful plot to show any dependence of the data value on a
    fiber's position in the focal plane.  The points for each fiber are
    color-coded according to their associated data value using the same scheme
    as :func:`by_fiber`.

    Args:
        xfocal(numpy.ndarray): A 1D array of x focal-plane positions, where
            the array index matches the fiber number and all fibers are
            included.
        yfocal(numpy.ndarray): A 1D array of y focal-plane positions, where
            the array index matches the fiber number and all fibers are
            included.
        data(numpy.ndarray): A 1D array of data values to plot, where the
            array index matches the fiber number and all fibers are included.
        mask(numpy.ndarray): An optional 1D array of boolean values with True
            values used to mask out values in the data array.  Masked values
            will not be plotted and will not be used to calculate the plot
            data range.
        subsets(dict): A dictionary of fiber subsets that will be separately
            identified in the plot.  Each dictionary must define values for
            two keys: 'options' and 'fibers'.  The options are a dictionary
            of arguments passed to :func:`matplotlib.pyplot.scatter` and used
            to style the subset.  The fibers value is used to index the data
            array to pick out the subset's data values.
        background(numpy.ndarray): An optional subset of fibers whose data
            values are used to fill the background using interpolation. The
            resulting background fill will only cover the convex hull of the
            subset, where interpolation is possible.
        numbered(numpy.ndarray): An optional subset of fibers that will be
            numbered in the generated plot.
        percentile_cut(float): Data will be clipped to this percentile value on
            both sides of its distribution.  Use a value of zero (the default)
            for no clipping.
        mesh_refinement(int): Smoothness of background fill interpolation to
            use. A value of zero (the default) corresponds to linear
            interpolation.
        plot_label(str): A label identifying this plot that will be displayed
            in the top-left corner.
        data_label(str): A label identifying the data values that will be
            used to label the y axis.
        show_background_mesh(bool): Draw the triangulation used for the
            background fill when this is True.
        number_color(str): Matplotlib color used to draw fiber numbers.
        default_options(dict): A dictionary of options passed to
            :func:`matplotlib.pyplot.scatter` that is used to draw data points.
            Options in a subset dictionary override any values here. Fibers not
            in any subset are drawn using these default options.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri

    plt.xlim(-1.01*rmax,+1.01*rmax)
    plt.ylim(-1.01*rmax,+1.01*rmax)
    plt.axis('off')

    # Plot the plate outline
    outline = plt.Circle((0,0), rmax, edgecolor='black', facecolor='none')
    plt.gcf().gca().add_artist(outline)

    # Find the limiting values to display.
    if mask is not None:
        valid = ~mask
    else:
        valid = np.ones_like(data, dtype=bool)
    vmin, vmax = np.percentile(
        data[valid],(percentile_cut, 100. - percentile_cut))

    if background is not None:
        # Draw a background triangulated mesh from the specified fibers.
        triangulation = matplotlib.tri.Triangulation(
            xfocal[valid & background], yfocal[valid & background])
        refiner = matplotlib.tri.UniformTriRefiner(triangulation)
        refined, zfield = refiner.refine_field(
            data[valid & background], subdiv=mesh_refinement)
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        levels = np.linspace(vmin, vmax, 100)
        plt.tricontourf(refined, zfield, levels=levels,
                        vmin=vmin, vmax=vmax, extend='both')
        if show_background_mesh:
            plt.triplot(triangulation, color='gray')

    # Draw each fiber using the specified marker.
    remaining = valid
    for subset_name, subset_config in subsets.iteritems():
        try:
            subset = remaining & subset_config['fibers']
        except Exception as e:
            raise ValueError(
                'Invalid subset configuration for {0}.'.format(subset_name))
        options = dict(default_options)
        options['label'] = subset_name
        if 'options' in subset_config:
            options.update(subset_config.get('options'))
        plt.scatter(xfocal[subset], yfocal[subset],
                    vmin=vmin, vmax=vmax, c=data[subset], **options)
        remaining = remaining & ~subset
    if subsets:
        plt.legend(loc='lower left', scatterpoints=1, markerscale=1.0)
    # Plot any remaining fibers using the default options.
    plt.scatter(xfocal[remaining], yfocal[remaining],
                vmin=vmin, vmax=vmax, c=data[remaining], **default_options)
    plt.gca().set_aspect(1.0)
    colors = plt.colorbar(pad=0.01)
    if data_label is not None:
        colors.set_label(data_label)

    # Add fiber-number annotations.
    if numbered is not None:
        for number, annotated in enumerate(numbered):
            if not annotated:
                continue
            xy = (xfocal[number], yfocal[number] + 10)
            label = '{:d}'.format(number + 1)
            plt.annotate(label, xy = xy, xytext=xy, xycoords='data',
                         textcoords='data', color=number_color,
                         horizontalalignment='center',
                         verticalalignment='bottom', fontsize=16)

    # Draw the plot label in the top-left corner.
    if plot_label is not None:
        xy = 0.01, 0.99
        plt.annotate(plot_label, xy=xy, xytext=xy, xycoords='axes fraction',
                     textcoords='axes fraction', fontsize=20, color='black',
                     horizontalalignment='left', verticalalignment='top')
