# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Support for querying the metadata associated with BOSS observations.
"""

from __future__ import division, print_function


def focal_plane(xfocal, yfocal, data, mask=None, subsets=dict(),
                background=None, numbered=None, percentile_cut=0.0,
                mesh_refinement=0, plot_label=None, data_label=None,
                show_background_mesh=False, number_color='red',
                default_options=dict(marker='o', lw=0.5, s=60), rmax=350.0):
    """
    Scatter plot of data values...
    """
    fig = plt.figure(figsize=(11,8))
    fig.set_tight_layout(True)
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
        levels = np.linspace(vmin, vmax, 100)
        plt.tricontourf(refined, zfield, levels=levels, vmin=vmin, vmax=vmax)
        if show_background_mesh:
            plt.triplot(triangulation, color='gray')

    # Draw each fiber using the specified marker.
    remaining = valid
    for subset_name, subset_config in subsets.iteritems():
        try:
            subset = remaining & subset_config['fibers']
        except Exception as e:
            raise ValueError('Invalid subset configuration for {0}.'.format(subset_name))
        print subset_name, np.count_nonzero(subset)
        options = dict(default_options)
        options['label'] = subset_name
        if 'options' in subset_config:
            options.update(subset_config.get('options'))
        plt.scatter(xfocal[subset], yfocal[subset],
                    vmin=vmin, vmax=vmax, c=data[subset], **options)
        remaining = remaining & ~subset
    if subsets:
        plt.legend(loc='lower left')
    # Plot any remaining fibers using the default options.
    plt.scatter(xfocal[remaining], yfocal[remaining],
                vmin=vmin, vmax=vmax, c=data[remaining], **default_options)
    plt.gca().set_aspect(1.0)
    colors = plt.colorbar()
    if data_label is not None:
        colors.set_label(data_label)

    # Add fiber-number annotations.
    if numbered is not None:
        for number, annotated in enumerate(numbered):
            if not annotated:
                continue
            xy = (xfocal[number], yfocal[number] + 10)
            label = '{:d}'.format(number + 1)
            plt.annotate(label, xy = xy, xytext=xy, xycoords='data', textcoords='data',
                         horizontalalignment='center', verticalalignment='bottom',
                         fontsize=16, color=number_color)

    # Draw the plot label in the top-left corner.
    if plot_label is not None:
        xy = 0.01, 0.99
        plt.annotate(plot_label, xy=xy, xytext=xy, xycoords='axes fraction',
                     textcoords='axes fraction',
                     horizontalalignment='left', verticalalignment='top',
                     fontsize=24, color='black')
