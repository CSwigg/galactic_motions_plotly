import plotly.graph_objects as go
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
from astropy import coordinates
from astropy.io import fits
from skimage import transform
import astropy.units as u
import numpy as np
import pandas as pd
from skimage import transform
import math

df_zucker = pd.read_csv('/Users/cam/Downloads/clouds_galactic.csv')

#we are going to downscale to make the file smaller so it's more easily saveable
data, header = fits.getdata("/Users/cam/Downloads/vergely_3D_Dust.fits",
                            header=True)


def extra_traces():
    traces = []
    R_1 = -6.1 * np.ones(1000)
    R_2 = -8.1 * np.ones(1000)
    R_3 = -10.1 * np.ones(1000)
    phi = np.linspace(-180, 180, 1000)
    z = np.zeros(1000)
    gc_line_1 = Galactocentric(rho=R_1 * u.kpc,
                               phi=phi * u.deg,
                               z=z * u.pc,
                               representation_type='cylindrical')
    gc_line_2 = Galactocentric(rho=R_2 * u.kpc,
                               phi=phi * u.deg,
                               z=z * u.pc,
                               representation_type='cylindrical')
    gc_line_3 = Galactocentric(rho=R_3 * u.kpc,
                               phi=phi * u.deg,
                               z=z * u.pc,
                               representation_type='cylindrical')
    gc_line_1 = gc_line_1.transform_to(Galactic)
    gc_line_2 = gc_line_2.transform_to(Galactic)
    gc_line_3 = gc_line_3.transform_to(Galactic)
    gc_line_1.representation_type = 'cartesian'
    gc_line_2.representation_type = 'cartesian'
    gc_line_3.representation_type = 'cartesian'
    gc_lines = coordinates.concatenate([gc_line_1, gc_line_2, gc_line_3])
    gc_line_ht = np.concatenate([['R = 6 kpc'] * len(gc_line_1),
                                 ['R = 8 kpc'] * len(gc_line_2),
                                 ['R = 10 kpc'] * len(gc_line_3)])

    scatter_zucker = go.Scatter3d(x=df_zucker.X,
                                  y=df_zucker.Y,
                                  z=df_zucker.Z,
                                  mode='markers',
                                  marker=dict(size=4,
                                              color='white',
                                              symbol='circle',
                                              opacity=.2),
                                  visible='legendonly',
                                  hovertext=df_zucker['name'].values,
                                  name='Zucker Clouds')
    scatter_gc_line_1 = go.Scatter3d(x=gc_lines.u.value * 1000,
                                     y=gc_lines.v.value * 1000,
                                     z=gc_lines.w.value * 1000,
                                     mode='markers',
                                     marker=dict(size=1.,
                                                 color='white',
                                                 symbol='circle',
                                                 opacity=1.),
                                     name='Constant GC Lines',
                                     hovertext=gc_line_ht)


    dust = vergely_dust()


    traces.append(scatter_zucker)
    traces.append(scatter_gc_line_1)
    traces.append(dust)


    return traces


def vergely_dust():

    data, header = fits.getdata("/Users/cam/Downloads/vergely_3D_Dust.fits",
                                header=True)
    downscale = 6.
    downcube = transform.pyramid_reduce(data,
                                        downscale=downscale,
                                        multichannel=False)
    x_origin, y_origin, z_origin = header['CRVAL1'], header['CRVAL2'], header[
        'CRVAL3']

    x_y_step = math.ceil(601 / downscale)
    z_step = math.ceil(81 / downscale)

    X, Y, Z = np.mgrid[x_origin:3000:x_y_step * 1j,
                       y_origin:3000:x_y_step * 1j, z_origin:400:z_step * 1j]
    X = X.flatten().astype(int)
    Y = Y.flatten().astype(int)
    Z = Z.flatten().astype(int)
    downcube = downcube.T.flatten()

    cut_value = 2500
    x_y_condition = np.where((X >= -cut_value) & (X <= cut_value) & (Y >= -cut_value)
                             & (Y <= cut_value))
    X_cut = X[x_y_condition]
    Y_cut = Y[x_y_condition]
    Z_cut = Z[x_y_condition]
    downcube_cut = downcube[x_y_condition]

    vol_plot = go.Volume(
        x=X_cut,
        y=Y_cut,
        z=Z_cut,
        value=downcube_cut,
        flatshading=True,
        opacity=.7,
        #isomin=np.percentile(downcube, 25),
        #isomax=np.percentile(downcube, 75),
        #isomin=150e-6,
        #isomax=1e-4,
        isomin=500e-6,
        showscale=False,
        colorscale='gray',
        opacityscale='max',
        #opacityscale = [[0, 1], [0.5, 0.2], [1, 1]],
        reversescale=True,
        surface=dict(show=True, count=8),
        spaceframe=dict(show=True),  #,
        contour=dict(show=False, width=8),
        hoverinfo='skip',
        visible = 'legendonly',
        name = 'Vergely+2022 Dust',
        showlegend = True
        )


    return vol_plot