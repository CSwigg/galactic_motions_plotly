import plotly.graph_objects as go
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
from astropy import coordinates
from astropy.io import fits
from skimage import transform
import astropy.units as u
import numpy as np
import pandas as pd

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


    # downscale = 3
    # downcube = transform.pyramid_reduce(data,
    #                                     downscale=downscale,
    #                                     multichannel=False)
    # x_origin, y_origin, z_origin = header['CRVAL1'], header['CRVAL2'], header[
    #     'CRVAL3']
    # X, Y, Z = np.mgrid[x_origin:x_origin + header['NAXIS1']:10j,
    #                 y_origin:y_origin + header['NAXIS2']:10j,
    #                 z_origin:z_origin + header['NAXIS3']:10j]
    # vol_plot = go.Volume(
    #     x=X.flatten().astype(int),
    #     y=Y.flatten().astype(int),
    #     z=Z.flatten().astype(int),
    #     value=downcube.T.flatten(),
    #     flatshading=True,
    #     opacity=0.5,
    #     #isomin=0.,
    #     #isomax=0.01,
    #     showscale=False,
    #     colorscale=[[0, 'white'], [1., 'gray']],
    #     opacityscale='max',
    #     reversescale=False,
    #     surface=dict(show=True, count=5),
    #     spaceframe=dict(show=True),  #,
    #     contour=dict(show=False, width=4),
    #     name = 'Vergely Dust'
    #     )
    #traces.append(vol_plot)

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
    traces.append(scatter_zucker)
    traces.append(scatter_gc_line_1)
    return traces
