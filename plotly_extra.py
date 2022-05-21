import plotly.graph_objects as go
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
from astropy import coordinates
import astropy.units as u
import numpy as np
import pandas as pd



def extra_traces(figure):
    df_zucker = pd.read_csv('/Users/cam/Downloads/clouds_galactic.csv')
    R_1 = -6.1*np.ones(1000)
    R_2 = -8.1*np.ones(1000)
    R_3 = -10.1*np.ones(1000)
    phi = np.linspace(-180,180,1000)
    z = np.zeros(1000)
    gc_line_1 = Galactocentric(rho=R_1*u.kpc, phi=phi*u.deg, z=z*u.pc, representation_type = 'cylindrical')
    gc_line_2 = Galactocentric(rho=R_2*u.kpc, phi=phi*u.deg, z=z*u.pc, representation_type = 'cylindrical')
    gc_line_3 = Galactocentric(rho=R_3*u.kpc, phi=phi*u.deg, z=z*u.pc, representation_type = 'cylindrical')
    gc_line_1 = gc_line_1.transform_to(Galactic)
    gc_line_2 = gc_line_2.transform_to(Galactic)
    gc_line_3 = gc_line_3.transform_to(Galactic)
    gc_line_1.representation_type = 'cartesian'
    gc_line_2.representation_type = 'cartesian'
    gc_line_3.representation_type = 'cartesian'
    gc_lines = coordinates.concatenate([gc_line_1,gc_line_2,gc_line_3])

    scatter_zucker = go.Scatter3d(x = df_zucker.X,
                            y = df_zucker.Y,
                            z = df_zucker.Z,
                            mode = 'markers',
                            marker = dict(size = 2,
                                            color = 'red',
                                            symbol = 'circle',
                                            opacity = 1.
                                        ),
                            visible = 'legendonly',
                            hovertext = df_zucker['name'].values,
                            name = 'Zucker Clouds'
                            )
    scatter_gc_line_1 = go.Scatter3d(x = gc_lines.u.value*1000,
                                    y = gc_lines.v.value*1000,
                                    z = gc_lines.w.value*1000,
                                    mode = 'markers',
                                    marker = dict(size = 1.5,
                                                color = 'white',
                                                symbol = 'circle',
                                                opacity = 1.
                                                ),
                                    name = 'Constant GC Lines'
                                )
    
    figure.add_trace(scatter_zucker)
    figure.add_trace(scatter_gc_line_1)
    return figure