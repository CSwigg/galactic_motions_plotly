import plotly.graph_objects as go
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
from astropy import coordinates
import astropy.units as u
import numpy as np
import pandas as pd



layout = go.Layout(updatemenus = [dict(type="buttons", 
                                       buttons = [dict(label="Play", method="animate", args = [None]),
                                                  dict(label="Pause", method="animate", args = [[None], 
                                                  {'frame':{'duration':0, 'redraw':False},
                                                   'mode':'immediate',
                                                   'transition':{'duration':0}}])],
                                       direction='left',
                                       pad={'r':10, 't':40},
                                       showactive=False,
                                       x=0.58,
                                       xanchor='right',
                                       y=-0.15,
                                       yanchor='top')],
                    scene = dict(
                                aspectmode = 'manual',
                                aspectratio = dict(x=.8, y=1, z=.08),
                                xaxis = dict(range = [-2000, 6000], 
                                             showgrid = False, 
                                             zeroline = False),
                                yaxis = dict(range = [-8000, 2000], 
                                             showgrid = False, 
                                             zeroline = False),
                                zaxis = dict(range = [-400, 400], 
                                             showgrid = False, 
                                             zeroline = False)
                               )
                  )


def extra_traces(figure):
    df_zucker = pd.read_csv('/Users/cam/Downloads/clouds_galactic.csv')
    R_1 = -6.1*np.ones(500)
    R_2 = -8.1*np.ones(500)
    R_3 = -10.1*np.ones(500)
    phi = np.linspace(-90,90,500)
    z = np.zeros(500)
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
                                            opacity = .2
                                        ),
                            visible = 'legendonly',
                            hovertext = df_zucker['name'].values,
                            name = 'Zucker Clouds'
                            )
    scatter_gc_line_1 = go.Scatter3d(x = gc_lines.u.value*1000,
                                    y = gc_lines.v.value*1000,
                                    z = gc_lines.w.value*1000,
                                    mode = 'markers',
                                    marker = dict(size = 1.,
                                                color = 'white',
                                                symbol = 'circle',
                                                opacity = 1.
                                                ),
                                    name = 'Constant GC Lines'
                                )
    
    figure.add_trace(scatter_zucker)
    figure.add_trace(scatter_gc_line_1)
    return figure