import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
import astropy.units as u
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.potential import plotRotcurve
from galpy.potential import NFWPotential, HernquistPotential, MiyamotoNagaiPotential
from astropy.io import fits
import plotly.graph_objects as go
from astropy import coordinates

from conversion_functions import convert_icrs_galactic
import plotly_layout


default_marker_properties = dict(color = 'white', size = 4)


class MovieStar:

    def __init__(self, name, time_to_integrate, df = None, key_dict = None, marker_properties = default_marker_properties):
        self.name = name
        self.key_dict = key_dict
        self.time_to_integrate = time_to_integrate
        self.marker_properties = marker_properties

        if df is not None:
            self.df = df.rename(columns = self.key_dict) # rename cols based on user-specification
            self.df = self.df[list(self.key_dict.values())] # limit dataframe to only necessary cols
        self.df_integrated = self.generate_orbit_dataframe()

    def __str__(self):
        return "Name : {}, Available Keys : {}".format(self.name, self.df.keys())

    def generate_orbit_dataframe(self):   

        if self.name == 'sun':
            orbit = Orbit([0 * u.deg, 0 * u.deg, 0 * u.pc, 0 * u.mas / u.yr, 0 * u.mas / u.yr, 0 * u.km / u.s],
                               lb=False, uvw=False, radec=True,
                               ro=Galactocentric.galcen_distance,
                               zo=Galactocentric.z_sun,
                               vo=Galactocentric.galcen_v_sun.d_y - GalacticLSR.v_bary.d_y,
                               solarmotion=u.Quantity([-GalacticLSR.v_bary.d_x, GalacticLSR.v_bary.d_y, GalacticLSR.v_bary.d_z]))
            
        else:
            icrs = SkyCoord(ra = self.df.ra.values*u.deg,
                            dec = self.df.dec.values*u.deg,
                            distance = (1000/self.df.parallax.values)*u.pc,
                            pm_ra_cosdec = self.df.pmra.values*(u.mas/u.yr),
                            pm_dec = self.df.pmdec.values*(u.mas/u.yr),
                            radial_velocity = self.df.radial_velocity.values*(u.km/u.s)
                        )
            orbit = Orbit(vxvv = icrs)

        orbit.integrate(self.time_to_integrate*u.Myr, pot = MWPotential2014)
        ra_int = orbit.ra(self.time_to_integrate*u.Myr).flatten()
        dec_int = orbit.dec(self.time_to_integrate*u.Myr).flatten()
        distance_int = orbit.dist(self.time_to_integrate*u.Myr).flatten()*1000

        if self.name == 'sun':
            t_list = []
            for t in self.time_to_integrate:
                t_list.append(t)
            df_integrated = pd.DataFrame({'t':t_list,
                    'ra':ra_int, 
                    'dec':dec_int, 
                    'distance':distance_int,
                    'symbol':['circle']*len(t_list)})

        else:
            group_name_list = []
            age_list = []
            t_list = []
            for name, age in zip(self.df.group_name.values, self.df.age.values):
                for t in self.time_to_integrate:
                    group_name_list.append(name)
                    age_list.append(age)
                    t_list.append(t)

            df_integrated = pd.DataFrame({'t':t_list,
                                'ra':ra_int, 
                                'dec':dec_int, 
                                'distance':distance_int,
                                'age':age_list,
                                'group_name':group_name_list,
                                'symbol':['circle']*len(t_list)
                                })
        df_integrated = convert_icrs_galactic(df_integrated) # adds cartesian coords to dataframe

        return df_integrated


class Movie:

    def __init__(self, movie_stars : dict, time : np.ndarray, movie_save_path : str):
        self.movie_stars = movie_stars # list of MovieStar instances
        self.movie_save_path = movie_save_path
        self.time = time

    def make_movie(self):
        '''
        Generates movie of MovieStar intances, centered on the Sun's motion
        '''
        self.create_sun_orbit()
        self.movie_stars.append(self.sun)

        figure = {
            'data': [],
            'layout': {},
            'frames': [],
            'config': {'scrollzoom': True}
        }
        sliders_dict = {
            'yanchor': 'top',
            'xanchor': 'left',
            'transition': {'duration': 2, 'easing': 'bounce-in'},
            'pad': {'b': 0, 't': 10, 'l':190, 'r':250},
            'len': 0.9,
            'currentvalue' : {"prefix": "Time (Myr): "},
            'x': 0.1,
            'y': 0,
            'active':0,
            'steps': []
        }

        #figure['layout'] = plotly_layout.layout

        for t in self.time:
            frame = {'data': [], 'name': str(t)}

            for movie_star in self.movie_stars:
                df = movie_star.df_integrated # integrated dataframe of MovieStar
                #assert df.t == self.time # make sure MovieStar instance has been integrated to match time-steps of movie
                df_t = df.loc[df.t == t]

                if movie_star.name != 'sun':
                    df_t.loc[df_t.age <= np.abs(t), 'symbol'] = 'circle-open'
                    hovertext = df_t['group_name'].values
                else:
                    hovertext = ['Sun']
                movie_star.marker_properties['symbol'] = df_t['symbol'].values

                trace_3d_frame = go.Scatter3d(x = df_t.x.values,
                                              y = df_t.y.values, 
                                              z = df_t.z.values,
                                              mode = 'markers',
                                              marker = movie_star.marker_properties,
                                              line = dict(width = 0),
                                              hovertext = hovertext,
                                              name = movie_star.name
                                             )
                frame['data'].append(trace_3d_frame)

            x_min = self.sun.df_integrated.loc[self.sun.df_integrated.t == t].x.iloc[0] - 1000
            x_max = self.sun.df_integrated.loc[self.sun.df_integrated.t == t].x.iloc[0] + 1000
            y_min = self.sun.df_integrated.loc[self.sun.df_integrated.t == t].y.iloc[0] - 1000
            y_max = self.sun.df_integrated.loc[self.sun.df_integrated.t == t].y.iloc[0] + 1000

            frame['layout'] = go.Layout(scene = dict(
                            aspectmode = 'manual',
                            aspectratio = dict(x=1., y=1., z=.4),
                            xaxis = dict(range = [x_min, x_max], 
                                            showgrid = False, 
                                            zeroline = False),
                            yaxis = dict(range = [y_min, y_max], 
                                            showgrid = False, 
                                            zeroline = False),
                            zaxis = dict(range = [-400, 400], 
                                            showgrid = False, 
                                            zeroline = False)
                            ))

            figure['frames'].append(go.Frame(frame))
            slider_step = {'args': [
                [t],
                {'frame': {'duration': 5, 'redraw': True},
                'mode': 'immediate',
            'transition': {'duration': 100}}
            ],
            'label': np.round(t,2),
            'name':'Time',
            'method': 'animate'
            }
            sliders_dict['steps'].append(slider_step)
        figure['data'] = figure['frames'][0]['data']
        figure['layout'] = figure['frames'][0]['layout']
        figure['layout']['sliders'] = [sliders_dict]
        figure['layout']['template'] = 'plotly_dark'
        fig = go.Figure(figure) 
        fig = plotly_layout.extra_traces(fig)
        fig.write_html(self.movie_save_path, auto_open = False)
        print('Movie Created')



    def create_sun_orbit(self):
        sun_marker_props = {'size' : 5, 'color' : 'yellow', 'opacity' : 1., 'symbol' : 'circle'}
        self.sun = MovieStar(name = 'sun', time_to_integrate = self.time, marker_properties=sun_marker_props)