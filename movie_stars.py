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
import plotly_extra


default_marker_properties = dict(color='white', size=4)


class MovieStar:
    """ Class which encapsulates stellar information and integrates orbits with galpy
    """

    def __init__(self, name, time_to_integrate, df=None, orbit=None, key_dict=None, marker_properties=default_marker_properties):
        """ Initialize MovieStar 

        Parameters
        ----------
        name : str
            Name of MovieStar instance to be represented as a trace in Plotly plot
        time_to_integrate : np.ndarray
            Time-step array used for orbit integration
        df : pd.DataFrame, optional
            DataFrame of stellar information
        orbit : galpy.Orbit, optional
            galpy.Orbit instance, by default None
        key_dict : dict, optional
            Dictionary for mapping input DataFrame columns to required column keys, by default None
        marker_properties : dict, optional
            Dictionary specifying marker properties for  go.Scatter3d, by default default_marker_properties
        """
        self.name = name
        self.key_dict = key_dict
        self.time_to_integrate = time_to_integrate
        self.marker_properties = marker_properties

        if df is None:
            self.df = df
            self.orbit = orbit
        else:
            # rename cols based on user-specification
            self.df = df.rename(columns=self.key_dict)
            # limit dataframe to only necessary cols
            self.df = self.df[list(self.key_dict.values())]
            self.orbit = None
        self.df_integrated = self.generate_orbit_dataframe()

    def __str__(self):
        return "Name : {}, Available Keys : {}".format(self.name, self.df.keys())

    def generate_orbit_dataframe(self):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """
        if self.orbit is None:
            icrs = SkyCoord(ra=self.df.ra.values*u.deg,
                            dec=self.df.dec.values*u.deg,
                            distance=(1000/self.df.parallax.values)*u.pc,
                            pm_ra_cosdec=self.df.pmra.values*(u.mas/u.yr),
                            pm_dec=self.df.pmdec.values*(u.mas/u.yr),
                            radial_velocity=self.df.radial_velocity.values *
                            (u.km/u.s)
                            )
            self.orbit = Orbit(vxvv=icrs)

        self.orbit.integrate(self.time_to_integrate*u.Myr, pot=MWPotential2014)
        ra_int = self.orbit.ra(self.time_to_integrate*u.Myr).flatten()
        dec_int = self.orbit.dec(self.time_to_integrate*u.Myr).flatten()
        distance_int = self.orbit.dist(
            self.time_to_integrate*u.Myr).flatten()*1000

        if (self.df is None) or ('group_name' not in self.df.columns):
            t_list = []
            df_integrated = pd.DataFrame({'t': self.time_to_integrate,
                                          'ra': ra_int,
                                          'dec': dec_int,
                                          'distance': distance_int})

        else:
            # TODO: make this for-loop more flexible using dictionary

            group_name_list = []
            age_list = []
            t_list = []

            group_name_list = []
            age_list = []
            t_list = []
            for name, age in zip(self.df.group_name.values, self.df.age.values):
                for t in self.time_to_integrate:
                    group_name_list.append(name)
                    age_list.append(age)
                    t_list.append(t)

            df_integrated = pd.DataFrame({'t': t_list,
                                          'ra': ra_int,
                                          'dec': dec_int,
                                          'distance': distance_int,
                                          'age': age_list,
                                          'group_name': group_name_list
                                          })

            if self.marker_properties['size'] == 'size-based':
                size_list = []
                for n in self.df.n.values:
                    if n > 45:
                        n = 100  # create max of size scale
                    if n < 8:
                        n = 2  # create min of size scale
                    else:
                        n = 5
                    for t in self.time_to_integrate:
                        size_list.append(n)

                df_integrated['size'] = size_list

        # adds cartesian coords to dataframe
        df_integrated = convert_icrs_galactic(df_integrated)
        return df_integrated

    def create_age_based_symbols(self):
        """ 
        """
        self.df_integrated['symbol'] = ['circle']*len(self.df_integrated)

        if self.marker_properties['size'] != 'size-based':
            self.df_integrated['size'] = [
                self.marker_properties['size']]*len(self.df_integrated)

        # self.df_integrated.loc[self.df_integrated['age'] <=
        #                       self.df_integrated['t'].abs(), 'symbol'] = 'circle-open'
        self.df_integrated.loc[self.df_integrated['age'] <=
                               self.df_integrated['t'].abs(), 'size'] = 0.00001


class Movie:
    """ Class which encapsulates and controls plotly animation
    """

    def __init__(self, movie_stars: dict, time: np.ndarray, movie_save_path: str, xyz_ranges=[1500, 1500, 400], center_star='Sun', annotations=None):
        """ Initalize Movie class

        Parameters
        ----------
        movie_stars : dict
            list of MovieStar instances
        time : np.ndarray
            should be same as time_to_integrate in MovieStar class
        movie_save_path : str
            path to save animation as html file
        """
        self.movie_stars = movie_stars  # list of MovieStar instances
        self.movie_save_path = movie_save_path
        self.time = time
        self.center_star = center_star
        self.annotations = annotations
        self.xyz_ranges = xyz_ranges

    def make_movie(self):
        '''
        Generates movie of MovieStar intances, centered on the Sun's motion
        '''
        self.sun = self.create_sun_orbit()
        self.movie_stars.append(self.sun)

        self.figure = {
            'data': [],
            'layout': {},
            'frames': [],
            'config': {'scrollzoom': True}
        }
        self.sliders_dict = {
            'yanchor': 'top',
            'xanchor': 'left',
            'transition': {'duration': 2, 'easing': 'bounce-in'},
            'pad': {'b': 0, 't': 10, 'l': 190, 'r': 250},
            'len': 0.9,
            'currentvalue': {"prefix": "Time (Myr): "},
            'x': 0.1,
            'y': 0,
            'active': 0,
            'steps': []
        }

        # generates frames for each timestep and calls 'generate_frame_layout'
        self.generate_frames()

        self.figure['data'] = self.figure['frames'][0]['data']
        self.figure['layout'] = self.figure['frames'][0]['layout']

        self.figure['layout']['updatemenus'] = [dict(type="buttons",
                                                     buttons=[dict(label="Play", method="animate", args=[None]),
                                                              dict(label="Pause", method="animate", args=[[None],
                                                                                                          {'frame': {'duration': 50, 'redraw': False},
                                                                                                           'mode': 'immediate',
                                                                                                           'transition': {'duration': 10}}])],
                                                     direction='left',
                                                     pad={'r': 10, 't': 40},
                                                     showactive=False,
                                                     x=0.58,
                                                     xanchor='right',
                                                     y=-0.15,
                                                     yanchor='top')]

        self.figure['layout']['scene_camera'] = dict(up=dict(x=0, y=1, z=0),
                                                     center=dict(
                                                         x=0, y=0, z=0),
                                                     eye=dict(x=0., y=-.5, z=0.8))
        self.figure['layout']['dragmode'] = 'turntable'
        self.figure['layout']['legend'] = dict(x=0,
                                               y=1.,
                                               font=dict(size=10,
                                                         family='Verdana'
                                                         ),
                                               itemsizing='constant',
                                               bgcolor='rgba(0,0,0,0)'
                                               )

        self.figure['layout']['sliders'] = [self.sliders_dict]
        self.figure['layout']['template'] = 'plotly_dark'
        fig = go.Figure(self.figure)
        fig = plotly_extra.extra_traces(fig)
        fig.write_html(self.movie_save_path, auto_open=False)
        print('Movie Created!')

    def generate_frames(self):

        for movie_star in self.movie_stars:
            # if movie_star.marker_properties['size'] == 'size-based':
            #     marker_props['size'] = df_t['size'].values/20
            if movie_star.marker_properties['symbol'] == 'age-based':
                movie_star.create_age_based_symbols()

        for t in self.time:
            frame = {'data': [], 'name': str(t)}
            for movie_star in self.movie_stars:

                df = movie_star.df_integrated  # integrated dataframe of MovieStar
                marker_props = movie_star.marker_properties.copy()
                df_t = df.loc[df.t == t]

                if movie_star.marker_properties['symbol'] == 'age-based':
                    marker_props['symbol'] = df_t['symbol'].values
                    marker_props['size'] = df_t['size'].values
                    if movie_star.marker_properties['size'] == 'size-based':
                        marker_props['size'] = marker_props['size']

                if 'group_name' in df_t.columns:
                    hovertext = df_t['group_name'].values
                else:
                    hovertext = None

                trace_3d_frame = go.Scatter3d(x=df_t.x.values,
                                              y=df_t.y.values,
                                              z=df_t.z.values,
                                              mode='markers',
                                              marker=marker_props,
                                              line=dict(width=0.0001),
                                              hovertext=hovertext,
                                              name=movie_star.name
                                              )
                frame['data'].append(trace_3d_frame)
                frame['layout'] = self.generate_frame_layout(t)

            self.figure['frames'].append(go.Frame(frame))
            slider_step = {'args': [
                [t],
                {'frame': {'duration': 5, 'redraw': True},
                 'mode': 'immediate',
                 'transition': {'duration': 100}}
            ],
                'label': np.round(t, 1),
                'name': 'Time',
                'method': 'animate'
            }
            self.sliders_dict['steps'].append(slider_step)

    def generate_frame_layout(self, t):

        x_box_half_size = self.xyz_ranges[0]/2
        y_box_half_size = self.xyz_ranges[1]/2
        z_box_half_size = self.xyz_ranges[2]/2

        largest_box_size = max(self.xyz_ranges)

        for ms in self.movie_stars:  # select the stars to be centered on
            if ms.name == self.center_star:
                df_integrated = ms.df_integrated

        layout_annotations = []
        if self.annotations is not None:
            for annotation in self.annotations:  # enter annotations
                layout_annotations.append(annotation)

        x_min = df_integrated.loc[df_integrated.t ==
                                  t].x.iloc[0] - x_box_half_size
        x_max = df_integrated.loc[df_integrated.t ==
                                  t].x.iloc[0] + x_box_half_size
        y_min = df_integrated.loc[df_integrated.t ==
                                  t].y.iloc[0] - y_box_half_size
        y_max = df_integrated.loc[df_integrated.t ==
                                  t].y.iloc[0] + y_box_half_size

        layout = go.Layout(
            scene=dict(
                aspectmode='manual',
                aspectratio=dict(x=x_box_half_size/largest_box_size, y=y_box_half_size /
                                 largest_box_size, z=1.5*z_box_half_size/largest_box_size),
                xaxis=dict(range=[x_min, x_max],
                           showgrid=False,
                           zeroline=False,
                           nticks=3,
                           showline=True,
                           linecolor='white',
                           linewidth=3,
                           title=dict(text='X (pc)')
                           ),
                yaxis=dict(range=[y_min, y_max],
                           showgrid=False,
                           zeroline=False,
                           nticks=3,
                           showline=True,
                           linecolor='white',
                           linewidth=3,
                           title=dict(text='Y (pc)')
                           ),
                zaxis=dict(range=[-z_box_half_size, z_box_half_size],
                           showgrid=False,
                           zeroline=False,
                           nticks=3,
                           showline=True,
                           linecolor='white',
                           linewidth=3,
                           title=dict(text='Z (pc)')
                           ),
                annotations=layout_annotations
            )
        )

        return layout

    def create_sun_orbit(self):

        sun_marker_props = {'size': 5, 'color': 'yellow',
                            'opacity': 1., 'symbol': 'circle'}

        orbit_sun = Orbit([0 * u.deg, 0 * u.deg, 0 * u.pc, 0 * u.mas / u.yr, 0 * u.mas / u.yr, 0 * u.km / u.s],
                          lb=False, uvw=False, radec=True,
                          ro=Galactocentric.galcen_distance,
                          zo=Galactocentric.z_sun,
                          vo=Galactocentric.galcen_v_sun.d_y - GalacticLSR.v_bary.d_y,
                          solarmotion=u.Quantity([-GalacticLSR.v_bary.d_x, GalacticLSR.v_bary.d_y, GalacticLSR.v_bary.d_z]))

        sun = MovieStar(name='Sun', orbit=orbit_sun,
                             time_to_integrate=self.time, marker_properties=sun_marker_props)

        return sun

    # def create_anotations(self):
