import yaml
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
import astropy.units as u
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.potential import NFWPotential, HernquistPotential, MiyamotoNagaiPotential, NonInertialFrameForce, SpiralArmsPotential
from astropy.io import fits
import plotly.graph_objects as go
from astropy import coordinates

from conversion_functions import convert_icrs_galactic, coordFIX_to_coordROT
import plotly_extra

default_marker_properties = dict(color='white', size=4)


class MovieStar:
    """ Class which encapsulates stellar information and integrates orbits with galpy
    """

    def __init__(self,
                 name,
                 time_to_integrate,
                 df=None,
                 orbit=None,
                 key_dict=None,
                 marker_properties=default_marker_properties,
                 marker_behavior='open-circle',
                 visible=True):
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
        self.t_mirrored = np.concatenate(
            [np.flip(-1 * self.time_to_integrate[1:]), self.time_to_integrate])
        self.marker_properties = marker_properties
        self.marker_behavior = marker_behavior
        self.visible = visible  # not being used currently

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

    @classmethod
    def from_yaml(cls, yaml_file_name, time_to_integrate):
        yaml_file_directory = 'traces_yaml/'
        yaml_file_name = yaml_file_directory + yaml_file_name

        with open(yaml_file_name, 'r') as stream:
            yaml_data = yaml.safe_load(stream)

        df = pd.read_csv(yaml_data['file_name'])

        return cls(df=df,
                   name=yaml_data['trace_name'],
                   key_dict=yaml_data['key_dict'],
                   marker_properties=yaml_data['marker_properties'],
                   marker_behavior=yaml_data['marker_behavior'],
                   time_to_integrate=time_to_integrate)

    def __str__(self):
        return "Name : {}, Available Keys : {}".format(self.name,
                                                       self.df.keys())

    def generate_orbit_dataframe(self):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """

        if self.orbit is None:
            icrs = SkyCoord(ra=self.df.ra.values * u.deg,
                            dec=self.df.dec.values * u.deg,
                            distance=(1000 / self.df.parallax.values) * u.pc,
                            pm_ra_cosdec=self.df.pmra.values * (u.mas / u.yr),
                            pm_dec=self.df.pmdec.values * (u.mas / u.yr),
                            radial_velocity=self.df.radial_velocity.values *
                            (u.km / u.s))
            self.orbit = Orbit(vxvv=icrs)

        #mp= MiyamotoNagaiPotential(a=0.5,b=0.0375,normalize=1.)
        #mp = HernquistPotential(normalize=1.)
        mp = MWPotential2014
        potential = mp
        ra_int_f, dec_int_f, distance_int_f, pmra_int_f, pmdec_int_f, rv_int_f = integrate_orbit(
            time=self.time_to_integrate,
            orbit=self.orbit,
            pot=potential,
            flip=False)
        ra_int_b, dec_int_b, distance_int_b, pmra_int_b, pmdec_int_b, rv_int_b = integrate_orbit(
            time=-1 * self.time_to_integrate[1:],
            orbit=self.orbit,
            pot=potential,
            flip=True)

        axis = len(np.shape(ra_int_f)) - 1
        ra_int = np.concatenate([ra_int_b, ra_int_f], axis=axis).flatten()
        dec_int = np.concatenate([dec_int_b, dec_int_f], axis=axis).flatten()
        distance_int = np.concatenate([distance_int_b, distance_int_f],
                                      axis=axis).flatten()
        pmra_int = np.concatenate([pmra_int_b, pmra_int_f],
                                  axis=axis).flatten()
        pmdec_int = np.concatenate([pmdec_int_b, pmdec_int_f],
                                   axis=axis).flatten()
        rv_int = np.concatenate([rv_int_b, rv_int_f], axis=axis).flatten()

        if (self.df is None):  # If Sun
            time = self.t_mirrored
            group_name = ['Sun'] * len(ra_int)
            ages = [4600] * len(ra_int)

        else:  # For all other traces
            # t = np.tile(self.t_mirrored, len(self.df))
            # age = np.tile(self.df.age.values, len(self.t_mirrored)).T
            # group_name = np.tile(self.df.group_name.values,
            #                      len(self.t_mirrored))
            group_name = []
            ages = []
            time = []
            for name, age in zip(self.df.group_name.values,
                                 self.df.age.values):
                for t in self.t_mirrored:
                    group_name.append(name)
                    ages.append(age)
                    time.append(t)

        df_integrated = pd.DataFrame({
            't': time,
            'ra': ra_int,
            'dec': dec_int,
            'distance': distance_int,
            'pmra': pmra_int,
            'pmdec': pmdec_int,
            'rv': rv_int,
            'age': ages,
            'group_name': group_name
        })

        # adds cartesian coords to dataframe
        df_integrated = convert_icrs_galactic(df_integrated)
        #df_integrated = coordFIX_to_coordROT(df_integrated) # NOTE: converts into rotating coordinate frame
        return df_integrated

    def create_age_based_symbols(self):
        """ 
        """
        self.df_integrated['symbol'] = ['circle'] * len(self.df_integrated)

        if self.marker_properties['size'] != 'size-based':
            self.df_integrated['size'] = [self.marker_properties['size']
                                          ] * len(self.df_integrated)

        if self.marker_behavior == 'disappear':
            self.df_integrated.loc[(self.df_integrated['t'] < 0.) & (
                self.df_integrated['age'] <= self.df_integrated['t'].abs()),
                                   'size'] = .000001

        elif self.marker_behavior == 'circle-open':
            self.df_integrated.loc[(self.df_integrated['t'] < 0.) & (
                self.df_integrated['age'] <= self.df_integrated['t'].abs()),
                                   'symbol'] = 'circle-open'

        elif self.marker_behavior == 'remain':
            pass


class Movie:
    """ Class which encapsulates and controls plotly animation
    """

    def __init__(self,
                 movie_stars: dict,
                 time: np.ndarray,
                 movie_save_path: str,
                 xyz_ranges=[1500, 1500, 400],
                 plot_sun=True,
                 center_star='Sun',
                 camera_follow=True,
                 center_galactic=False,
                 annotations=None,
                 add_extra_traces=True):
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
        self.time_mirror = np.concatenate(
            [np.flip(-1 * self.time[1:]), self.time])
        self.plot_sun = plot_sun
        self.center_star = center_star
        self.annotations = annotations
        self.xyz_ranges = xyz_ranges
        self.camera_follow = camera_follow
        self.center_galactic = center_galactic
        self.add_extra_traces = add_extra_traces

        # read from yaml files
        layout_file = open_yaml('layout_yaml/3d_layout.yaml')
        self.layout_dict = layout_file['3d_galactic_layout']
        self.figure = layout_file['initial_figure']
        self.slider_dict = open_yaml(
            'layout_yaml/sliders.yaml')['3d_galactic_slider']
        self.camera = open_yaml(
            'layout_yaml/cameras.yaml')['3d_galactic_camera']

        # setup initial layout/figure properties
        #self.layout_dict['scene_camera'] = self.camera
        self.figure['layout'] = self.layout_dict

    def make_movie(self):
        '''
        Generates movie of MovieStar intances, centered on the Sun's motion
        '''
        if self.plot_sun:
            self.sun = self.create_sun_orbit()
            self.movie_stars.append(self.sun)

        self.slider_dict['active'] = int(len(self.time_mirror) / 2)

        # generates frames for each timestep and calls 'generate_frame_layout'
        self.generate_frames()

        self.figure['data'] = self.figure['frames'][int(
            len(self.time_mirror) / 2)]['data']  # sets initial figure frame
        self.figure['layout'] = self.figure['frames'][int(
            len(self.time_mirror) /
            2)]['layout']  # sets initial figure frame layout
        self.figure['layout']['scene_camera'] = self.camera  # sets camera
        self.figure['layout']['sliders'] = [self.slider_dict]  # sets slider

        if self.movie_save_path is not None:
            fig = go.Figure(self.figure)
            fig.write_html(self.movie_save_path,
                           auto_open=False,
                           auto_play=False)
        else:
            return self.figure

        print('Movie Created')

    def create_scatter(self, movie_star, df, marker_props, hovertext):

        scatter = go.Scatter3d(x=df.x.values,
                               y=df.y.values,
                               z=df.z.values,
                               mode='markers',
                               marker=marker_props,
                               line=dict(width=0.0001),
                               hovertext=hovertext,
                               name=movie_star.name)
        return scatter

    def generate_frames(self):

        for movie_star in self.movie_stars:
            if movie_star.marker_properties['symbol'] == 'age-based':
                movie_star.create_age_based_symbols()

        for t in self.time_mirror:
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
                    hovertext = df_t['group_name']
                else:
                    hovertext = None

                trace_frame = self.create_scatter(movie_star, df_t,
                                                  marker_props, hovertext)

                frame['data'].append(trace_frame)
                frame['layout'] = self.generate_frame_layout(t)

            if t == 0.:
                if self.add_extra_traces:
                    extra_traces = plotly_extra.extra_traces()
                    for trace in extra_traces:
                        frame['data'].append(trace)

            self.figure['frames'].append(go.Frame(frame))
            slider_step = {
                'args': [[t], {
                    'frame': {
                        'duration': 5,
                        'redraw': True
                    },
                    'mode': 'immediate',
                    'transition': {
                        'duration': 100
                    }
                }],
                'label':
                np.round(t, 1),
                'name':
                'Time',
                'method':
                'animate'
            }
            self.slider_dict['steps'].append(slider_step)

    def generate_frame_layout(self, t):

        x_box_half_size = self.xyz_ranges[0] / 2
        y_box_half_size = self.xyz_ranges[1] / 2
        z_box_half_size = self.xyz_ranges[2] / 2

        largest_box_size = max(self.xyz_ranges)

        layout_annotations = []
        if self.annotations is not None:
            for annotation in self.annotations:  # enter annotations
                layout_annotations.append(annotation)

        if self.camera_follow:
            for ms in self.movie_stars:  # select the stars to be centered on
                if ms.name == self.center_star:
                    df_integrated = ms.df_integrated

            x_min = df_integrated.loc[df_integrated.t ==
                                      t].x.iloc[0] - x_box_half_size
            x_max = df_integrated.loc[df_integrated.t ==
                                      t].x.iloc[0] + x_box_half_size
            y_min = df_integrated.loc[df_integrated.t ==
                                      t].y.iloc[0] - y_box_half_size
            y_max = df_integrated.loc[df_integrated.t ==
                                      t].y.iloc[0] + y_box_half_size
            z_min = -z_box_half_size
            z_max = z_box_half_size
        else:
            if self.center_galactic:
                x_offset = 8300
            else:
                x_offset = 0.
            x_min = -x_box_half_size + x_offset
            x_max = x_box_half_size + x_offset
            y_min = -y_box_half_size
            y_max = y_box_half_size
            z_min = -z_box_half_size
            z_max = z_box_half_size

        # Update for the user-input aspectratio
        layout_dict_t = copy.deepcopy(self.layout_dict)
        layout_dict_t['scene']['aspectratio'][
            'x'] *= x_box_half_size / largest_box_size
        layout_dict_t['scene']['aspectratio'][
            'y'] *= y_box_half_size / largest_box_size
        layout_dict_t['scene']['aspectratio'][
            'z'] *= z_box_half_size / largest_box_size

        # Update for the user-input ranges
        layout_dict_t['scene']['xaxis']['range'] = [x_min, x_max]
        layout_dict_t['scene']['yaxis']['range'] = [y_min, y_max]
        layout_dict_t['scene']['zaxis']['range'] = [z_min, z_max]

        layout_dict_t['scene']['annotations'] = layout_annotations
        layout = go.Layout(layout_dict_t)

        return layout

    def create_sun_orbit(self):

        sun_marker_props = {
            'size': 3,
            'color': 'yellow',
            'opacity': 1.,
            'symbol': 'circle'
        }

        orbit_sun = Orbit([
            0 * u.deg, 0 * u.deg, 0 * u.pc, 0 * u.mas / u.yr, 0 * u.mas / u.yr,
            0 * u.km / u.s
        ],
                          lb=False,
                          uvw=False,
                          radec=True,
                          ro=Galactocentric.galcen_distance,
                          zo=Galactocentric.z_sun,
                          vo=Galactocentric.galcen_v_sun.d_y -
                          GalacticLSR.v_bary.d_y,
                          solarmotion=u.Quantity([
                              -GalacticLSR.v_bary.d_x, GalacticLSR.v_bary.d_y,
                              GalacticLSR.v_bary.d_z
                          ]))

        sun = MovieStar(name='Sun',
                        orbit=orbit_sun,
                        time_to_integrate=self.time,
                        marker_properties=sun_marker_props)

        return sun


# ---- Helper Functions -----
def open_yaml(filepath):
    with open(filepath, 'r') as stream:
        yaml_data = yaml.safe_load(stream)
    return yaml_data


def integrate_orbit(time, orbit, pot, flip=False):
    time = time * u.Myr

    orbit.integrate(time, pot=pot)
    ra = orbit.ra(time)
    dec = orbit.dec(time)
    distance = orbit.dist(time) * 1000
    pmra = orbit.pmra(time)
    pmdec = orbit.pmdec(time)
    rv = orbit.vlos(time)

    if flip:  # flips the time order of the returned ra, dec, distance values
        axis = len(np.shape(ra)) - 1
        ra = np.flip(ra, axis=axis)
        dec = np.flip(dec, axis=axis)
        distance = np.flip(distance, axis=axis)
        pmra = np.flip(pmra, axis=axis)
        pmdec = np.flip(pmdec, axis=axis)
        rv = np.flip(rv, axis=axis)

    return (ra, dec, distance, pmra, pmdec, rv)