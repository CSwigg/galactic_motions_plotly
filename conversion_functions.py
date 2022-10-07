import math
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, Galactic, GalacticLSR, Galactocentric
import astropy.units as u

def convert_icrs_galactic(df):
    icrs = ICRS(df.ra.values*u.deg, df.dec.values*u.deg, (df.distance.values)*u.pc)
    gal = icrs.transform_to(Galactic)
    gal.representation_type = 'cartesian'
    df['x'] = gal.u.value
    df['y'] = gal.v.value
    df['z'] = gal.w.value
    return df

def convert_icrs_galactocentric(df):
    icrs = ICRS(df.ra.values*u.deg, df.dec.values*u.deg, (df.distance.values)*u.pc)
    gal = icrs.transform_to(Galactocentric)
    gal.representation_type = 'cartesian'
    df['x_gc'] = gal.x.value
    df['y_gc'] = gal.y.value
    df['z_gc'] = gal.z.value
    return df

def coordFIX_to_coordROT(df_gal, r_sun=8.122, v_sun=220):
    w0 = v_sun / r_sun
    r_sun = r_sun * 1e3  # in pc!
    w1 = w0 / 10
    t1 = df_gal['t'] * 0.01022

    df_gal = convert_icrs_galactocentric(df_gal)    

    r = np.sqrt(df_gal['x_gc']**2 + df_gal['y_gc']**2)
    theta = np.arctan2(df_gal['x_gc'], df_gal['y_gc'])
    # rLSR Cartesian coordinates (rotating LSR)
    x_gal_rot = r_sun - r * np.cos(theta - w1 * t1 + math.pi / 2)
    y_gal_rot = r * np.sin(theta - w1 * t1 + math.pi / 2)
    z_gal_rot = df_gal['z_gc']

    # NOTE: next 3 lines replace the normal Gal x,y,z coordinates in the df
    df_gal['x'] = x_gal_rot
    df_gal['y'] = y_gal_rot
    df_gal['z'] = z_gal_rot

    return(df_gal)