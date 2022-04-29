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