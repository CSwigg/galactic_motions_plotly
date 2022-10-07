import numpy as np
import math
import pandas as pd


def coordFIX_to_coordROT(df_gal, r_sun=8.5, v_sun=220):
    w0 = v_sun / r_sun
    r_sun = r_sun * 1e3  # in pc!
    w1 = w0 / 10
    t1 = df_gal['t'] * 0.01022

    r = np.sqrt(df_gal['X_gal']**2 + df_gal['Y_gal']**2)
    theta = np.arctan2(df_gal['X_gal'], df_gal['Y_gal'])
    # rLSR Cartesian coordinates (rotating LSR)
    x_gal_rot = r_sun - r * np.cos(theta - w1 * t1 + math.pi / 2)
    y_gal_rot = r * np.sin(theta - w1 * t1 + math.pi / 2)
    z_gal_rot = df_gal['Z_gal']

    df_gal['X_gal_rot'] = x_gal_rot
    df_gal['Y_gal_rot'] = y_gal_rot
    df_gal['Z_gal_rot'] = z_gal_rot

    return(df_gal)


# def carteFIX_to_carteROT(
#         carte_gal,
#         gal_labels,
#         r_sun=8.5,
#         v_sun=220,
#         carte_labels=['x_carte_rLSR', 'y_carte_rLSR', 'z_carte_rLSR']):
#     '''
#     From Cartesian non-rotating in kpc and km/s
#     returns Cartesian rotating in kpc and km/s
#     INPUT:
#     carte_gal : dataframe with galactic cartesian coordinates
#     gal_labels : column names of the coordinates. IMPORTANT: ORDER MATTERS x, y, z, time
#     '''

#     # constants
#     w0 = v_sun / r_sun
#     r_sun = r_sun * 1e3  # in pc!
#     w1 = w0 / 10
#     t1 = carte_gal[gal_labels[3]] * 0.01022
#     #######################################################################
#     #   POSITIONS
#     #######################################################################
#     # Cilindric Galactic
#     r = np.sqrt(carte_gal[gal_labels[0]]**2 + carte_gal[gal_labels[1]]**2)
#     theta = np.arctan2(carte_gal[gal_labels[0]], carte_gal[gal_labels[1]])
#     # rLSR Cartesian coordinates (rotating LSR)
#     x_carte = r_sun - r * np.cos(theta - w1 * t1 + math.pi / 2)
#     y_carte = r * np.sin(theta - w1 * t1 + math.pi / 2)
#     z_carte = carte_gal[gal_labels[2]]
#     # Output dataframes
#     carte_coords = pd.DataFrame(data=np.stack([x_carte, y_carte, z_carte]).T,
#                                 columns=carte_labels)
#     return carte_coords