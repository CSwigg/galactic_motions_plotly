import pandas as pd
import numpy as np
import movie_stars

time_to_integrate = np.linspace(0, -30, 200)

# ------------ Open Clusters --------------
df_open_clusters = pd.read_csv('/Users/cam/Downloads/open_clusters_3d.csv')
df_open_clusters['age_myr'] = df_open_clusters['age_myr']/1e6
df_oc_young = df_open_clusters.loc[df_open_clusters.age_myr <= 20.]
df_oc_old = df_open_clusters.loc[(df_open_clusters.age_myr > 20.) & (df_open_clusters.age_myr < 40)]

open_clusters_key_dict = {'ra':'ra', 
                          'dec':'dec', 
                          'parallax':'parallax', 
                          'pmra':'pmra',
                          'pmdec':'pmdec',
                          'rv':'radial_velocity',
                          'age_myr':'age', 
                          'cluster':'group_name'}
oc_marker_props = {'size' : 4, 
                             'color' : 'cyan',
                             'opacity' : 0.8, 
                             'symbol' : 'circle'
                             }
oc_marker_props_old = oc_marker_props.copy()
oc_marker_props_old['color'] = 'red'

open_clusters_young = movie_stars.MovieStar(name = 'Open Clusters < 20 Myr', 
                                      df = df_oc_young, 
                                      key_dict = open_clusters_key_dict,
                                      time_to_integrate = time_to_integrate,
                                      marker_properties=oc_marker_props
                                      )

open_clusters_old = movie_stars.MovieStar(name = 'Open Clusters > 20 Myr', 
                                      df = df_oc_old, 
                                      key_dict = open_clusters_key_dict,
                                      time_to_integrate = time_to_integrate,
                                      marker_properties=oc_marker_props_old
                                      )




# ------------ Kerr+2021 --------------


movie_save_path = '/Users/cam/Desktop/astro_research/radcliffe/galactic_motions_movie_1.html'
data_list = [open_clusters_young, open_clusters_old]
movie = movie_stars.Movie(movie_stars = data_list, time = time_to_integrate, movie_save_path=movie_save_path)
movie.make_movie()