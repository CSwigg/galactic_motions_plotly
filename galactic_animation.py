import pandas as pd
import numpy as np
from movie_stars import *
from astropy.io import fits
import yaml
import warnings
warnings.filterwarnings("ignore")

time_to_integrate = np.linspace(0., 60, 100)

radcliffe_wave = MovieStar.from_yaml(yaml_file_name = 'radcliffe_wave.yaml', time_to_integrate = time_to_integrate)
maria_young_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_young_clusters.yaml', time_to_integrate = time_to_integrate)
maria_older_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_older_clusters.yaml', time_to_integrate = time_to_integrate)
#maria_old_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_oldest_clusters.yaml', time_to_integrate = time_to_integrate)
maria_heiles_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_heiles_clusters.yaml', time_to_integrate = time_to_integrate)
maria_expanding_clusters_2 = MovieStar.from_yaml(yaml_file_name = 'maria_expanding_clusters_2.yaml', time_to_integrate = time_to_integrate)
maria_expanding_clusters_1 = MovieStar.from_yaml(yaml_file_name = 'maria_expanding_clusters_1.yaml', time_to_integrate = time_to_integrate)
maria_split = MovieStar.from_yaml(yaml_file_name = 'maria_split_clusters.yaml', time_to_integrate = time_to_integrate)
maria_cepheus = MovieStar.from_yaml(yaml_file_name = 'maria_cepheus_clusters.yaml', time_to_integrate = time_to_integrate)
quillen = MovieStar.from_yaml(yaml_file_name = 'quillen.yaml', time_to_integrate = time_to_integrate)

# df_integrated1 = maria_expanding_clusters_1.df_integrated
# df_integrated2 = maria_expanding_clusters_2.df_integrated
# df_integrated_heiles = maria_heiles_clusters.df_integrated
# df_integrated_rw = radcliffe_wave.df_integrated
# df_integrated_split = maria_split.df_integrated
# df_integrated_cepheus_spur = maria_cepheus.df_integrated
# df_integrated_older= maria_older_clusters.df_integrated
# df_integrated_young= maria_young_clusters.df_integrated


# df_integrated1.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/maria_expanding_clusters_1_int.csv', index = False)
# df_integrated2.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/maria_expanding_clusters_2_int.csv', index = False)
# df_integrated_heiles.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/maria_heiles_clusters_int.csv', index = False)
# df_integrated_rw.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/rw_int.csv', index = False)
# df_integrated_split.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/split_int.csv', index = False)
# df_integrated_cepheus_spur.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/cepheus_spur_int.csv', index = False)
#df_integrated_older.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/maria_older_int.csv', index = False)
#df_integrated_young.to_csv('/Users/cam/Desktop/astro_research/radcliffe/clusters_integrated/maria_young_int.csv', index = False)

#traces_list = [radcliffe_wave, maria_split, maria_cepheus, maria_young_clusters, maria_older_clusters, maria_heiles_clusters, maria_expanding_clusters_1, maria_expanding_clusters_2]
#traces_list = [radcliffe_wave, maria_split, maria_heiles_clusters, maria_expanding_clusters_1, maria_expanding_clusters_2]
traces_list = [radcliffe_wave, maria_split, maria_young_clusters, maria_older_clusters, maria_heiles_clusters, maria_expanding_clusters_1, maria_expanding_clusters_2, quillen]

annotations = [dict(x=-62, y=-335, z=-17, text='Vela OB2 (t=0)', xanchor='center', opacity=1., font=dict(size=12, color='white'), showarrow=True),
               dict(x=125, y=-46, z=47, text='Sco-Cen (t=0)', xanchor='center',
                    opacity=1., font=dict(size=12), showarrow=True),
               ]


movie_save_path = '/Users/cam/Desktop/astro_research/radcliffe/movies/gal_mov_large.html'
movie = Movie(movie_stars=traces_list,
                          time=time_to_integrate,
                          movie_save_path=movie_save_path,
                          center_star='Sun',
                          annotations=None,
                          xyz_ranges=[25000, 25000, 600],
                          camera_follow=False,
                          center_galactic=True,
                          add_extra_traces=True
                          )
movie.make_movie()
