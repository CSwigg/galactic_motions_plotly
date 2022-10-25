import pandas as pd
import numpy as np
from movie_stars import *
import warnings
from astropy.io import fits
import yaml
warnings.filterwarnings("ignore")

time_to_integrate = np.linspace(0., 50., 100)

radcliffe_wave = MovieStar.from_yaml(yaml_file_name = 'radcliffe_wave.yaml', time_to_integrate = time_to_integrate)
maria_young_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_young_clusters.yaml', time_to_integrate = time_to_integrate)
maria_older_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_older_clusters.yaml', time_to_integrate = time_to_integrate)
maria_feedback_clusters = MovieStar.from_yaml(yaml_file_name = 'maria_expanding_clusters.yaml', time_to_integrate = time_to_integrate)
maria_split = MovieStar.from_yaml(yaml_file_name = 'maria_split_clusters.yaml', time_to_integrate = time_to_integrate)
maria_cepheus = MovieStar.from_yaml(yaml_file_name = 'maria_cepheus_clusters.yaml', time_to_integrate = time_to_integrate)

traces_list = [radcliffe_wave, maria_split, maria_cepheus, maria_young_clusters, maria_older_clusters, maria_feedback_clusters]

movie_save_path = '/Users/cam/Desktop/astro_research/radcliffe/movies/gal_mov_test_refactored.html'
movie = Movie(movie_stars=traces_list,
                          time=time_to_integrate,
                          movie_save_path=movie_save_path,
                          center_star='Sun',
                          annotations=None,
                          xyz_ranges=[4000, 4000, 700],
                          camera_follow=True,
                          center_galactic=False,
                          add_extra_traces=True)
movie.make_movie()
