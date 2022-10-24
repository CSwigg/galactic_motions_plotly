import movie_stars
import pandas as pd
import numpy as np

df_apogee = pd.read_csv('/Users/cam/Downloads/apogee_dr17.csv')
time_to_integrate = np.linspace(0, -60, 100)

def traces():
    
    points_disappear = False

    #df_maria_rw = pd.read_csv('/Users/cam/Downloads/maria_rw_clusters.csv')
    df_maria_rw = pd.read_csv('/Users/cam/Downloads/RW_clusters_selection4.csv')
    df_maria_rw['age_myr'] = (10**df_maria_rw.age.values)/1e6
    df_maria_rw = df_maria_rw.drop(columns = ['age'])
    maria_clusters_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'rv_weighted_med': 'radial_velocity',
        'age_myr': 'age',
        'cluster': 'group_name'
    }

    maria_marker_props = {
        'size': 11.,
        'color': '#FF4500',
        'opacity': .7,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }

    maria_clusters_rw = movie_stars.MovieStar(name='Radcliffe Wave',
                                            df=df_maria_rw,
                                            key_dict=maria_clusters_key_dict,
                                            time_to_integrate=time_to_integrate,
                                            marker_properties=maria_marker_props,
                                            disappear=points_disappear,
                                            visible=True)

    df_maria = pd.read_csv(
        '/Users/cam/Downloads/maria_catalog_september_with_vlsr.csv')
    ntot = df_maria.N_tot.values
    #df_maria = df_maria.drop(columns=['radial_velocity', 'age'])
    df_maria = df_maria.sample(frac=1).reset_index(drop=True)
    age_cut_upper = 60
    age_cut_lower = 25
    df_maria = df_maria.loc[df_maria.ref != 'Liu+2019']

    df_maria = df_maria.loc[(df_maria.rv_weighted.notnull())
                            & (df_maria.rv_weighted != 0.) &
                            (df_maria.N_rv_weighted > 5) &
                            (df_maria.rv_weighted.abs().between(-100, 100))]

    df_maria_older = df_maria.loc[(df_maria.age_myr <= age_cut_upper)
                                & (df_maria.age_myr > age_cut_lower)]
    df_maria_older.to_csv('/Users/cam/Downloads/maria_older.csv', index = False)
    maria_marker_props = {
        'size': 5.,
        'color': 'cyan',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }
    maria_clusters = movie_stars.MovieStar(
        name='Maria Catalog: {} Myr < Age <= {} Myr'.format(
            age_cut_lower, age_cut_upper),
        df=df_maria_older,
        key_dict=maria_clusters_key_dict,
        time_to_integrate=time_to_integrate,
        marker_properties=maria_marker_props,
        disappear=points_disappear,
        visible=True)
    
    df_maria_young = df_maria.loc[df_maria.age_myr <= age_cut_lower]
    maria_marker_props_copy = maria_marker_props.copy()
    maria_marker_props_copy['size'] = 5.
    maria_marker_props_copy['color'] = 'red'
    maria_clusters_young = movie_stars.MovieStar(name='Maria Catalog: Age <= {} Myr'.format(age_cut_lower),
                                       df=df_maria_young,
                                       key_dict=maria_clusters_key_dict,
                                       time_to_integrate=time_to_integrate,
                                       marker_properties=maria_marker_props_copy,
                                       disappear=points_disappear,
                                       visible=True)
    names = [
    'ASCC_32', 'OC-0450', 'CWNU_1162', 'CWNU_1169', 'UPK_473', 'UPK_421',
    'CWNU_1178', 'CWNU_1206', 'CWNU_1162', 'CWNU_1178', 'CWNU_1204',
    'Gulliver_10', 'CWNU_1159', '2288', 'OC-0411', 'OC-0399', 'OC-0410',
    'OC-0407', 'CWNU_1013', 'OC-0408', 'UPK_500', '2395', '2271', 'IC_2395',
    'UPK_523', '2388', 'Collinder_132', 'UPK_477', 'UBC_7', 'UPK_469',
    'UPK_482', 'NGC_2451B', 'UPK_540', 'CWNU_1024', '2436', '2386',
    'CWNU_1046', 'UPK_478', 'LP_2383', '2383', 'Theia_35', 'UPK_451', '2397',
    'Collinder_135', 'Collinder_140', 'Ruprecht_31', 'CWNU_1101', 'NGC_2547',
    'CWNU_1101', 'CWNU_1082', 'UPK_514', 'CWNU_528', '2400', 'UPK_496',
    'UPK_535', 'CWNU_1144', 'UPK_489', 'UPK_471', 'Theia_97', '2271',
    'OC-0395', 'OC-0406', 'NGC_2362', 'LP_2270'
]
    feature_name = 'Expanding Clusters'

    #names = ['NGC_2547']
    #feature_name = 'NGC_2547'

    df_maria_fb = df_maria.loc[df_maria.cluster.isin(names)]
    df_maria_fb.to_csv('/Users/cam/Downloads/maria_clusters_feedback.csv', index=False)
    maria_marker_props_fb = {
        'size': 5.,
        'color': '#7B68EE',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }

    maria_clusters_fb = movie_stars.MovieStar(
        name=feature_name,
        df=df_maria_fb,
        key_dict=maria_clusters_key_dict,
        time_to_integrate=time_to_integrate,
        marker_properties=maria_marker_props_fb,
        disappear=points_disappear,
        visible=True)


    data_list = [
        maria_clusters, maria_clusters_rw
    ]
    return (data_list, time_to_integrate)