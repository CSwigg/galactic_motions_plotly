import movie_stars
import pandas as pd
import numpy as np

df_apogee = pd.read_csv('/Users/cam/Downloads/apogee_dr17.csv')
time_to_integrate = np.linspace(0, -60, 100)

def traces():
    ### Sigma ScoCen
    df_sc = pd.read_csv('/Users/cam/Downloads/sc_sebastian_full_new.csv')
    df_sc = df_sc.loc[df_sc.n_rv > 10]
    sco_cen_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'rv_chosen': 'radial_velocity',
        'label': 'group_name',
        'age': 'age'
    }
    sco_cen_marker_props = {
        'size': 5,
        'opacity': 1.,
        'symbol': 'age-based',
        'color': '#00FA9A',
        'line': dict(width=0.)
    }

    sco_cen = movie_stars.MovieStar(name='SiGMA Sco-Cen',
                                    time_to_integrate=time_to_integrate,
                                    df=df_sc,
                                    key_dict=sco_cen_key_dict,
                                    marker_properties=sco_cen_marker_props,
                                    disappear=True)



    ### Kerr+2021
    df_kerr = pd.read_csv('../kerr_sos_edr3.csv')
    df_kerr = df_kerr.rename(columns={'EOM': 'eom'})
    #(df_kerr.tlc != 22)
    #df_kerr = df_kerr.loc[(df_kerr.eom >= 0)]
    df_kerr = pd.merge(left=df_kerr,
                    right=df_apogee,
                    how='left',
                    left_on='source_id',
                    right_on='gaiaedr3_source_id')
    df_kerr = df_kerr.loc[(df_kerr.parallax > 0.)]

    df_kerr_apogee = df_kerr.loc[df_kerr['vhelio_avg'].notnull()]
    df_kerr_apogee['rv_best'] = df_kerr_apogee['vhelio_avg']
    df_kerr_apogee['rv_best_err'] = df_kerr_apogee['verr']

    # df_kerr_gaia = df_kerr.loc[(df_kerr['vhelio_avg'].isnull()) & (df_kerr['radial_velocity'].notnull())]
    # df_kerr_gaia['rv_best'] = df_kerr_gaia['radial_velocity']

    df_kerr_sos = df_kerr.loc[(df_kerr['RVcor_merged'].notnull())
                            & (df_kerr['vhelio_avg'].isnull())]
    df_kerr_sos['rv_best'] = df_kerr_sos['RVcor_merged']
    df_kerr_sos['rv_best_err'] = df_kerr_sos['ERRVcor_merged']

    df_kerr = pd.concat([df_kerr_apogee, df_kerr_sos])
    df_kerr['rv_best_rel_err'] = np.abs(df_kerr['rv_best_err'].values /
                                        df_kerr['rv_best'].values)
    df_kerr = df_kerr.loc[df_kerr['rv_best_rel_err'] < .3]

    #df_kerr = df_kerr.drop(columns = ['radial_velocity'])

    df_kerr = df_kerr.groupby(['eom']).median().reset_index()
    df_kerr = df_kerr.loc[df_kerr.Age.notnull()]

    kerr_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'rv_best': 'radial_velocity',
        'eom': 'group_name',
        'Age': 'age'
    }
    kerr_marker_props = {
        'size': 5,
        'opacity': 1.,
        'symbol': 'age-based',
        'color': '#4169E1',
        'line': dict(width=0.)
    }

    kerr = movie_stars.MovieStar(name='Kerr+2021',
                                time_to_integrate=time_to_integrate,
                                df=df_kerr,
                                key_dict=kerr_key_dict,
                                marker_properties=kerr_marker_props,
                                disappear=True)


    df_nuria = pd.read_csv('/Users/cam/Downloads/table_YLAs_equ_family2.csv')
    df_nuria['age'] = np.ones(len(df_nuria)) * 35

    nuria_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'RV': 'radial_velocity',
        'age': 'age',
        'Group': 'group_name'
    }
    nuria_marker_props = {
        'size': 5,
        'color': 'orange',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }

    nuria = movie_stars.MovieStar(name='Nuria Moving Groups',
                                time_to_integrate=time_to_integrate,
                                df=df_nuria,
                                key_dict=nuria_key_dict,
                                marker_properties=nuria_marker_props)


    df_orion = pd.read_csv(
        '/Users/cam/Desktop/astro_research/orion/orion_full_groups_with_age.csv')
    df_orion = df_orion.loc[df_orion.apogee2_radial_velocity.notnull()]
    df_orion = df_orion.groupby(['label']).median().reset_index()
    df_orion = df_orion.loc[df_orion.probability > .2]

    orion_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'apogee2_radial_velocity': 'radial_velocity',
        'age': 'age',
        'label': 'group_name'
    }
    orion_marker_props = {
        'size': 7,
        'color': '#FF8C00',
        'opacity': .7,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }
    orion = movie_stars.MovieStar(name='Swiggum+21 Orion',
                                time_to_integrate=time_to_integrate,
                                df=df_orion,
                                key_dict=orion_key_dict,
                                marker_properties=orion_marker_props,
                                disappear=True)
    orion.df_integrated.to_csv('../clusters_integrated/orion.csv', index=False)



    df_maria_rw = pd.read_csv('/Users/cam/Downloads/maria_rw_clusters.csv')
    df_maria_rw = df_maria_rw.loc[(df_maria_rw.age_myr <= 50)
                                & (df_maria_rw.rv_weighted.notnull()) &
                                (df_maria_rw.rv_weighted != 0.) &
                                (df_maria_rw.N_rv_weighted > 3) &
                                (df_maria_rw.rv_weighted.abs().between(-50, 50))]
    df_maria_rw = df_maria_rw.drop(columns=['radial_velocity', 'age'])

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
        'size': 4.,
        'color': 'red',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }

    maria_clusters_rw = movie_stars.MovieStar(name='Maria Catalog: Radcliffe Wave',
                                            df=df_maria_rw,
                                            key_dict=maria_clusters_key_dict,
                                            time_to_integrate=time_to_integrate,
                                            marker_properties=maria_marker_props,
                                            disappear=True,
                                            visible=True)
    maria_clusters_rw.df_integrated.to_csv(
        '../clusters_integrated/maria_clusters_rw.csv', index=False)


    df_maria = pd.read_csv('/Users/cam/Downloads/maria_clusters_new.csv')
    df_maria = df_maria.loc[df_maria.ref != 'Liu+2019']
    df_maria = df_maria.loc[(df_maria.age_myr <= 50) &
                            (df_maria.rv_weighted.notnull()) &
                            (df_maria.rv_weighted != 0.) &
                            (df_maria.N_rv_weighted > 3) &
                            (df_maria.rv_weighted.abs().between(-100, 100))
                            #(1000/df_maria.parallax.values <= 2500)
                            ]
    df_maria = df_maria.drop(columns=['radial_velocity', 'age'])
    df_maria = df_maria.sample(frac=1).reset_index(drop=True)

    maria_marker_props = {
        'size': 3.,
        'color': 'cyan',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }

    maria_clusters = movie_stars.MovieStar(name='Maria Catalog < 50 Myr',
                                        df=df_maria,
                                        key_dict=maria_clusters_key_dict,
                                        time_to_integrate=time_to_integrate,
                                        marker_properties=maria_marker_props,
                                        disappear=True,
                                        visible=True)
    maria_clusters.df_integrated.to_csv(
        '../clusters_integrated/maria_clusters.csv', index=False)


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
        'UPK_535', 'CWNU_1144', 'UPK_489', 'UPK_471', 'Theia_97', '2271', 'OC-0395'
    ]
    feature_name = 'Maria Catalog: Extended Vela Feedback Event'

    #names = ['NGC_2547']
    #feature_name = 'NGC_2547'

    df_maria_fb = df_maria.loc[df_maria.cluster.isin(names)]

    maria_marker_props_fb = {
        'size': 5.,
        'color': 'magenta',
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
        disappear=True,
        visible=True)



    df_laura = pd.read_csv('/Users/cam/Downloads/CamerenTable_CrASubgroups.csv')

    laura_marker_props = {
        'size': 7.,
        'color': '#B2AC88',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }

    laura_clusters_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'radial_velocity': 'radial_velocity',
        'age': 'age',
        'name': 'group_name'
    }

    laura_clusters = movie_stars.MovieStar(name='Laura Corona Australis',
                                        df=df_laura,
                                        key_dict=laura_clusters_key_dict,
                                        time_to_integrate=time_to_integrate,
                                        marker_properties=laura_marker_props,
                                        disappear=False)



    df_elias = pd.read_csv('/Users/cam/Downloads/quillen.csv')
    elias_props = {
        'size': 5.,
        'color': '#66CDAA',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }
    elias_clusters_key_dict = {
        'ra': 'ra',
        'dec': 'dec',
        'parallax': 'parallax',
        'pmra': 'pmra',
        'pmdec': 'pmdec',
        'radial_velocity': 'radial_velocity',
        'age_myr': 'age',
        'Name': 'group_name'
    }

    elias = movie_stars.MovieStar(name='Quillen+2020 OCs',
                                df=df_elias,
                                key_dict=elias_clusters_key_dict,
                                time_to_integrate=time_to_integrate,
                                marker_properties=elias_props,
                                disappear=True)



    df_cg_vela = pd.read_csv('/Users/cam/Downloads/cantat_2019a_vela_grouped.csv')
    cg_vela_props = {
        'size': 5.,
        'color': '#FFB6C1',
        'opacity': 1.,
        'symbol': 'age-based',
        'line': dict(width=0.)
    }
    cg_vela_clusters_key_dict = {
        'RA_ICRS': 'ra',
        'DE_ICRS': 'dec',
        'Plx': 'parallax',
        'pmRA': 'pmra',
        'pmDE': 'pmdec',
        'RV': 'radial_velocity',
        'age_myr': 'age',
        'Pop': 'group_name'
    }

    cg_vela = movie_stars.MovieStar(name='CG+2019a Vela',
                                    df=df_cg_vela,
                                    key_dict=cg_vela_clusters_key_dict,
                                    time_to_integrate=time_to_integrate,
                                    marker_properties=cg_vela_props,
                                    disappear=True)

    # df_hmxb = pd.read_csv('/Users/cam/Downloads/hmxb_google.csv')
    # df_hmxb = df_hmxb.loc[df_hmxb.age_sn != -9999]
    # df_hmxb['parallax'] = 1/df_hmxb.distance
    
    # hmxb_props = {
    #     'size': 7.,
    #     'color': '#FFB6C1',
    #     'opacity': 1.,
    #     'symbol': 'age-based',
    #     'line': dict(width=0.)
    # }
    # hmxb_key_dict = {
    #     'ra': 'ra',
    #     'dec': 'dec',
    #     'parallax': 'parallax',
    #     'pmra': 'pmra',
    #     'pmdec': 'pmdec',
    #     'rv': 'radial_velocity',
    #     'age_sn': 'age',
    #     'hmxb': 'group_name'
    # }
    
    # hmxb = movie_stars.MovieStar(name='HMXBs',
    #                                 df=df_hmxb,
    #                                 key_dict=hmxb_key_dict,
    #                                 time_to_integrate=time_to_integrate,
    #                                 marker_properties=hmxb_props,
    #                                 disappear=points_disappear)
    

    data_list = [
        maria_clusters, maria_clusters_rw, maria_clusters_fb, sco_cen, kerr, orion,
        laura_clusters, nuria, elias, cg_vela
    ]
    return (data_list, time_to_integrate)