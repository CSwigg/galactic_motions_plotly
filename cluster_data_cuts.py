from astropy.io import fits
import pandas as pd
import numpy as np
from astropy.coordinates import ICRS, Galactic, SkyCoord

df_maria = pd.read_csv(
    '/Users/cam/Downloads/maria_catalog_september_with_vlsr.csv')
out_path = '/Users/cam/Desktop/astro_research/radcliffe/galactic_motions_plotly/data/'


def maria_clusters_initial_cuts(df_maria, frac_sample):
    #df_maria = df_maria.loc[df_maria.ref != 'Liu+2019']
    #df_maria = df_maria.loc[df_maria.ref != 'Hao+2022']

    print(len(df_maria))
    df_maria = df_maria.loc[(df_maria.rv_weighted.notnull())
                            & (df_maria.N_tot >= 30)
                            & (df_maria.rv_weighted != 0.) &
                            (df_maria.N_rv_weighted >= 10) &
                            (df_maria.rv_weighted.abs().between(-100, 100))]

    df_maria = df_maria.sample(frac=frac_sample).reset_index(drop=True)

    return df_maria


def maria_rw_joao_selection():
    df_maria_rw = pd.read_csv(
        '/Users/cam/Downloads/RW_clusters_selection4.csv')
    df_maria_rw['age_myr'] = (10**df_maria_rw.age.values)/1e6
    df_maria_rw = df_maria_rw.drop(columns=['age'])
    df_maria_rw.to_csv(out_path + 'maria_clusters_rw_joao.csv', index=False)


def maria_rw_cam_selection(df_maria):

    age_cut = 30  # Myr
    dist_from_model_cut = 150  # pc

    df_rw_model = pd.read_csv(
        '/Users/cam/Downloads/Best_Fit_Wave_Model.csv'
    )
    clusters_coords = SkyCoord(u=df_maria.X.values,
                               v=df_maria.Y.values,
                               w=df_maria.Z.values,
                               frame=Galactic,
                               representation_type='cartesian')
    rw_coords = SkyCoord(u=df_rw_model.x.values,
                         v=df_rw_model.y.values,
                         w=df_rw_model.z.values,
                         frame=Galactic,
                         representation_type='cartesian')
    idx, d2d, d3d = clusters_coords.match_to_catalog_3d(rw_coords)
    df_maria['d_to_rw'] = d3d.value

    df_maria_rw = df_maria.loc[(df_maria['d_to_rw'] <= dist_from_model_cut)
                               & (df_maria['age_myr'] <= age_cut)]
    df_maria_rw.to_csv(out_path + 'maria_clusters_rw_cam.csv', index=False)


def maria_clusters_age_selection(df_maria, age_lower_cut, age_upper_cut):
    df_maria_young = df_maria.loc[df_maria.age_myr <= age_lower_cut]
    df_maria_older = df_maria.loc[df_maria.age_myr.between(age_lower_cut,
                                                           age_upper_cut,
                                                           inclusive='right'
                                                                      )]
    df_maria_young.to_csv(out_path + 'maria_clusters_young.csv', index=False)
    df_maria_older.to_csv(out_path + 'maria_clusters_older.csv', index=False)


def maria_clusters_heiles_shell(df_maria):
    names = [
        'ASCC_32', 'OC-0450', 'CWNU_1162', 'CWNU_1169', 'UPK_473', 'UPK_421',
        'CWNU_1178', 'CWNU_1206', 'CWNU_1162', 'CWNU_1178', 'CWNU_1204',
        'Gulliver_10', 'CWNU_1159', '2288', 'OC-0411', 'OC-0399', 'OC-0410',
        'OC-0407', 'CWNU_1013', 'OC-0408', 'UPK_500', '2395', '2271',
        'IC_2395', 'UPK_523', '2388', 'Collinder_132', 'UPK_477', 'UBC_7',
        'UPK_469', 'UPK_482', 'NGC_2451B', 'UPK_540', 'CWNU_1024', '2436',
        '2386', 'CWNU_1046', 'UPK_478', 'LP_2383', '2383', 'Theia_35',
        'UPK_451', '2397', 'Collinder_135', 'Collinder_140', 'Ruprecht_31',
        'CWNU_1101', 'NGC_2547', 'CWNU_1101', 'CWNU_1082', 'UPK_514',
        'CWNU_528', '2400', 'UPK_496', 'UPK_535', 'CWNU_1144', 'UPK_489',
        'UPK_471', 'Theia_97', '2271', 'OC-0395', 'OC-0406', 'NGC_2362',
        'LP_2270'
    ]
    df_maria_fb = df_maria.loc[df_maria.cluster.isin(names)]
    print(len(df_maria_fb))
    df_maria_fb.to_csv(out_path + 'maria_clusters_heiles_shell.csv',
                       index=False)


def maria_clusters_expanding_group_1(df_maria):
    names = ['CWNU_1236', 'UPK_113', 'UPK_61', 'Theia_98', 'CWNU_1084',
       'UPK_78', 'Stephenson_1', 'CWNU_534', 'UPK_653', 'Theia_116',
       'UPK_92', 'CWNU_1075', 'PHOC_41', 'UPK_83', 'NGC_7058',
       'Alessi_19', 'Col', 'UPK_101', 'CWNU_1113', 'Main', 'Car',
       'Theia_232', 'CWNU_1143', 'CWNU_1032', 'IC_4665', 'Melotte_20',
       'Platais_8', 'chiFor', 'CWNU_1025', 'UPK_37', 'IC_2602', '32Ori',
       'UPK_35', 'UBC_26', 'Collinder_359', 'Alessi_13', 'CWNU_1015',
       'UPK_640', 'UBC_9', 'UPK_632', 'CWNU_1183', 'UPK_622', 'UPK_592',
       'UPK_120']
    df_maria_expanding_1 = df_maria.loc[df_maria.cluster.isin(names)]
    print(len(df_maria_expanding_1))
    df_maria_expanding_1.to_csv(out_path + 'maria_clusters_expanding_1.csv',
                       index=False)

def maria_clusters_expanding_group_2(df_maria):
    names = [
        'Alessi_5', 'UPK_502', 'NGC_3228', 'CWNU_1085', 'ASCC_58', 'Platais_9',
        'UPK_507', 'UPK_481', 'UPK_554', 'UPK_538', 'UPK_479', 'UPK_438',
        'UPK_481', 'Teutsch_38', 'CWNU_1198', 'CWNU_1037', 'CWNU_1039', 'Theia_58',
        'UPK_488', 'CWNU_1102', 'BH_23', 'UPK_515', 'Theia_246', 'CWNU_1016',
        'CWNU_280', 'Trumpler_10', 'UPK_511', 'CWNU_1034', 'CWNU_1044', 'CWNU_515',
        'NGC_2451A', 'UPK_449', 'CWNU_1020', 'OC-0423', 'CWNU_1069', 'IC_2391',
        'CWNU_1134', 'CWNU_1160', 'CWNU_1134', 'CWNU_1078', 'CWNU_1167',
        'CWNU_289', 'BH_164', 'UPK_600', 'CWNU_1078', 'UPK_569', 'CWNU_1008',
        'CWNU_1019', 'UPK_569', 'CWNU_1173', 'CWNU_1016'
    ]
    df_maria_expanding_2 = df_maria.loc[df_maria.cluster.isin(names)]
    print(len(df_maria_expanding_2))
    df_maria_expanding_2.to_csv(out_path + 'maria_clusters_expanding_2.csv',
                       index=False)


def maria_clusters_split_extended(df_maria):
    #df_maria_split = pd.read_csv('/Users/cam/Downloads/maria_young_split_names') # <-- Glue selection

    names = ['Mamajek_1', 'OC-0068', 'UPK_573', 'Main', 'UPK_41', 'Pozzo_1',
       'CWNU_1083', 'OC-0408', 'Upper_Sco', 'CWNU_1261', 'UPK_648',
       'UPK_523', 'CWNU_1048', 'OC-0478', 'OC-0407', 'CWNU_1004',
       'CWNU_1135', 'UPK_36', 'OC-0479', 'LP_2270', 'OC-0062', 'TWHya',
       'UPK_504', 'LP_2435', 'UPK_473', 'CWNU_1096', 'UPK_513',
       'Collinder_132', 'Theia_31', 'UPK_38', 'UPK_519', 'OC-0459',
       'OC-0056', 'CWNU_1136', 'LP_2387', 'OC-0411', 'Distrib',
       'Chameleon_II', 'etaCHA', 'OC-0410', 'Gulliver_10', 'Feigelson_1',
       'OC-0470', 'UPK_503', 'CWNU_1146', 'OC-0666', 'CWNU_1163',
       'OC-0399', 'OC-0451', 'UPK_655', 'OC-0690', 'Haffner_13',
       'CWNU_338', 'OC-0579', 'UPK_474', 'epsCHA', 'PHOC_39', 'farSouth',
       'UPK_39', 'UPK_606', 'CWNU_1224', 'CWNU_1226', 'OC-0578',
       'CWNU_1204']
    df_maria_split = df_maria.loc[df_maria.cluster.isin(names)]
    df_maria_split.to_csv(out_path + 'maria_clusters_split.csv', index = False)




########################################################################################################################

df_maria = maria_clusters_initial_cuts(df_maria, frac_sample=1.)
maria_rw_joao_selection()
maria_rw_cam_selection(df_maria)
maria_clusters_age_selection(df_maria, 25, 70)
maria_clusters_heiles_shell(df_maria)
maria_clusters_expanding_group_1(df_maria)
maria_clusters_expanding_group_2(df_maria)
maria_clusters_split_extended(df_maria)