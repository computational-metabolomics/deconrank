import collections

adduct_groups = collections.OrderedDict([('1 165.078',
                                          {'11215': ['[4M+K]+', '11215', '[866][M]+'],
                                           '3': ['[M+H-NH3-HCOOH]+', '3', ''],
                                           '328': ['[M+H]+', '328', ''],
                                           '3977': ['[2M+H]+', '3977', '[340][M]+'],
                                           '4695': ['[2M+Na]+', '4695', '[380][M]+'],
                                           '508': ['[M+Na]+', '508', '[52][M]+'],
                                           '5125': ['[2M+K]+', '5125', ''],
                                           '52': ['[M+H-HCOOH]+', '52', ''],
                                           '5297': ['[2M+2Na-H]+', '5297', ''],
                                           '704': ['[M+K]+', '704', '[69][M]+'],
                                           '781': ['[M+2Na-H]+', '781', '[81][M]+'],
                                           '8596': ['[3M+Na]+', '8596', '[629][M]+'],
                                           '8934': ['[3M+K]+', '8934', '[662][M]+'],
                                           '9044': ['[3M+2Na-H]+', '9044', '[673][M]+'],
                                           '92': ['[M+H-NH3-H2O]+', '92', '[8][M]+']}),
                                         ('1 103.047',
                                          {'340': ['[M+NH4+HCOOH]+', '340', ''],
                                           '343': ['[M+NH4+HCOOH]+', '343', '[35][M]+'],
                                           '57': ['[M+H+NH3]+', '57', ''],
                                           '58': ['[M+H+NH3]+', '58', '[4][M]+'],
                                           '6': ['[M+H]+', '6', ''],
                                           '6265': ['[4M+H]+', '6265', '']}),
                                         ('1 305.155',
                                          {'10373': ['[4M+2H]2+', '10373', ''],
                                           '10663': ['[2M+Na]+', '10663', ''],
                                           '3185': ['[M+H]+', '3185', ''],
                                           '841': ['[M+H-NH3-CO2-CH2O]+', '841', '']}),
                                         ('1 747.189',
                                          {'11649': ['[M+H+NH3]+', '11649', ''],
                                           '1510': ['[M+3H]3+', '1510', '']}),
                                         ('1 323.068',
                                          {'1880': ['[M+H-C2H4O2]+', '1880', '[176][M]+'],
                                           '1974': ['[M+H-C3H4O]+', '1974', ''],
                                           '3756': ['[M+H]+', '3756', '']}),
                                         ('1 566.129',
                                          {'5747': ['[M+H-C6H8O6]+', '5747', ''],
                                           '9929': ['[M+H+NH3]+', '9929', '[756][M]+']}),
                                         ('1 443.156',
                                          {'6617': ['[M+H-CH3]+', '6617', ''],
                                           '6984': ['[M+H]+', '6984', '']}),
                                         ('2 343.124',
                                          {'1469': ['[M+H-C4H6-COCH2]+', '1469', ''],
                                           '1840': ['[M+H-(H2O)3-CO]+', '1840', ''],
                                           '2648': ['[M+H-(H2O)3]+', '2648', ''],
                                           '2907': ['[M+H-HCOOH]+', '2907', ''],
                                           '4919': ['[M+H+NH3]+', '4919', '[392][M]+'],
                                           '929': ['[M+H-C4H6-C4H6O]+', '929', '']}),
                                         ('2 334.148',
                                          {'3609': ['[M+H-O]+', '3609', ''],
                                           '4112': ['[M+H]+', '4112', '']}),
                                         ('2 321.025',
                                          {'4223': ['[M+H+NH3]+', '4223', ''],
                                           '4904': ['[M+Na+NH3]+', '4904', '']}),
                                         ('2 338.052',
                                          {'4223': ['[M+H]+', '4223', ''],
                                           '4904': ['[M+Na]+', '4904', '']}),
                                         ('2 443.251',
                                          {'4434': ['[M+H-NH3-C3H4-COCH2]+', '4434', ''],
                                           '6994': ['[M+H]+', '6994', ''],
                                           '7529': ['[M+Na]+', '7529', ''],
                                           '7839': ['[M+K]+', '7839', '']}),
                                         ('3 181.073',
                                          {'120': ['[M+H-HCOOH]+', '120', ''],
                                           '173': ['[M+H-NH3-H2O]+', '173', '[18][M]+'],
                                           '316': ['[M+H-NH3]+', '316', '[32][M]+'],
                                           '39': ['[M+H-H2O-HCOOH]+', '39', ''],
                                           '45': ['[M+H-NH3-HCOOH]+', '45', '[2][M]+'],
                                           '464': ['[M+H]+', '464', '[47][M]+'],
                                           '64': ['[M+H-NH3-COCH2]+', '64', '[5][M]+'],
                                           '706': ['[M+Na]+', '706', '']}),
                                         ('3 274.209',
                                          {'1822': ['[M+H-CH2]+', '1822', ''],
                                           '83': ['[M+2H-H20]2+', '83', '']}),
                                         ('3 119.042',
                                          {'1279': ['[4M+2H]2+', '1279', ''],
                                           '128': ['[M+H+NH3]+', '128', ''],
                                           '129': ['[M+H+NH3]+', '129', '[12][M]+']}),
                                         ('3 317.077', {'1263': ['[M+H-H2O-H2O-C2H4O', '1263', '']}),
                                         ('3 246.039', {'1882': ['[M+H+NH3]+', '1882', '[177][M]+']}),
                                         ('3 343.107',
                                          {'1882': ['[M+H-H2O-H2O-C2H4O', '1882', '[177][M]+']})])

isotope_groups = collections.OrderedDict([('[4]', {'58': ['[4][M]+', '58'], '61': ['[4][M+1]+', '61']}),
                                          ('[8]', {'92': ['[8][M]+', '92'], '95': ['[8][M+1]+', '95']}),
                                          ('[35]',
                                           {'343': ['[35][M]+', '343'],
                                            '349': ['[35][M+1]+', '349'],
                                            '354': ['[35][M+2]+', '354']}),
                                          ('[52]',
                                           {'508': ['[52][M]+', '508'], '519': ['[52][M+1]+', '519']}),
                                          ('[69]',
                                           {'704': ['[69][M]+', '704'],
                                            '718': ['[69][M+1]+', '718'],
                                            '736': ['[69][M+2]+', '736']}),
                                          ('[81]',
                                           {'781': ['[81][M]+', '781'], '790': ['[81][M+1]+', '790']}),
                                          ('[176]',
                                           {'1880': ['[176][M]+', '1880'],
                                            '1900': ['[176][M+1]+', '1900']}),
                                          ('[340]',
                                           {'3977': ['[340][M]+', '3977'],
                                            '4013': ['[340][M+1]+', '4013'],
                                            '4056': ['[340][M+2]+', '4056']}),
                                          ('[380]',
                                           {'4695': ['[380][M]+', '4695'],
                                            '4712': ['[380][M+1]+', '4712']}),
                                          ('[448]',
                                           {'5883': ['[448][M]+', '5883'],
                                            '5908': ['[448][M+1]+', '5908']}),
                                          ('[534]',
                                           {'7243': ['[534][M]2+', '7243'],
                                            '7261': ['[534][M+1]2+', '7261']}),
                                          ('[625]',
                                           {'8544': ['[625][M]2+', '8544'],
                                            '8554': ['[625][M+1]2+', '8554']}),
                                          ('[629]',
                                           {'8596': ['[629][M]+', '8596'],
                                            '8618': ['[629][M+1]+', '8618']}),
                                          ('[662]',
                                           {'8934': ['[662][M]+', '8934'],
                                            '8965': ['[662][M+1]+', '8965'],
                                            '8975': ['[662][M+2]+', '8975']}),
                                          ('[673]',
                                           {'9044': ['[673][M]+', '9044'],
                                            '9073': ['[673][M+1]+', '9073']}),
                                          ('[680]',
                                           {'9143': ['[680][M]2+', '9143'],
                                            '9155': ['[680][M+1]2+', '9155']}),
                                          ('[702]',
                                           {'9376': ['[702][M]+', '9376'],
                                            '9400': ['[702][M+1]+', '9400']}),
                                          ('[748]',
                                           {'9837': ['[748][M]+', '9837'],
                                            '9860': ['[748][M+1]+', '9860']}),
                                          ('[756]',
                                           {'9929': ['[756][M]+', '9929'],
                                            '9946': ['[756][M+1]+', '9946']}),
                                          ('[866]',
                                           {'11215': ['[866][M]+', '11215'],
                                            '11226': ['[866][M+1]+', '11226']}),
                                          ('[392]',
                                           {'4919': ['[392][M]+', '4919'],
                                            '4946': ['[392][M+1]+', '4946']}),
                                          ('[758]',
                                           {'9947': ['[758][M]+', '9947'],
                                            '9962': ['[758][M+1]+', '9962'],
                                            '9982': ['[758][M+2]+', '9982']}),
                                          ('[2]', {'45': ['[2][M]+', '45'], '50': ['[2][M+1]+', '50']}),
                                          ('[5]', {'64': ['[5][M]+', '64'], '67': ['[5][M+1]+', '67']}),
                                          ('[12]',
                                           {'129': ['[12][M]+', '129'], '138': ['[12][M+1]+', '138']}),
                                          ('[18]',
                                           {'173': ['[18][M]+', '173'], '180': ['[18][M+1]+', '180']}),
                                          ('[32]',
                                           {'316': ['[32][M]+', '316'],
                                            '324': ['[32][M+1]+', '324'],
                                            '338': ['[32][M+2]+', '338']}),
                                          ('[47]',
                                           {'464': ['[47][M]+', '464'], '473': ['[47][M+1]+', '473']}),
                                          ('[177]',
                                           {'1882': ['[177][M]+', '1882'],
                                            '1904': ['[177][M+1]+', '1904']})])


header = ['peakID', 'grpid', 'mzmed', 'mzmin', 'mzmax', 'rtmed', 'rtmin', 'rtmax', 'npeaks', 'blank', 'lcms', 'lcmsms', 'Blank_P_WAX_1_LCMS-FC_Phenyl_pos_inj1', 'Blank_P_WAX_1_LCMS-FC_Phenyl_pos_inj2', 'Blank_P_WAX_1_LCMS-FC_Phenyl_pos_inj3', 'Blank_P_WAX_1_LCMS-FC_Phenyl_pos_inj4', 'Blank_P_WAX_1_LCMS-FC_Phenyl_pos_inj5', 'Daph_P_WAX_1_LCMS-FC_Phenyl_pos_inj1', 'Daph_P_WAX_1_LCMS-FC_Phenyl_pos_inj2', 'Daph_P_WAX_1_LCMS-FC_Phenyl_pos_inj3', 'Daph_P_WAX_1_LCMS-FC_Phenyl_pos_inj4', 'Daph_P_WAX_1_LCMS-FC_Phenyl_pos_inj5', 'Daph_P_WAX_1_LCMSMS-FC_Phenyl_pos_inj1__dynamicExclu_NoLists', 'Daph_P_WAX_1_LCMSMS-FC_Phenyl_pos_inj2__incl2__1', 'Daph_P_WAX_1_LCMSMS-FC_Phenyl_pos_inj3__incl2___2', 'Daph_P_WAX_1_LCMSMS-FC_Phenyl_pos_inj4__incl1__', 'blank_median_I', 'lcms_median_I', 'lcmsms_median_I', 'sample_median_I', 'blank_RSD_I', 'lcms_RSD_I', 'lcmsms_RSD_I', 'blank_coverage', 'lcms_coverage', 'lcmsms_coverage', 'blank_RSD_RT', 'lcms_RSD_RT', 'lcmsms_RSD_RT', 'rsd_all_RT', 'blank_valid', 'lcms_valid', 'lcmsms_valid', 'all_sample_valid', 'isotopes', 'adduct', 'pcgroup', 'msms_count']


grouped_features = collections.OrderedDict([(1,
              ['9073',
               '704',
               '5297',
               '11226',
               '8934',
               '4695',
               '790',
               '3',
               '508',
               '736',
               '8596',
               '328',
               '8618',
               '8965',
               '3977',
               '718',
               '11215',
               '5125',
               '92',
               '95',
               '9044',
               '4056',
               '781',
               '4712',
               '4013',
               '52',
               '519',
               '8975']),
             (2, ['6265', '58', '57', '6', '61', '354', '340', '343', '349']),
             (3, ['41']),
             (4, ['47']),
             (5, ['3185', '10373', '841', '10663']),
             (6, ['967']),
             (7, ['1468']),
             (8, ['11649', '1510']),
             (9, ['1974', '1900', '1880', '3756']),
             (10, ['1934']),
             (11, ['1956']),
             (12, ['2459']),
             (13, ['3035']),
             (14, ['3236']),
             (15, ['3770']),
             (16, ['9929', '9946', '5747']),
             (17, ['5883', '5908']),
             (18, ['6158']),
             (19, ['6984', '6617']),
             (20, ['7179']),
             (21, ['7261', '7243']),
             (22, ['7508']),
             (23, ['7852']),
             (24, ['7898']),
             (25, ['8544', '8554']),
             (26, ['9155', '9143']),
             (27, ['9400', '9376']),
             (28, ['9860', '9837']),
             (29, ['11530']),
             (30, ['98']),
             (31, ['724']),
             (32, ['857']),
             (33, ['2648', '4946', '2907', '4919', '1840', '929', '1469']),
             (34, ['966']),
             (35, ['1102']),
             (36, ['1417']),
             (37, ['2030']),
             (38, ['2137']),
             (39, ['3403']),
             (40, ['3609', '4112']),
             (41, ['4904', '4223']),
             (42, ['6994', '4434', '7529', '7839']),
             (43, ['6130']),
             (44, ['6162']),
             (45, ['6374']),
             (46, ['8630']),
             (47, ['9982', '9947', '9962']),
             (48,
              ['39',
               '338',
               '464',
               '706',
               '45',
               '316',
               '50',
               '120',
               '473',
               '180',
               '64',
               '324',
               '67',
               '173']),
             (49, ['83', '1822']),
             (50, ['128', '129', '138', '1279']),
             (51, ['234']),
             (52, ['405']),
             (53, ['712']),
             (54, ['827']),
             (55, ['1123']),
             (56, ['1263']),
             (57, ['1285']),
             (58, ['1904', '1882'])])

grouped_features_nm = collections.OrderedDict([(1, ['165.078']),
             (2, ['103.047']),
             (3, []),
             (4, []),
             (5, ['305.155']),
             (6, []),
             (7, []),
             (8, ['747.189']),
             (9, ['323.068']),
             (10, []),
             (11, []),
             (12, []),
             (13, []),
             (14, []),
             (15, []),
             (16, ['566.129']),
             (17, []),
             (18, []),
             (19, ['443.156']),
             (20, []),
             (21, []),
             (22, []),
             (23, []),
             (24, []),
             (25, []),
             (26, []),
             (27, []),
             (28, []),
             (29, []),
             (30, []),
             (31, []),
             (32, []),
             (33, ['343.124']),
             (34, []),
             (35, []),
             (36, []),
             (37, []),
             (38, []),
             (39, []),
             (40, ['334.148']),
             (41, ['321.025', '338.052']),
             (42, ['443.251']),
             (43, []),
             (44, []),
             (45, []),
             (46, []),
             (47, []),
             (48, ['181.073']),
             (49, ['274.209']),
             (50, ['119.042']),
             (51, []),
             (52, []),
             (53, []),
             (54, []),
             (55, []),
             (56, ['317.077']),
             (57, []),
             (58, ['246.039', '343.107'])])


score_rules = {'[3M+2Na+2Li-H]3+': '9', '[2M-2H+Al3+]+': '4', '[3M+Li]+': '4', '[2M+3Na-H]2+': '7', '[2M+2H]2+': '3', '[M+3H-COCH2]3+': '5', '[M+H-C4H6O]+': '3', '[M+H-C3H9N-C2H4O2]+': '3', '[2M+Na+2K+Li-H]3+': '8', '[M+2H-NH3]2+': '4', '[2M+2K-H]+': '4', '[3M+H]+': '3', '[2M+K+Li-H]+': '6', '[M+H-CH3OH]+': '3', '[M+H-NH3-CO-CO]+': '3', '[M+Na+3K-H]3+': '7', '[2M+Na+K+2Li-H]3+': '8', '[M+2K]2+': '3', '[M+H-C6H10O5]+': '3', '[3M+2K]2+': '5', '[2M+H+Na]2+': '3', '[2M+Na+K+Li-H]2+': '7', '[M+K+2Li]3+': '4', '[2M+3H]3+': '4', '[3M+K]+': '4', '[M+H-C3H9N]+': '3', '[4M-2H+Fe3+]+': '6', '[M+Na+Li]2+': '3', '[M+H+Li]2+': '2', '[M+H-COCH2-C4H8]+': '3', '[M+3Na+Li-H]3+': '7', '[M+Na+2K]3+': '4', '[M+2Na+K-H]2+': '6', '[3M+2Li]2+': '5', '[2M+K]+': '3', '[M+H-HCOOH]+': '3', '[3M+Na+2K-H]2+': '8', '[M+Li+HCOOH]+': '3', '[M+3H-C6H10O4]3+': '5', '[4M+2K-H]+': '6', '[2M+Na+K]2+': '4', '[3M+Na]+': '4', '[M+H-phenylacetyl]+': '3', '[M+H+K]2+': '2', '[M+Na+K+2Li-H]3+': '7', '[M+3H-HCOOH]3+': '5', '[M+Na+3Li-H]3+': '7', '[M+K+3Li-H]3+': '7', '[M+H-H2O-C2H2O2]+': '3', '[4M+2Na]2+': '6', '[M+2H-CO]2+': '4', '[4M+H]+': '4', '[3M+Na+2Li]3+': '6', '[M+H-C8H6O-H2O]+': '3', '[2M+2Na+2Li-H]3+': '8', '[2M+2Na-H]+': '4', '[M+H-O]+': '3', '[M+H-Leu]+': '3', '[M+2H-HCOOH]2+': '4', '[3M+Na+2Li-H]2+': '8', '[M+H-C3H4O]+': '3', '[M+2Na-H]+': '3', '[3M+NH4]+': '3', '[3M+K+Li]2+': '5', '[M+H-Val]+': '3', '[M+H-NH3-CO2-CH2O]+': '3', '[M+H-SO3-H2O]+': '3', '[M+H-S-NH3-HCOOH]+': '3', '[2M+2K]2+': '4', '[2M+2Na+K]3+': '5', '[M+H-C4H8-C4H6]+': '3', '[M+H-CH3]+': '3', '[3M+H+K]2+': '4', '[3M+Na+2K+Li-H]3+': '9', '[M+3H-C4H8]3+': '5', '[2M+3K-H]2+': '7', '[3M+K+2Li]3+': '6', '[M+NH4+HCOOH]+': '2', '[3M+Na+K+2Li-H]3+': '9', '[M+H-NH3-C3H4]+': '3', '[M+H+(CH3)2CO-H2O]+': '4', '[2M+K+2Li]3+': '5', '[3M+2Na+Li-H]2+': '8', '[M+H-H2O-CO2]+': '3', '[M+H-C4H6]+': '3', '[M+2H-C3H2O3]2+': '4', '[2M+2K+2Li-H]3+': '8', '[M+2Li]2+': '3', '[M+H+KCl]+': '2', '[M+Na+2Li]3+': '4', '[M+H-NH3-CO2]+': '3', '[2M+2K+Li]3+': '5', '[M+2H-CH3]2+': '4', '[2M+3Li-H]2+': '7', '[2M+Na+Li]2+': '4', '[3M+3K-H]2+': '8', '[M+Na+2K-H]2+': '6', '[M+H-NH3-NH3-C3H4]+': '3', '[3M+2Na+K]3+': '6', '[M+H-NH3-CO2-C3H4O]+': '3', '[4M+H+Na]2+': '5', '[M+H-CH4]+': '3', '[M+H-C8H6O]+': '3', '[M+3Na+K-H]3+': '7', '[M+2H-CO2]2+': '4', '[M+H-NH3-C3H4-COCH2]+': '3', '[3M+3Na+K-H]3+': '9', '[2M+H]+': '2', '[M+H-C2H4O2]+': '3', '[M+3K+Li-H]3+': '7', '[M+H-H2O-H2O-C2H4O (McLafferty)]+': '3', '[M+H-(H2O)3-CO]+': '3', '[M+NH4]+': '1', '[M+H-NH3-C2H6]+': '3', '[M+Li]+': '2', '[3M+2Na+K+Li-H]3+': '9', '[M+H-C3H6]+': '3', '[M+H-C2H4]+': '3', '[M+H-CO]+': '3', '[2M-2H+Fe3+]+': '4', '[M+H-CO2-C3H6]+': '3', '[M+H-H2O-HCOOH]+': '3', '[M+2H-C6H8O6]2+': '4', '[3M+2K-H]+': '5', '[M+K]+': '2', '[M+H+(NaCl)2]+': '3', '[M+2H-CH2]2+': '4', '[M+H-NH3-COCH2]+': '3', '[M+H-gluc+H2O]+': '3', '[M+H-S]+': '3', '[2M+3Na+K-H]3+': '8', '[M+H-gluc]+': '3', '[4M+Na]+': '5', '[2M+K+3Li-H]3+': '8', '[3M+K+2Li-H]2+': '8', '[M+H-C3H6O]+': '3', '[M+Li+NaCOOH]+': '3', '[3M+Na+Li]2+': '5', '[M+Na+HCOOH]+': '2', '[M+H-C4H6-H2O]+': '3', '[M+H-C3H2O3]+': '3', '[2M+Na+2Li]3+': '5', '[M+H-C4H6-COCH2]+': '3', '[M+2Li-H]+': '5', '[M+H-NH3]+': '3', '[M+H-C2H4O2-CH3OH]+': '3', '[2M+2Na+Li-H]2+': '7', '[4M+2K]2+': '6', '[4M+3H]3+': '6', '[2M+2Na+K+Li-H]3+': '8', '[4M-2H+Al3+]+': '6', '[M+3H-CH4]3+': '5', '[3M+2K+Li-H]2+': '8', '[M+2Na+Li-H]2+': '6', '[2M+Na+3Li-H]3+': '8', '[3M+2K+Li]3+': '6', '[2M+2Na+K-H]2+': '7', '[M+H-C3H4O-C4H6]+': '3', '[3M+H+Na]2+': '4', '[M+H-gluc-(H2O)3]+': '3', '[2M+Na]+': '3', '[M+H-SO3-H2O-NH3]+': '3', '[M+2H-C2H4]2+': '4', '[M+K+NH3]+': '3', '[M+H-NH3-C8H6O-CH2]+': '3', '[M+K+Li]2+': '3', '[3M+Na+3Li-H]3+': '9', '[M+H+Na]2+': '2', '[M+H-CO2]+': '3', '[M+Na+2Li-H]2+': '6', '[M+3H-H20]3+': '5', '[3M+2Na+K-H]2+': '8', '[3M+K+3Li-H]3+': '9', '[M+3H-C5H8O4]3+': '5', '[M+3H-CO2]3+': '5', '[M+H+HCOOH]+': '2', '[M+H-C4H8]+': '3', '[M+3H-C3H2O3]3+': '5', '[M+Li+NH3]+': '3', '[3M+Na+K-H]+': '7', '[M+2Na+2Li-H]3+': '7', '[M+H-C8H6O-NH3]+': '3', '[4M+K]+': '5', '[M+H-CH3S]+': '3', '[3M+3H]3+': '5', '[M+3Na-H]2+': '6', '[3M+2K+2Li-H]3+': '9', '[M+H-HCOOH-HCOOH]+': '3', '[2M+2Li]2+': '4', '[M+Na+2K+Li-H]3+': '7', '[M+H-C2H4-HCOOH]+': '3', '[M+H-gluc-(H2O)3-CO]+': '3', '[2M+H+Li]2+': '3', '[M+H-H2O]+': '3', '[M+H-NH3-CO2-NH3-H2O]+': '3', '[M+3Li-H]2+': '6', '[M+K+NaCOOH]+': '3', '[3M+Na+K+Li-H]2+': '8', '[2M+Na+2K]3+': '5', '[2M+NH4]+': '3', '[2M+Na+Li-H]+': '6', '[M+2K-H]+': '3', '[M+Na+NH3]+': '3', '[M+H-(H2O)2]+': '3', '[M+H-NH3-CO-COCH2]+': '3', '[3M+H+Li]2+': '5', '[2M+Na+3K-H]3+': '8', '[M+2Na+2K-H]3+': '7', '[3M+3Li-H]2+': '8', '[M+H-C4H8O2]+': '3', '[M+Na+K+Li]3+': '4', '[3M+Na+Li-H]+': '7', '[M+H-C5H10]+': '3', '[2M+H+K]2+': '3', '[M+H+CHOONa]+': '2', '[M+H-C4H6-C4H6O]+': '3', '[M+H-NH3-HCOOH-CH3OH]+': '3', '[M+2H-H20]2+': '4', '[M+H-NH3-H2O-H2O]+': '3', '[M+3H-CO]3+': '5', '[M+H-C6H10O4]+': '3', '[M+2H-C6H10O5]2+': '4', '[M+H-CH2O]+': '3', '[M+2Na+K+Li-H]3+': '7', '[M+2K+2Li-H]3+': '7', '[M+H-NH3-H2O]+': '3', '[3M+Na+K+Li]3+': '6', '[M+H-C3H6O-CH3OH]+': '3', '[3M+2Na+2K-H]3+': '9', '[3M+2Na]2+': '5', '[3M-2H+Al3+]+': '5', '[M+3H-NH3]3+': '5', '[M+2H-C4H8]2+': '4', '[2M+Na+K-H]+': '6', '[3M+Na+3K-H]3+': '9', '[M+H-C4H6-NH3-H2O]+': '3', '[M+H-C2H2]+': '3', '[M+Na+Li-H]+': '5', '[4M+2Na-H]+': '6', '[M+2H-C6H10O4]2+': '4', '[2M+2Li-H]+': '6', '[2M+Na+2K-H]2+': '7', '[M+H-NH3-HCOOH]+': '3', '[M+Na+K]2+': '3', '[3M+2H]2+': '4', '[M+K+HCOOH]+': '2', '[M+2H]2+': '3', '[M+H+NH3]+': '1', '[4M+H+K]2+': '5', '[M+3H-C6H10O5]3+': '5', '[M+2Na+Li]3+': '4', '[M+2Na]2+': '3', '[M+H-C5H8O4]+': '3', '[2M+Na+2Li-H]2+': '7', '[2M+3K+Li-H]3+': '8', '[M+H-COCH2]+': '3', '[2M+2Na+2K-H]3+': '8', '[3M+2Na-H]+': '5', '[M+H-C4H6-C2H4]+': '3', '[M+2H-COCH2]2+': '4', '[M+3H]3+': '3', '[M+3H-C6H8O6]3+': '5', '[M+3K-H]2+': '6', '[M+H-CH2]+': '3', '[3M+Na+K]2+': '5', '[3M+2Na+Li]3+': '6', '[M+K+Li-H]+': '5', '[3M+3K+Li-H]3+': '9', '[M+H-C5H8O]+': '3', '[2M+K+Li]2+': '4', '[M+H-C2H4O]+': '3', '[M+H-C2H2O2]+': '3', '[M+H-gluc-(H2O)2]+': '3', '[M+H-C2H4-CO2]+': '3', '[M+H-C5H8]+': '3', '[2M+H+NH4]2+': '3', '[2M+2Na+Li]3+': '5', '[2M+2Na]2+': '4', '[4M+Na+K]2+': '6', '[M+3H-C2H4]3+': '5', '[M+H-C3H4]+': '3', '[M+Na+K+Li-H]2+': '6', '[2M+K+2Li-H]2+': '7', '[M+2K+Li-H]2+': '6', '[3M+3Na+Li-H]3+': '9', '[M+H-C6H12]+': '3', '[3M+H+NH4]2+': '5', '[3M-2H+Fe3+]+': '5', '[2M+Li]+': '3', '[M+3H-CH2]3+': '5', '[M+3H-CH3]3+': '5', '[M+H+NaCl]+': '2', '[3M+2Li-H]+': '7', '[M+2H-C5H8O4]2+': '4', '[2M+2K+Li-H]2+': '7', '[M+H+NH4]2+': '2', '[M+H-C2H6]+': '3', '[M+H-NH3-CO2-C5H8]+': '3', '[M+H-NH3-CO-COCH2-C4H6O]+': '3', '[M+H-C6H8O6]+': '3', '[3M+Na+2K]3+': '6', '[2M+Na+K+Li]3+': '5', '[M+H-(H2O)3]+': '3', '[M+H]+': '1', '[4M+2H]2+': '5', '[M+H-SO3]+': '3', '[M+H-H20]+': '3', '[M+Na+K-H]+': '3', '[M+Na]+': '2', '[M+2Na+K]3+': '4', '[M+H-C3H4O-C4H8O2]+': '3', '[M+2K+Li]3+': '4', '[M+Na+NaCOOH]+': '3', '[2M+3Na+Li-H]3+': '8', '[M+2H-CH4]2+': '4', '[M+K+2Li-H]2+': '6', '[3M+K+Li-H]+': '7', '[3M+3Na-H]2+': '8'}

bas = [1.0, 1.0, 8.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0,
 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0,
11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0,
 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0]


totalS = [78.921, 55.196, 39.0, 36.662, 34.732, 32.598, 32.306, 32.241, 32.239, 31.517, 30.953, 30.809, 30.765, 30.759, 30.745, 24.817, 15.016, 11.337, 11.091, 10.144, 9.857, 9.851,
          9.813, 9.802, 9.754, 9.751, 9.747, 9.219, 9.204, 9.127, 9.095, 9.091, 9.084, 9.072, 9.069, 9.066, 9.049, 9.039, 9.038, 9.038, 9.03, 9.027, 9.027, 9.026, 9.023, 9.022, 9.017,
          9.014, 9.013, 9.009, 9.008, 9.008, 9.008, 9.008, 9.007, 9.001, 9.0, 9.0, 7.477, 4.233, 3.567, 2.823, 2.028, 1.612, 1.519, 1.073, 0.975, 0.784, 0.715, 0.712, 0.669, 0.587,
          0.416, 0.366, 0.346, 0.342, 0.339, 0.332, 0.314, 0.309, 0.283, 0.263, 0.257, 0.207, 0.202, 0.197, 0.195, 0.194, 0.149, 0.149, 0.117, 0.108, 0.094, 0.092, 0.083, 0.083,
          0.078, 0.074, 0.073, 0.065, 0.062, 0.061, 0.058, 0.053, 0.047, 0.047, 0.04, 0.039, 0.037, 0.036, 0.032, 0.03, 0.029, 0.028, 0.027, 0.027, 0.027, 0.025, 0.025, 0.024,
          0.023, 0.018, 0.015, 0.013, 0.011, 0.011, 0.009, 0.009, 0.009, 0.008, 0.008, 0.008, 0.008, 0.007, 0.005, 0.005, 0.003, 0.003, 0.003, 0.002]


excluded = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1
, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


excluded_final_dims = [0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1]

excluded_final_dims_purity_na = [0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1]