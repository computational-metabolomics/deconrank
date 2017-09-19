from __future__ import (
    print_function
)
import csv
import re
import os
import numpy as np



def initial_filter(d_table, grouped_feature_d, stier_perc=0.5, pthr=None, irm=None ):

    # order based on total score (descending) i.e. Best score at the top 100, 50, 20 etc .....
    d_table = d_table[d_table['totalS'].argsort()[::-1]]

    #==============================================
    # Remove duplicate mz peaks
    #==============================================
    # remove mz values for inclusion list that are within 5 decimal places of one another

    mzstrt = np.round(d_table['mz'], 5)

    uni, unidx = np.unique(mzstrt, return_index=True)

    pre = d_table.shape[0]
    d_table = d_table[unidx]
    post = d_table.shape[0]


    if pre is not post:
        print("Removed duplicates from inclusion list, total of: ", post-pre)


    # Have to reorder
    d_table = d_table[d_table['totalS'].argsort()[::-1]]



    ###############################################
    # Dynamic peak list creation
    ###############################################
    # add a column to identify if peak is to be exclded based on global peak properties (pthr, snr, irm etc)
    # NOTE: this is now already created when we generate the d_table


    #==============================================
    # Only keep peaks greater than n purity
    #==============================================
    if pthr:
        d_table['excluded'][d_table['medianPurity'] <= pthr] = 1

    #==============================================
    # Only keep peaks greater than n SNR
    #==============================================
    # d_table[d_table[:, cn['snr']].astype(np.float) <= snr, cn['excluded']] = 1

    #==============================================
    # Remove isotopes
    #==============================================
    ic = 0
    for i in d_table['isotopes']:
        for ir in irm:
            if re.match(".*"+re.escape(ir)+".*", i):
                d_table['excluded'][ic] = 1
                break
        ic += 1




    #==============================================
    #  Remove n percentage of the second tier peaks for each peak cluster
    #==============================================
    c = 0
    pids = d_table['peakID']
    for peaks in d_table:

        # Get peaks from the peak cluster

        # skipped second tier peaks
        if peaks['groupid']=='NA':
            continue

        cids = grouped_feature_d[int(peaks['groupid'])]

        # Only do this for peak clusters
        if len(cids)<2:
            continue

        # get peaks in the cluster
        mask = np.in1d(pids, cids)

        peakclust = d_table[mask]

        # Get only the second tier peaks
        peakclust = peakclust[peakclust['groupid'] == 'NA']

        # Only include included peaks that have not already been excluded
        peakclust = peakclust[peakclust['excluded']<1]

        # Order by descending score
        peakclust = peakclust[peakclust['totalS'].argsort()[::-1]]

        # Flag the last 50% peaks
        stier_lim = round(peakclust.shape[0] * stier_perc, 0)
        peakclust = peakclust[int(stier_lim):peakclust.shape[0]+1]
        mask = np.in1d(pids, peakclust['peakID'])

        d_table['excluded'][mask] = 1

    # excludedFinal is used later on but we initiate here to avoid confusion later (hopefully!)
    d_table['excludedFinal'] = d_table['excluded']


    # with open('test_scores.csv', 'wb') as csvfile:
    #     w = csv.writer(csvfile, delimiter=',')
    #     w.writerow(d_table.dtype.names)
    #     # print(d_table.dtype.names)
    #     for d in d_table:
    #         w.writerow(d)

    return d_table
