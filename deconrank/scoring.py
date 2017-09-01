from __future__ import (
    print_function
)
import csv
import re
import os
import numpy as np


def read_in_rules(fileIn, ruledir=None, pol=None, rule_pth=None):
    # Checks the polarity of the file, then turns the scored adduct rule file into a dictionary

    print("###getting rules for adduct scoring (polarity:",)

    if not rule_pth:
        if ruledir:
            rule_dir = ruledir
        else:
            rule_dir = os.path.dirname(os.path.realpath(__file__))

        if re.match(".*pos.*", fileIn, re.IGNORECASE) or pol == "pos":
            polarity = "POS"
            print("positive)###")
            # rules_pth = os.path.join(rule_pth, "CAMERA_rules_pos.csv")
            rule_pth = os.path.join(rule_dir, "CAMERA_rules_PosFinal_PlusLi.csv")

        elif re.match(".*neg.*", fileIn, re.IGNORECASE) or pol == "neg":
            print("negative)###")
            polarity = "NEG"
            # rules_pth = os.path.join(rule_pth, "CAMERA_rules_neg.csv")
            rule_pth = os.path.join(rule_dir, "CAMERA_rules_NegFinal_PlusLi.csv")
        else:
            print("Can not determine if file input is positive or negative, the file name should include "
                  "either the text 'pos', or 'neg' ")
            quit()

    # Turn rules file into simple dictionary of adduct name and score
    try:
        with open(rule_pth, 'rb') as csvfile:

            r = csv.reader(csvfile, delimiter=',', quotechar='"')

            h = next(r, None)

            f_score = h.index('frag_score')

            name = h.index('name')

            scores = {}
            for row in r:
                scores[row[name]] = row[f_score]

    except IOError as e:
        print("Unable to open rule file")
        quit()

    return scores, polarity

def get_mz_idx(header):
    if 'mz' in header:
        return header.index('mz') - 1
    elif 'mzmed' in header:
        return header.index('mzmed') - 1
    else:
        print("mz column not found, need either [mz, mzmed] in columns")
        quit()

def score_adducts(grouped_features, grouped_features_nm, features, header, scores):
    # Run through the grouped features dictionary and assign adduct and cluster size scores
    print("\tscoring")
    full_list = []

    # Initialise the output file.
    print("\tgrouping features")

    missing_scored = []

    try:
        adduct_idx = header.index('adduct') - 1
    except ValueError as e:
        print("Adduct not in column names, error :{}".format(e))
        quit()

    if 'i' in header:
        i_idx = header.index('i') - 1
    elif 'intensity' in header:
        i_idx = header.index('intensity') - 1
    elif 'sample_median_I' in header:
        i_idx = header.index('sample_median_I') - 1
    else:
        print("Intensity column not found, need either [i, intensity, sample_median_I] in columns")
        quit()

    mz_idx = get_mz_idx(header)

    for g, grp_f in grouped_features.iteritems():
        temp = []

        # Neutral mass of group
        nm = grouped_features_nm[g]

        tracebackmz = []
        tracebacki = []
        tracebackaddf = []

        for idx in grp_f:

            feature_info = features[idx]

            # first get the adduct from the feature dictionary
            adduct_full = feature_info[adduct_idx]
            mz = float(feature_info[mz_idx])
            inten = float(feature_info[i_idx])

            tracebackmz.append(mz)
            tracebacki.append(inten)
            tracebackaddf.append(adduct_full)

            items = adduct_full.split(" ")

            # remove that pesky McLafferty
            items = [x for x in items if x not in ['(McLafferty)]+', '(McLafferty)']]

            # If there is a no adduct description then feature is scored high (8)
            if not adduct_full:
                temp.append(
                    [str(idx), inten, 8])
            else:
                # There might be multiple adduct explanations of the feature. Need to choose the
                # explanation that fits the grouped_feature neutral mass (nm)
                for i in range(0, len(items), 2):

                    # Stores, as a tuple, the adduct form and the the nm value for that adduct.
                    adduct_sml = items[i:i + 2]

                    # Adduct has to match one of the nm of the feature group (can be >1 in some cases)
                    if adduct_sml[1] in nm:

                        # Get the score out of the scores dictionary, if there is either no adduct
                        # annotation or the adduct has not been scored in the ruleset file, the feature
                        # will be scored 0
                        try:
                            score = float(scores[adduct_sml[0]])
                        except (KeyError, ValueError) as e:
                            # if either not scored or not present in ruleset file, or just no adduct score
                            # as the max
                            score = 6
                            # Record any adducts that were not scored in the ruleset files to show user
                            missing_scored.append(adduct_sml[0])

                        temp.append(
                            [str(idx), inten, score])

        # Sort grouped features, first by adduct type (first tier), then by intensity note the minus sign
        # means descending (highest to lowest!!), so x[2] adduct (ascending) and the -x[1] intensity (descending)
        # This means if we get adducts of the same score we use the most intense for the first tier adduct
        temp.sort(key=lambda x: (x[2], -x[1]))

        clustn = len(tracebackmz)

        # tracebackmz = [round(i, 4) for i in tracebackmz]
        # tracebacki = [round(i, 0) for i in tracebacki]
        #
        # tracebackmz = " ".join(map(str, tracebackmz))
        # tracebacki = " ".join(map(str, tracebacki))
        # gf = " ".join(grouped_features[g])

        # weighted length score tmpe[0][4] is the adduct score of the best match
        # wscore = clustn/float(temp3[0][4])

        # Add clustn score and traceback details
        flrow = temp[0] + [clustn, g]

        # Add most intense/best adduct form of feature to the full list file
        full_list.append(flrow)

    storedidx = [row[0] for row in full_list]
    # Add the features to end that were not the picked adducts and give the lowest possible score
    for fid, feature in features.iteritems():
        # Check if not in the full_list of features
        if fid not in storedidx:
            full_list.append([str(fid),
                              feature[i_idx],
                              11,
                              1,
                              'NA'])

    missing_scored = filter(None, missing_scored)

    if set(missing_scored):
        print("The following adduct annotation did not have associated frag_score so "
              "were given the default score of 6:", set(missing_scored))

    return full_list

def create_score_table(scored_adduct_list, features, header):
    # loop through scored_adduct_list
    # extract the feature associated with each score
    # get the values used for scores: medianPurity (but potentially could add more)

    mz_idx = get_mz_idx(header)
    l = []
    for s in scored_adduct_list:
        if 'medianPurity' in header:
            medianPurity = features[s[0]][header.index('medianPurity') - 1]
        else:
            medianPurity = 0

        mz = features[s[0]][mz_idx]

        isotope = features[s[0]][header.index('isotopes') - 1]
        adduct = features[s[0]][header.index('adduct') - 1]

        peakID = s[0]
        intensity = s[1]
        bestAdductScore = s[2]
        clustn = s[3]
        groupid = s[4]

        l.append(tuple(
            [peakID, mz, intensity, isotope, adduct, groupid, bestAdductScore, clustn, medianPurity] + [0] * 11))

    d_table =  load_score_table(l)

    return d_table




def load_score_table(l):
    d_table = np.array(l, dtype= [('peakID', 'S8'),
                                  ('mz', 'f8'),
                                  ('intensity', 'f8'),
                                  ('isotopes', 'S15'),
                                  ('adduct', 'S100'),
                                  ('groupid', 'S4'),
                                  ('bestAdductScore', 'f8'),
                                  ('clustn', 'f8'),
                                  ('medianPurity', 'f8'),
                                  ('intp', 'f8'),
                                  ('addp', 'f8'),
                                  ('clustnp', 'f8'),
                                  ('pp', 'f8'),
                                  ('intensityw', 'f8'),
                                  ('adductw', 'f8'),
                                  ('clustnw', 'f8'),
                                  ('purityw', 'f8'),
                                  ('totalS', 'f8'),
                                  ('excluded', 'i4'),
                                  ('excludedFinal', 'i4')
                                  ])
    return d_table

def min_max_scale(x):
    if np.unique(x).size == 1:
        # if all the same number then output as list of zero
        return np.zeros(len(x))
    else:
        return (x - min(x)) / (max(x) - min(x))


def normalise_scores(d_table, weights):

    #Adduct score
    #Convert adduct score into percentage
    #(e.g. if worst rank is 11. Then rank 11 = 0%,  rank 10 = 9% rank 9 = 18%... rank 1 = 100%, )
    adds = d_table['bestAdductScore']
    adds = adds.astype(np.float)

    m = adds.max()

    addp = np.zeros(adds.size)

    for i in xrange(adds.size):
        # First reverse the rank score so highest score is the best. e.g. so that 11 is best and 1 is worst
        r = abs(adds[i] - (m + 1))

        # get percentage
        addp[i] = ((r - 1) / (m - 1)) * 100

    # clustnm score
    clustn = d_table['clustn'].astype(np.float)
    clustnp = min_max_scale(clustn) * 100

    # intensity score
    ints = d_table['intensity'].astype(np.float)
    intp = min_max_scale(ints) * 100

    # purity score (already bounded between 0 and 1)
    pp = d_table['medianPurity'].astype(np.float) * 100

    purityw = (pp) * weights['purity']
    adductw = (addp) * weights['adduct']
    intensityw = (intp) * weights['intensity']
    clustnw = (clustnp) * weights['clustn']
    scorew = purityw + adductw + intensityw + clustnw

    d_table['pp'] = pp
    d_table['addp'] = addp
    d_table['intp'] = intp
    d_table['clustnp'] = clustnp
    d_table['purityw'] = purityw
    d_table['adductw'] = adductw
    d_table['intensityw'] = intensityw
    d_table['clustnw'] = clustnw
    d_table['totalS'] = scorew

    idxf = d_table['totalS'].astype(np.float)

    idxd = idxf.argsort()[::-1]

    d_table = d_table[idxd]

    return d_table
