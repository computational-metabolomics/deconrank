from __future__ import (
    print_function
)
import csv
import collections

def group_peaks(fn):
    """
    Information regarding the pcgroup is contained within a hybrid dict key for the adduct data
    Information regarding isotope groups is inherently constrained to given pcgroups by numerical tag for M+ and M+1.
    """
    print("###grouping peaks###")

    with open(fn, 'rb') as csvfile:

        # read in the CAMERA-annotate object for deconvolution (.csv format).
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')

        # Get header from the input file: the first row in the file.
        header = next(reader, None)

        # Record column number for column headed as 'isotopes'
        cln_iso = header.index('isotopes')

        # ignore everything after rsd column (generated from the R workflow script not XCMS/CAMERA)
        # cln_msms_count = header.index('msms_count')

        # Split the input data in to 'pcgroups' - clusters of peaks in similar retention windows.
        groups = collections.OrderedDict()

        # Initialise the AdductGroups dictionary.
        AdductGroups = collections.OrderedDict()

        # Initialise the IsotopeGroups dictionary.
        IsotopeGroups = collections.OrderedDict()

        # Initialise the output dictionary, containing features of interest.
        Features = collections.OrderedDict()

        # Loop through the iterator
        for row in reader:


            # Added by T. N. Lawson - ignores the additional information in final columns for picking of features etc.
            # row = row[:cln_msms_count + 1]

            # Row index number as key (peakID), values comprise all remaining columns
            # Features[str(row[0])] = row[1:]
            Features[row[0]] = row[1:]

            # Check to see if the 'pcgroup' is already in the groups dictionary.
            if row[cln_iso + 2] not in groups:

                # 'pcgroup' number taken as groups key. Values: row index# taken as key in values structure, plus storing values for the isotope and adduct columns
                groups[str(row[cln_iso + 2])] = {str(row[0]): [row[cln_iso], row[cln_iso + 1]]}

            # If 'pcgroup' already in dict 'groups', append values for current row.
            else:
                groups[row[cln_iso + 2]][row[0]] = [row[cln_iso], row[cln_iso + 1]]

            # If Adduct column contains an entry:
            if " " in str(row[cln_iso + 1]):

                # Split the adduct information column on the blank space.
                items = row[cln_iso + 1].split(" ")

                # For those weird adducts calls e.g. McLafferty....
                items = [x for x in items if x not in ['(McLafferty)]+', '(McLafferty)', '(McLafferty)]-']]

                # Chunks the adduct row information in to sequential pairs of adduct label and proposed neutral mass.
                for i in range(0, len(items), 2):

                    # Stores, as a tuple, the adduct form and the the m/z value for that adduct.

                    adduct = items[i:i + 2]

                    # k is being constructed as a string containing both the 'pcgroup' and the neutral mass based on adduct info.
                    k = "%s %s" % (row[cln_iso + 2], adduct[1])

                    # Check to see if the string, k, is a key in OrderedDict 'AdductGroups'
                    if k not in AdductGroups:
                        # Hybrid string used as the key, comprising both the 'pcgroup' for the feature and it's neutral mass.
                        # The values are stored with the index value of the row (smart) and with the adduct form as 'value'
                        AdductGroups[k] = {row[0]: [adduct[0], row[1], row[cln_iso]]}
                    else:
                        # Where multiple adducts are proposed, but all with the same neutral mass, add these to the same 'value' field.
                        AdductGroups[k][row[0]] = [adduct[0], row[1], row[cln_iso]]

                # If the Isotope column is not empty:
                if row[cln_iso] != "":

                    # Split based upon string and take the isotope group number.
                    isoNum = row[cln_iso].split("[M", 2)[0]
                    # If group is not in IsotopeGroups dict., add it as key and insert.
                    if str(isoNum) not in IsotopeGroups:

                        # IsotopeGroups key is the isotope grouping number. This contains dictionary with key = index # and value = mz.
                        IsotopeGroups[isoNum] = {row[0]: [row[cln_iso], row[1]]}
                    else:
                        IsotopeGroups[isoNum][row[0]] = [row[cln_iso], row[1]]

            # Check whether there is information in the isotope group column.
            elif row[cln_iso] != "":

                # Split based upon string and take the isotope group number.
                isoNum = row[cln_iso].split("[M", 2)[0]

                # If group is not in IsotopeGroups dict., add it as key and insert
                if str(isoNum) not in str(IsotopeGroups.keys()):
                    IsotopeGroups[isoNum] = {row[0]: [row[cln_iso], row[1]]}

                else:
                    IsotopeGroups[isoNum][row[0]] = [row[cln_iso], row[1]]

    return Features, AdductGroups, IsotopeGroups, header

def combine(features, AdductGroups, IsotopeGroups):

    """For the refinement steps we should be exploiting the 'groups' dict created earlier to ensure we
    combined co-eluting features and not isomeric species that elute later in the run"""

    print("###refining###")

    #List initialisation
    rts = []
    rts_unfil = []

    #Initialise the dictionary 'grouped_features'
    grouped_features = collections.OrderedDict()

    #Initialise the dictionary 'grouped_features'
    grouped_features_nm = collections.OrderedDict()

    #Begin group count
    group = 1

    #Initialise exclusion list
    excluded = []

    #For key in Features dictionary (i.e. peakID #)
    for peakID in features:

        #Append the fourth value in the dict[key] values object
        rts_unfil.append(float(features[peakID][3]))

        #If the peakID has not already been excluded due to prior grouping
        if str(peakID) not in excluded:

            #Initialise list for related peakID values
            temp = []
            nm = []

            #For key in AdductGroups dict (combination of index # and neutral mass).
            for k, adduct_val in AdductGroups.iteritems():
                if peakID in adduct_val.keys():
                    #Get keys within the dictionary
                    #Structure is dict[key] = {key: [values]} Gets the second key here.
                    #References the unique key index for the metabolic feature.
                    temp.extend(list(AdductGroups[k].keys()))
                    # Get the associated nm for this groupe (could be multiple for cases where there are
                    # alternative hypoth), e.g. an AdductGroup of:
                    #  121.10119 [M+H]+ 120.096 [2M+3H]3+ 180.149
                    #  181.15840 [3M+2H]2+ 120.096 [M+H]+ 180.149
                    # There is v. occasionally bug with situations where the group is a mixture of 2 different
                    # types of adduct forms. (not a massive problem but should be aware)
                    nm.append(k.split(' ')[1])

            #Add peakID to temp list
            temp.append(peakID)

            #Remove redundancy in the group list
            temp2 = set(temp)

            #For all related adducts (with same neutral mass and retention time range - termed pcgroup i.e. peak correlation group), find all associated isotopes for each adduct form.
            for i in IsotopeGroups:

                #Check if peakID from Features dict. matches the key in IsotopeGroups dict.

                for k in IsotopeGroups[i].keys():

                    #for every value in the group list
                    for t in temp2:

                        #check if peakID value is also present in the IsotopeGroups dictionary
                        if str(t) == str(k):

                        #If yes, extend the temporary list to group adducts and isotopes of related features, together. We are adding the unique row index #.
                            temp.extend(list(IsotopeGroups[i].keys()))
                        else:
                            continue

            #Remove redundancy in the group list
            temp2 = set(temp)

            #Add group list to grouped_features dictionary

            if len(grouped_features.keys()) == 0 and len(temp2) >= 1:
                grouped_features[group] = list(set(temp))
                grouped_features_nm[group] = nm
                group += 1

            elif len(grouped_features.keys()) >= 1 and len(temp2) == 1:
                grouped_features[group] = list(set(temp))
                grouped_features_nm[group] = nm
                group += 1

            elif len(temp2) > 1:
                grouped_features[group] = list(temp2)
                grouped_features_nm[group] = nm
                group += 1

            #Add all grouped features to the exclusion list in order to prevent false splitting of features
            excluded.extend(temp2)

    return grouped_features, grouped_features_nm
