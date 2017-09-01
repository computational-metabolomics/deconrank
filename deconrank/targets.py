from __future__ import (
    print_function
)
import csv
import os
import numpy as np

def create_dims_targets(d_table, max_time=1800, min_time=120, max_cid_time=300, peak_time_hcd=10,
                        peak_time_cid=12, delay_time=0.24, cid_perc=10):
    #==============================================
    # Create the dynamic peak lists for msn and ms
    #==============================================
    max_time = max_time  # 1800 seconds (30 mins) (plus delay time of 24 seconds)
    min_time = min_time  # 120 seconds (2 mins)

    max_cid_time = max_cid_time # (5 mins)

    peak_time_hcd = peak_time_hcd
    peak_time_cid = peak_time_cid

    # all peaks (after excluding percentage of second tier)
    # We only want to use "included" peaks for this
    maskE = d_table['excluded'] < 1
    d_table_small = d_table[maskE]

    total_peak_nm = d_table_small.shape[0]

    # Get the length of time for the cid section
    delay_time = delay_time

    cid_peak_nm = round(total_peak_nm * cid_perc, 0) # rounded to represent 1 peak
    cid_total_time = cid_peak_nm * peak_time_cid

    # we have a max time limit for the CID section
    if cid_total_time > max_cid_time:
        cid_total_time = max_cid_time

    # Get the max time for the HCD section
    max_time_hcd = max_time - cid_total_time

    # Get the maximum number of peaks that can be assessed for HCD within the HCD time limit
    peak_max_limit = max_time_hcd / peak_time_hcd

    # If the total number of peaks is > tham the max that can be fragmented, cut down the peaklist
    # else get a peaklist size based on the number of peaks available.

    # Final inclusion/exclusion column based on time constraints
    # NOTE: now already done at the d_table creation stage

    pids = d_table['peakID']

    if total_peak_nm > peak_max_limit:
        hcd_total_time = max_time_hcd+delay_time

        # only get the top peaks
        npa_small = d_table_small[:peak_max_limit, ]

        # Final flag for the full peak list
        mask = np.in1d(pids, d_table_small['peakID'])
        d_table['excludedFinal'][mask==False] = 1

    else:
        hcd_total_time = (total_peak_nm * peak_time_hcd)+delay_time

    end_time = hcd_total_time+cid_total_time

    # Ensure that a min time is always used
    if end_time < min_time:
        end_time = min_time

    hcd_total_time_min = str(round(hcd_total_time/60.0, 2))
    end_time_min = str(round(end_time/60.0, 2))

    # Create the required columns in suitable format for the target list file
    mz = d_table_small['mz']
    mz = mz.astype(float)
    mz = np.round(mz, 5)
    space = np.repeat("\t", len(mz))
    s = np.repeat(0, len(mz))
    f = np.repeat(end_time_min, len(mz))
    c = np.core.defchararray.add('peak_', d_table_small['peakID'])

    # if no peaks then return
    if d_table_small.shape[0]==0:
        return 'NO_PEAKS'

    targets = np.column_stack((mz, s, f, space, space, c))

    return targets, end_time_min, hcd_total_time_min



def write_out_dims_targets(suffix, out_dir, end_time_min, hcd_total_time_min, targets, pol):
    #==============================================
    # Write out
    #==============================================
    nm_target = suffix + "_target.txt"

    # saveand target list
    tdir = os.path.join(out_dir, "targets")
    if not os.path.exists(tdir):
        os.makedirs(tdir)
    nt_target = os.path.join(tdir, nm_target)

    xcalibur_auto_pth = os.path.join(out_dir, "XcaliburAutoInput.txt")

    # Get file pats and methods for XcaliburAutoInput file
    if pol=="POS":
        meth_template = "C:\\Documents and Settings\\XPMUser\\Desktop\\XcaliburAuto\\templates\DMA_Pos_MS_MSMS_MSn_Final_v3.meth"
    else:
        meth_template = "C:\\Documents and Settings\\XPMUser\\Desktop\\XcaliburAuto\\templates\DMA_nESI_Neg_Polar_FINAL_MSMS_MSn_3e6_TEMPLATE.meth"

    target_ex = "C:\\Documents and Settings\\XPMUser\Desktop\\XcaliburAuto\\targets\\"+nm_target
    method_ex = "C:\\Documents and Settings\\XPMUser\\Desktop\\XcaliburAuto\\methods\\"+suffix+"n.meth"

    # xcalibur auto input string
    xcalibur_auto_txt = [meth_template, target_ex, end_time_min, hcd_total_time_min, method_ex]

    # add to string to bottom of file (if file exists, if nt the create file and add to first line)
    try:
        with open(xcalibur_auto_pth, "a") as file:
            file.write(' '.join(xcalibur_auto_txt))
    except IOError:
        with open(xcalibur_auto_pth, "a") as file:
            file.write(' '.join(xcalibur_auto_txt))

    # Write out target
    np.savetxt(nt_target, targets, fmt="%s", delimiter="\t")

    return xcalibur_auto_txt