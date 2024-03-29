# coding: utf-8
from __future__ import (
    print_function
)
import os
import csv
import argparse
import textwrap

from grouping import group_peaks, combine
from scoring import read_in_rules, score_adducts, create_score_table, normalise_scores
from filter import initial_filter
from targets import write_out_dims_targets, create_dims_targets

class Deconrank(object):

    def __init__(self, in_file, rule_pth=None, polarity=None, out_dir='.', delim=','):
        self.features = None
        self.adduct_groups = None
        self.isotope_groups = None
        self.header =  None
        self.grouped_features = None
        self.grouped_features_nm = None
        self.scores = None
        self.polarity = None
        self.scored_list = None
        self.delim = delim
        if delim == ',':
            self.file_ending = 'csv'
        else:
            self.file_ending = 'tsv'

        self.load_score_rules(in_file, rule_pth, polarity)
        self.suffix = os.path.basename(in_file).split('.')[0]

        self.wd = {'adduct': 0.3, 'intensity': 0.3, 'purity': 0.2, 'clustn': 0.2}

        self.in_file = in_file
        self.out_dir = out_dir

    def load_score_rules(self, in_file, rule_pth, polarity):
        # Checks the polarity of the input file, then turns the scored adduct rule file into a dictionary
        self.scores, self.polarity = read_in_rules(fileIn=in_file, rule_pth=rule_pth, pol=polarity)

    def group(self):
        features, adduct_groups, isotope_groups, header = group_peaks(self.in_file, self.delim)
        grouped_features, grouped_features_nm = combine(features, adduct_groups, isotope_groups)

        self.features =  features
        self.adduct_groups = adduct_groups
        self.isotope_groups = isotope_groups
        self.header = header
        self.grouped_features = grouped_features
        self.grouped_features_nm = grouped_features_nm

    def score(self, weights=None):
        if weights:
            self.wd = {'adduct': weights[0], 'intensity': weights[1], 'purity': weights[2], 'clustn': weights[3]}
        else:
            self.wd = {'adduct': 0.3, 'intensity': 0.3, 'purity': 0.2, 'clustn': 0.2}

        self.scored_adduct_list = score_adducts(self.grouped_features, self.grouped_features_nm,
                                         self.features, self.header, self.scores)

        self.d_table = create_score_table(self.scored_adduct_list, self.features, self.header)

        self.d_table = normalise_scores(self.d_table, self.wd)


    def filter(self, irm=None, stp=0, pthr=None):
        if irm:
            self.irm = irm
        else:
            self.irm = ['[M+1]+','[M+2]+','[M+1]2+','[M+2]2+','[M+1]-','[M+2]-','[M+1]2-','[M+2]2-']

        self.stp = stp
        self.pthr = pthr

        self.d_table = initial_filter(self.d_table, self.grouped_features, irm=self.irm,
                                      stier_perc=self.stp, pthr=self.pthr)


    def dims_targets(self, max_time=1800, min_time=120, max_cid_time=300, peak_time_hcd=10, peak_time_cid=12,
                     delay_time=0.24, cid_perc=10, target_name=None, method_template_name=None):
        print('create dims targets')
        print(self.d_table.shape[0])
        if self.d_table.shape[0]==0:
            return 0
        

        targets, end_time_min, hcd_total_time_min, d_table = create_dims_targets(self.d_table, max_time=max_time, min_time=min_time,
                                                                       max_cid_time=max_cid_time, peak_time_hcd=peak_time_hcd,
                                                                       peak_time_cid=peak_time_cid, delay_time=delay_time,
                                                                       cid_perc=cid_perc)
        self.d_table = d_table
        print('write out targets')
        write_out_dims_targets(out_dir=self.out_dir,
                               end_time_min=end_time_min, hcd_total_time_min=hcd_total_time_min,
                               targets=targets, pol=self.polarity, target_name=target_name,
                               method_template_name=method_template_name)

        self.targets = targets
        self.max_time = max_time
        self.min_time = min_time
        self.max_cid_time=max_cid_time
        self.peak_time_hcd=peak_time_hcd
        self.peak_time_cid=peak_time_cid
        self.delay_time=delay_time
        self.cid_perc=cid_perc

    def write_out_scores(self, modify_original=False, full=False):
        if modify_original:
            outname = os.path.join(self.out_dir, self.suffix+'_scores.'+self.file_ending)
        else:
            outname = os.path.join(self.out_dir, 'scores.' + self.file_ending)

        with open(outname, 'wb') as csvfile:
            w = csv.writer(csvfile, delimiter=self.delim)
            d_table_names = list(self.d_table.dtype.names)
            if full:
                colnames = [d_table_names[0]] + self.header[1:len(self.header)] + d_table_names[6:len(d_table_names)]
            else:
                colnames = d_table_names
            w.writerow(colnames)

            for d in self.d_table:
                if full:
                    d = list(d)
                    w.writerow([d[0]]+self.features[d[0]]+d[6:len(d)])
                else:
                    w.writerow(d)


    def write_out_traceback(self, modify_original=False):
        # Need to add adducts etc
        # print(self.header)

        if modify_original:
            outname = os.path.join(self.out_dir, self.suffix+'_traceback.'+self.file_ending)
        else:
            outname = os.path.join(self.out_dir, 'traceback.' + self.file_ending)

        with open(outname, 'wb') as csvfile:
            w = csv.writer(csvfile, delimiter=self.delim)

            w.writerow(['full_group_idd']+self.header[1:len(self.header)])
            for d in self.d_table:

                if not d['groupid'] == 'NA':
                    gfids = self.grouped_features[int(d['groupid'])]
                    # print(gfids)
                    for gfid in gfids:
                        f = self.features[gfid]
                        # print(f)
                        w.writerow([d['groupid']] + f)


def main():
    import datetime
    import numpy as np
    p = argparse.ArgumentParser(prog='PROG',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''Perform deconvolution and precursor ranking''',
                                epilog=textwrap.dedent('''
                            -------------------------------------------------------------------------

                            Example Usage

                            python deconrank.py -i [in dir] -o [out dir] -p [polarity OPTIONAL] -w [list of weights: adduct,intensity,purity,clustn ]

                            python deconrank.py -i /path/2/camera_out.csv -o /path/2/out_dir/ -w 0.3,0.3,0.2,0.2

                            '''))

    p.add_argument('-i', dest='camera_peaklist_pth', help='path to the camera output file', required=True)
    p.add_argument('-o', dest='out_dir', help='out folder', required=True)
    p.add_argument('--pol', dest='pol', help='polarity [pos, neg], will assign automatically based on'
                                            'input file name if not defined', required=False)
    p.add_argument('--tech', dest='tech', help='Technology used [dims, lcms] default dims', default='dims', required=False)
    p.add_argument('--rp', dest='rp', help='rule path', required=False)
    p.add_argument('--w', dest='w', help='Weights for scoring variables (adduct, intensity, purity, clustern)', required=False)
    p.add_argument('--pthr', dest='pthr', help='Purity threshold', required=False)
    p.add_argument('--stp', dest='stp', help=textwrap.dedent('''
                        Second tier percentage, percent of second tier that go into the target list.
                        Should be in decimal format i.e. 0.1 = 10 percent. Default is to remove all second tier peaks
                        (i.e. set to 0.0)
                        '''), required=False, default=0.0)
    p.add_argument('--irm', dest='irm', help='Isotopes to remove from fragmentation', required=False, default='[M+1]+,[M+2]+,[M+3]+')

    p.add_argument('--max_time', dest='max_time', required=False, default=1800)
    p.add_argument('--min_time', dest='min_time', required=False, default=120)
    p.add_argument('--max_cid_time', dest='max_cid_time', required=False, default=300)
    p.add_argument('--peak_time_hcd', dest='peak_time_hcd', required=False, default=10)
    p.add_argument('--peak_time_cid', dest='peak_time_cid', required=False, default=12)
    p.add_argument('--percentage_cid', dest='percentage_cid', required=False, default=0.3333)
    p.add_argument('--delay_time', dest='delay_time', required=False, default=24)
    p.add_argument('--delim', dest='delim', required=False, default='comma')
    p.add_argument('--modify_name', dest='modify_name', action='store_true')
    p.add_argument('--target_name', dest='target_name', required=False)
    p.add_argument('--method_template_name', dest='method_template_name', required=False)
    p.add_argument('--full_output', dest='full_output', action='store_true')

    st = datetime.datetime.now()
    print("###start time:", st.strftime("%A%d%B%Y_%I%M"), "###")
    args = p.parse_args()

    if not os.path.exists(args.camera_peaklist_pth) or os.path.getsize(args.camera_peaklist_pth)==0:
        print('NO INPUT FILE EXISTS!', args.camera_peaklist_pth)
        return 0

    if args.pol:
        polarity = args.pol
    else:
        polarity = ""

    if args.irm:
        irm_s = args.irm.replace('__ob__', '[')
        irm_s = irm_s.replace('__cb__', ']')
        irm = irm_s.split(',')
    else:
        irm = None
    print(irm)
    
    if args.w:
        weights = args.w.split(',')
        weights = [float(i) for i in weights]
        if sum(weights) > 1:
            print("Weights need to add up to 1")
            return
    else:
        weights = None

    if args.pthr:
        pthr = np.float(args.pthr)
    else:
        pthr = 0

    if args.delim=='comma':
        delim = ','
    elif args.delim=='tab':
        delim = '\t'
    else:
        print("delim needs to be either 'comma' or 'tab'")




    stp = np.float(args.stp)
    dr = Deconrank(in_file=args.camera_peaklist_pth, out_dir =args.out_dir, polarity=polarity, rule_pth=args.rp,
                   delim=delim)
    dr.group()
    dr.score(weights=weights)
    dr.filter(irm=irm, stp=stp, pthr=pthr)

    if not args.tech == 'dims':
        dr.write_out_scores(args.modify_name, full=args.full_output)

    dr.write_out_traceback(args.modify_name)


    if args.tech=='dims':

        if args.target_name:
            target_name = args.target_name
        else:
            target_name = ''

        dr.dims_targets(max_time=float(args.max_time),
                        min_time=float(args.min_time),
                        max_cid_time=float(args.max_cid_time),
                        peak_time_hcd=float(args.peak_time_hcd),
                        peak_time_cid=float(args.peak_time_cid),
                        delay_time=float(args.delay_time),
                        cid_perc=float(args.percentage_cid),
                        target_name= target_name,
                        method_template_name=args.method_template_name)

        dr.write_out_scores(args.modify_name, full=args.full_output)

    ft = datetime.datetime.now()
    d = ft - st
    print("###end time:", ft.strftime("%A, %d. %B %Y %I:%M%p"), "###")
    print("###Time took (min, sec):", divmod(d.days * 86400 + d.seconds, 60), "###")




if __name__ == '__main__':
    main()
    # dr = Deconrank('../tests/data/peaklist_positive.csv', out_dir='../tests/data/')
    # dr = Deconrank('../tests/data/A01_Polar_Daph_WAX1_Phenyl_LCMS_Neg_DIMS_annotated.csv', out_dir='../tests/data/')
    # dr.group()
    # dr.score()
    # dr.filter()
    # dr.write_out_scores()
    # dr.write_out_traceback()
    # dr.dims_targets()
    # dr.lcms_targets()
