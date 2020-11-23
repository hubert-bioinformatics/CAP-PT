#!/bin/python
from __future__ import print_function
import argparse
import os
import openpyxl
import sys
import pandas as pd
from collections import OrderedDict


def argument_parse():
    """ """
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-v',
            '--input-vcf',
            dest='invcf',
            help='A snpeff annotated vcf',
            required=True
            )
    parser.add_argument(
            '-c',
            '--cap-excel',
            dest='c_exc',
            help='A CAP-PT excel file',
            required=True
            )
    parser.add_argument(
            '-r',
            '--recal-bam',
            dest='recal',
            help='A recal bam for making a screenshot of IGV',
            required=True
            )
    parser.add_argument(
            '-b',
            '--bed',
            dest='bed',
            help='A bed region of calling variant',
            required=True
            )    
    args = parser.parse_args()
    invcf = args.invcf
    c_exc = args.c_exc
    recal = args.recal
    bed = args.bed
    return invcf, c_exc, recal, bed


def c_exc_parse(c_exc):
    """
    Parsing CAP-PT entry excel file.
    """
    df_exc = pd.read_excel(c_exc, header=1, encoding='utf-8')
    return df_exc


def bed_parse(bed):
    """ """
    bed_dict = {}
    with open(bed, 'r') as rd_bed:
        for line in rd_bed:
            if not line: break
            if line.startswith('#'): # header line
                continue
            chrom = ''
            start = ''
            end = ''
            S = line.split('\t')
            chrom = S[0].replace('chr', '')
            start = S[1]
            end = S[2]
            if chrom in bed_dict:
                bed_dict[chrom].append((start, end))
            else:
                bed_dict[chrom] = [(start, end)]
    return bed_dict


def invcf_parse(invcf):
    """
    Parsing snpeff annotated vcf file.
    """
    invcf_dict = {}
    vcf_header_list = []
    sample_list = []
    with open(invcf, 'r') as read_vcf:
        for line in read_vcf:
            chrom = ''
            pos = ''
            ref = ''
            alt = ''
            var_key = ''
            info = ''
            form_idx = 0
            form_list = []
            form_val_list = []
            sample_val_dict = {'sample':{}}
            if not line: break
            if line.startswith('##'): continue # vcf description line
            if line.startswith('#CHROM'): # vcf header line
                S = line.rstrip().split('\t')
                vcf_header_list = S
                continue
            S = line.rstrip().split('\t')
            chrom = S[vcf_header_list.index('#CHROM')]
            pos = S[vcf_header_list.index('POS')]
            ref = S[vcf_header_list.index('REF')]
            alt = S[vcf_header_list.index('ALT')]
            var_key = '{}:{}'.format(chrom.replace('chr', ''), pos) # vcf dict key: "chr1:7201362"
            info = S[vcf_header_list.index('INFO')]
            ann_dict = snpeff_parse(info)
            form_idx = vcf_header_list.index('FORMAT')
            form_list = S[form_idx].rstrip().split(':')
            form_val_list = S[form_idx+1:]
            sample_list = vcf_header_list[form_idx+1:]
            for sample in range(0, len(form_val_list)):
                sample_name = ''
                tmp_list = []
                gt = ''
                ad = ''
                dp = ''
                sample_name = vcf_header_list[form_idx+sample+1]
                tmp_list = form_val_list[sample].rstrip().split(':')
                for key, value in zip(form_list, tmp_list):
                    if key == 'GT':
                        gt = value
                        if gt == '0/1': gt = 'HET'
                        elif gt == '0/0' or gt == '1/1': gt = 'HOM'
                        elif gt == '.': gt = '.'
                        else: pass 
                            #raise Exception('Wrong genotype!') # in a case of 1/2, 0/3, etc
                    if key == 'AD': ad = value
                    if key == 'DP': dp = value
                sample_val_dict['sample'][sample_name] = {'GT':gt, 'AD':ad, 'DP':dp}
            invcf_dict[var_key] = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt}
            invcf_dict[var_key].update(ann_dict)
            invcf_dict[var_key].update(sample_val_dict)
    return invcf_dict, sample_list


def snpeff_parse(info):
    """
    Parsing snpeff annotation information.
    """
    ann_dict = {'snpeff':{}}
    S = info.rstrip().split(';')
    for each_info in S:
        if each_info.startswith('ANN='):
            S_eff = each_info.rstrip().split(',')
    eff_cnt = len(S_eff)
    for each_eff in S_eff:
        eff_var_dict = {}
        SS = each_eff.rstrip().split('|')
        eff_func = ''
        eff_gene = ''
        eff_nm = ''
        eff_ccha = ''
        eff_pcha = ''
        eff_key = ''
        eff_func = SS[1]
        eff_gene = SS[3]
        eff_nm = SS[6]
        eff_ccha = SS[9]
        eff_pcha = SS[10]
        if eff_pcha == '':
            eff_pcha = 'p.?'
        else:
            head = eff_pcha.split('.')[1][:3]
            tail = eff_pcha.split('.')[1][-3:]
            num = eff_pcha.split('.')[1][3:-3]
            if head == tail:
                eff_pcha = 'p.{}{}='.format(head, num)
        eff_key = '{}/{}/{}/{}'.format(eff_func, eff_ccha, eff_pcha, eff_nm)
        eff_var_dict['gene'] = eff_gene
        eff_var_dict['nm'] = eff_nm
        eff_var_dict['cchange'] = eff_ccha
        eff_var_dict['snpeff_pchange'] = eff_pcha
        ann_dict['snpeff'][eff_key] = eff_var_dict
    if eff_cnt == len(ann_dict['snpeff']): return ann_dict
    else: raise Exception('SnpEff parsing is wrong!') # checkpoint: snpeff ann count


def make_answer(df_exc, invcf_dict, bed_dict, sample_list):
    """
    Making final result.
    """

    def distinct_nm(cap_nm, vcf_nm):
        """
        Comparing cap-pt nm and vcf nm,
        then sorting out them by NM-matching into one of three types
        : whole match(2) or partial match(1) or unmatch(0)
        """
        nm_criteria = 0 # 0: unmatch / 1: partial match / 2: whole match
        if cap_nm == vcf_nm:
            nm_criteria = 2
            return nm_criteria
        else:
            cap_nm_body = cap_nm.split('.')[0]
            vcf_nm_body = vcf_nm.split('.')[0]
            if cap_nm_body == vcf_nm_body:
                nm_criteria = 1
                return nm_criteria
            else:
                return nm_criteria

    def make_snpeff_final_list(semi_final_list, sample_list, key, invcf_dict, cap_nm):
        """
        Parsing annotation value,
        then assigning it to one of two lists along with nm criteria
            whole match(2): 'cDNA change, protein change' to semi_final_list
            partial match(1): 'cDNA change, protein change, vcf_nm' to semi_final_list
            unmatch(0): 'cDNA change, protein change, vcf_nm' to unmatch_list
        """
        unmatch_list = []
        for ann_key in invcf_dict[key]['snpeff']:
            nm_criteria = 0 # unmatch
            vcf_nm = ann_key.split('/')[-1]
            nm_criteria = distinct_nm(cap_nm, vcf_nm) # comparing cap-pt nm and vcf nm
            if nm_criteria == 2: # whole match
                semi_final_list.append(', '.join(ann_key.split('/')[1:3]))
            elif nm_criteria == 1: # partial match
                semi_final_list.append(', '.join(ann_key.split('/')[1:4]))
            else:
                unmatch_list.append(', '.join(ann_key.split('/')[1:4]))
        semi_final_list.append('|'.join(unmatch_list))
        for sample in sample_list:
            semi_final_list.append(':'.join(invcf_dict[key]['sample'][sample].values()))
        return semi_final_list

    def match_region(cap_chr, cap_start, cap_end, bed_dict):
        """ """
        match_flag = 0
        whole_list = []
        partial_list = []
        cap_rng = range(int(cap_start), int(cap_end)+1)
        cap_len = len(cap_rng)
        for analysis_region in bed_dict[cap_chr]:
            if len(set(cap_rng) & set(range(int(analysis_region[0]), int(analysis_region[1])+1))) == cap_len:
                whole_list.append('{}:{}-{}'.format(cap_chr, analysis_region[0], analysis_region[1]))
            elif len(set(cap_rng) & set(range(int(analysis_region[0]), int(analysis_region[1])+1))) == 0:
                pass
            else:
                partial_list.append('{}:{}-{}'.format(cap_chr, analysis_region[0], analysis_region[1]))
        if len(whole_list) > 0:
            match_flag = 2 # whole match
            return match_flag, whole_list
        elif len(partial_list) > 0:
            match_flag = 1 # partial match
            return match_flag, partial_list
        else:
            return match_flag, []

    final = ''
    final += '\t'.join(list(df_exc) + [
            'match_flag (2: whole, 1:partial, 0:none)',
            'match_region',
            'vcf_Chr',
            'vcf_Pos',
            'vcf_Ref',
            'vcf_Alt',
            'vcf_matched_desc',
            'vcf_unmatched_desc',
            ] + sample_list) + '\n'
    bed_final = ''
    cap_cnt = 1
    for cnt in range(0, len(df_exc)):
        print ('\t{}/{}: {}%'.format(cap_cnt, len(df_exc), round(cap_cnt/1.0/len(df_exc)*100, 2)), end='\r')
        semi_final_list = []
        cap_chr = ''
        cap_start = 0
        cap_end = 0
        cap_nm = ''
        cap_rng = df_exc.loc[:, 'Position/Interval'].iloc[cnt]
        cap_chr = str(cap_rng.split(':')[0])
        if '-' in cap_rng: # interval region
            cap_start = int(cap_rng.split(':')[1].split('-')[0])
            cap_end = int(cap_rng.split(':')[1].split('-')[1])
        else:
            cap_start = cap_rng.split(':')[1]
            cap_end = int(cap_start) + 1
        bed_final += '{}\t{}\t{}\t{}\n'.format(cap_chr, cap_start, cap_end, cnt) # for IGV screenshot
        cap_nm = df_exc.loc[:, 'Transcript '].iloc[cnt]
        match_flag, match_list = match_region(cap_chr, cap_start, cap_end, bed_dict)
        semi_final_list.append(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}'
                .format(
                        df_exc.loc[:, :].iloc[cnt][0],
                        df_exc.loc[:, :].iloc[cnt][1],
                        df_exc.loc[:, :].iloc[cnt][2],
                        df_exc.loc[:, :].iloc[cnt][3],
                        df_exc.loc[:, :].iloc[cnt][4],
                        match_flag,
                        ','.join(match_list).rstrip()
                        )
                )
        if '-' in cap_rng: # interval region
            for var_key in invcf_dict:
                vcf_chr = ''
                vcf_pos = 0
                vcf_chr = str(var_key.split(':')[0])
                vcf_pos = int(var_key.split(':')[1])
                if cap_chr == vcf_chr:
                    if (vcf_pos >= cap_start) and (cap_end >= vcf_pos): # valid variant
                        semi_final_list.append(invcf_dict[var_key]['chrom'])
                        semi_final_list.append(invcf_dict[var_key]['pos'])
                        semi_final_list.append(invcf_dict[var_key]['ref'])
                        semi_final_list.append(invcf_dict[var_key]['alt'])
                        semi_final_list = make_snpeff_final_list(semi_final_list, sample_list, var_key, invcf_dict, cap_nm)
        else: # a single base position
            if cap_rng in invcf_dict:
                semi_final_list.append(invcf_dict[cap_rng]['chrom'])
                semi_final_list.append(invcf_dict[cap_rng]['pos'])
                semi_final_list.append(invcf_dict[cap_rng]['ref'])
                semi_final_list.append(invcf_dict[cap_rng]['alt'])
                semi_final_list = make_snpeff_final_list(semi_final_list, sample_list, cap_rng, invcf_dict, cap_nm)
            else:
                pass
        cap_cnt += 1
        final += '\t'.join(semi_final_list) + '\n'
        print ('\t'.format(), end='\r')
    print('')
    with open('cap-pt.bed', 'w') as wt_bed:
        wt_bed.write(bed_final)
    return final


def write_output(final, invcf):
    """
    Make tsv format output file.
    """
    with open('{}'.format(invcf.replace('.vcf', '.tsv')), 'w') as wt_fianl:
        wt_fianl.write(final)
    return True


def make_screenshot(recal):
    """Make a screenshot of IGV
    
    Make a screenshot of all examples in cap-pt excel file

    Args:
        recal(file): a recal bam file

    Returns:
        None

    """
    
    snapshot_script = '/clinix/Analysis/Projects/RnD/Certification/CAP-PT/CAP-PT_script/make_IGV_snapshots.py'
    IGV = '/clinix/Tools/Dev/IGV_2.4.4/igv.jar'
    screenshot_cmd = (
            'python {} {} -g hg19 -o ./snapshot -bin {} -r cap-pt.bed'
            .format(snapshot_script, recal, IGV)
            )
    os.system(screenshot_cmd)
    return True


def main():
    """
    A script for parsing snpeff annotated vcf,
    and finding answer to CAP-PT questions.

    Usage:
        python cap_vcf_parsing.py
        -v [vcf]
        -c [excel_file]
    
    Args:
        vcf(file): snpeff annotated vcf file
        excel_file(file): cap-pt entry excel file

    Return:
        tsv(file): tsv format file with correct answer
    
    """
    invcf, c_exc, recal, bed = argument_parse()
    print('--- I. Parse: CAP-PT Excel ---')
    df_exc = c_exc_parse(c_exc)
    print('--- II. Parse: Analysis bed ---')
    bed_dict = bed_parse(bed)
    print('--- III. Parse: Annotated vcf ---')
    invcf_dict, sample_list = invcf_parse(invcf)
    print('--- IV. Process: Answer tsv result ---')
    final = make_answer(df_exc, invcf_dict, bed_dict, sample_list)
    if write_output(final, invcf):
        pass
    print('--- V. Process: IGV Screenshot ---')
    if make_screenshot(recal):
        pass
    print('--- Completed ---')


if __name__=='__main__':
    main()
