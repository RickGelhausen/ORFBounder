#!/usr/bin/env python
import os,sys
import re
import argparse
import collections
import csv

import pandas as pd
import numpy as np

from pathlib import Path

import lib.misc as misc
import lib.expression as expr
import lib.io as io

import matplotlib.pyplot as plt
from statsmodels.graphics.gofplots import qqplot
from scipy.stats import shapiro
from scipy.stats import anderson
from scipy.stats import normaltest
import scipy
import argparse

from sklearn.preprocessing import StandardScaler
import scipy.stats

def read_data(data_file):
    with open(data_file, "r") as f:
        data = [float(x.rstrip("\n")) for x in f.readlines()]

    return data

def get_common_dists():
    dist_names = ["weibull_min","norm","weibull_max","beta",\
              "invgauss","uniform","gamma","expon",\
              "lognorm","pearson3","triang"]
    return dist_names

def get_all_dists():
    dist_names = ["weibull_min","norm","weibull_max","beta",\
              "invgauss","uniform","gamma","expon",\
              "lognorm","pearson3","triang"]

    all = ["alpha", "anglit", "arcsine", "betaprime", "bradford", "burr", "cauchy", "chi", "chi2", "cosine", "dgamma", \
    "dweibull", "exponweib", "exponpow", "f", "fatiguelife", "fisk", "foldcauchy", "foldnorm",\
    "gausshyper", "genexpon", "genextreme", "gengamma", "genhalflogistic", "genlogistic", "genpareto", "gilbrat", \
    "gompertz", "gumbel_l", "gumbel_r", "halfcauchy", "halflogistic", "halfnorm", "hypsecant", "invgamma", "invweibull", \
    "johnsonsb", "johnsonsu", "laplace", "logistic", "loggamma", "loglaplace", "lomax", "maxwell", "mielke", "nakagami","ncx2", \
    "ncf", "nct", "pareto", "powerlaw", "powerlognorm", "powernorm", "rdist", "reciprocal", "rayleigh", "rice", "recipinvgauss", \
    "semicircular", "t", "truncexpon", "truncnorm", "tukeylambda", "vonmises", "wald", "wrapcauchy", "ksone", "kstwobign"]

    # problematic = ["erlang", ]
    dist_names.extend(all)
    return dist_names
    # distributions = []
    # for this in dir(scipy.stats):
    #     if "fit" in eval("dir(scipy.stats." + this + ")"):
    #         distributions.append(this)
    # return distributions

def chi_square_method(data, size, dist_names, output_path, card):
    """
    """

    chi_square = []
    p_values = []

    # 11 equi-distant bins of observed Data
    percentile_bins = np.linspace(0,100,51)
    percentile_cutoffs = np.percentile(data, percentile_bins)
    observed_frequency, bins = (np.histogram(data, bins=percentile_cutoffs))
    cum_observed_frequency = np.cumsum(observed_frequency)

    # Loop through candidate distributions
    for distribution in dist_names:
        # Set up distribution and get fitted distribution parameters
        dist = getattr(scipy.stats, distribution)
        param = dist.fit(data)

        # p-value calculation
        p = scipy.stats.kstest(data, distribution, args=param)[1]
        p = np.around(p, 5)
        p_values.append(p)


        # Get expected counts in percentile bins
        # cdf of fitted sistrinution across bins
        cdf_fitted = dist.cdf(percentile_cutoffs, *param[:-2], loc=param[-2], scale=param[-1])
        expected_frequency = []
        for bin in range(len(percentile_bins)-1):
            expected_cdf_area = cdf_fitted[bin+1] - cdf_fitted[bin]
            expected_frequency.append(expected_cdf_area)

        # Chi-square Statistics
        expected_frequency = np.array(expected_frequency) * size
        cum_expected_frequency = np.cumsum(expected_frequency)
        # ss = sum ((( cum_observed_frequency -cum_expected_frequency) ** 2) / cum_observed_frequency)
        ss = scipy.stats.chisquare(cum_observed_frequency, cum_expected_frequency)
        chi_square.append(ss)


    #Sort by minimum ch-square statistics
    results = pd.DataFrame()
    results['Distribution'] = dist_names
    results['chi_square'] = chi_square
    results['p_values'] = p_values
    results.sort_values(['chi_square'], inplace=True)

    with open(os.path.join(output_path, "%s_top_distributions.tsv" % card), "w") as f:
        results.to_csv(f, sep="\t")

    return results

def standardize_data(y):
    """
    Standardize data
    """
    sc=StandardScaler()
    yy = y.reshape(-1,1)
    sc.fit(yy)
    y_std =sc.transform(yy)
    y_std = y_std.flatten()
    return y_std

def plot_best_n_distributions(y, results, bestn, output_path, card):
    """
    Create a plot with the n best distributions fitted to the histogram
    """

    x = np.arange(len(y))

    number_of_bins = 100
    bin_cutoffs = np.linspace(np.percentile(y,0), np.percentile(y,99),number_of_bins)

    h = plt.hist(y, bins = bin_cutoffs, color='0.75')
    dist_names = results['Distribution'].iloc[0:bestn]

    parameters = []

    for dist_name in dist_names:
        dist = getattr(scipy.stats, dist_name)
        param = dist.fit(y)
        parameters.append(param)

        pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])
        scale_pdf = np.trapz(h[0], h[1][:-1]) / np.trapz(pdf_fitted, x)
        pdf_fitted *= scale_pdf

        plt.plot(pdf_fitted, label=dist_name)
        plt.xlim(0,np.percentile(y, 99))

    plt.legend()
    plt.savefig(os.path.join(output_path, "%s_top_distributions.pdf" % card))
    plt.close()
    # Store distribution paraemters in a dataframe (this could also be saved)
    dist_parameters = pd.DataFrame()
    dist_parameters['Distribution'] = (results['Distribution'].iloc[0:bestn])
    dist_parameters['Distribution parameters'] = parameters

    # Print parameter results
    with open(os.path.join(output_path, "%s_top_distributions.txt" % card), "w") as f:
        f.write('Distribution parameters:\n')
        f.write('------------------------\n')

        for index, row in dist_parameters.iterrows():
            f.write("Distribution:" + row[0] + "\n")
            f.write("Parameters:" + ",".join([str(v) for v in list(row[1])]) + "\n")

def qq_pp_plot(y_std, results, bestn, size, output_path, card):
    """
    Create a
    """
    data = y_std.copy()
    data.sort()

    dist_names = results['Distribution'].iloc[0:bestn]

    # Loop through selected distributions (as previously selected)

    for distribution in dist_names:
        # Set up distribution
        dist = getattr(scipy.stats, distribution)
        param = dist.fit(y_std)

        # Get random numbers from distribution
        norm = dist.rvs(*param[0:-2],loc=param[-2], scale=param[-1],size = size)
        norm.sort()

        # Create figure
        fig = plt.figure(figsize=(8,5))

        # qq plot
        ax1 = fig.add_subplot(121) # Grid of 2x2, this is suplot 1
        ax1.plot(norm,data,"o")
        min_value = np.floor(min(min(norm),min(data)))
        max_value = np.ceil(max(max(norm),max(data)))
        ax1.plot([min_value,max_value],[min_value,max_value],'r--')
        ax1.set_xlim(min_value,max_value)
        ax1.set_xlabel('Theoretical quantiles')
        ax1.set_ylabel('Observed quantiles')
        title = 'qq plot for ' + distribution +' distribution'
        ax1.set_title(title)

        # pp plot
        ax2 = fig.add_subplot(122)

        # Calculate cumulative distributions
        bins = np.percentile(norm,range(0,101))
        data_counts, bins = np.histogram(data,bins)
        norm_counts, bins = np.histogram(norm,bins)
        cum_data = np.cumsum(data_counts)
        cum_norm = np.cumsum(norm_counts)
        cum_data = cum_data / max(cum_data)
        cum_norm = cum_norm / max(cum_norm)

        # plot
        ax2.plot(cum_norm,cum_data,"o")
        min_value = np.floor(min(min(cum_norm),min(cum_data)))
        max_value = np.ceil(max(max(cum_norm),max(cum_data)))
        ax2.plot([min_value,max_value],[min_value,max_value],'r--')
        ax2.set_xlim(min_value,max_value)
        ax2.set_xlabel('Theoretical cumulative distribution')
        ax2.set_ylabel('Observed cumulative distribution')
        title = 'pp plot for ' + distribution +' distribution'
        ax2.set_title(title)

        # Display plot
        plt.tight_layout(pad=4)
        plt.savefig(os.path.join(output_path, "%s_qq_pp_plots_%s.pdf" % (card, distribution)))
        plt.close()

def read_input_table(input_table):
    """
    Read .xlsx table and return a dataframe
    """

    return pd.read_excel(input_table, sheet_name=None)["CDS"]

# def create_idr_format_data_frame(xlsx_df, replicate):
#     """
#     create a dataframe for a given replicate, dest="genome"
#     """
#
#     replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
#     non_header = ["chromosome", "start", "stop", "Name", "score", "strand", "signalValue", "pvalue", "qvalue", "summit"]
#     nTuple = collections.namedtuple('Pandas', non_header)
#
#     rows = []
#     for row in xlsx_df.itertuples(index=False, name=None):
#         genome, start, stop, strand = row[2], int(row[3]), int(row[4]), row[5]
#         peak_height = row[replicate_index]
#         if np.isnan(peak_height) or peak_height == 0:
#             continue
#
#         # TODO MISSING TIS/TTS SPECIFICATIONS
#         if strand == "+":
#             out_start, out_stop = start-1, start+1
#         else:
#             out_start, out_stop = stop-3, stop-1
#
#         result = [genome, out_start, out_stop, ".", 0, strand, float(peak_height), -1, -1, -1]
#
#         rows.append(nTuple(*result))
#
#     return pd.DataFrame.from_records(rows, columns=non_header)
#
# def create_gff_format_data_frame(xlsx_df, replicate):
#     """
#     create a dataframe for a given replicate
#     """
#
#     replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
#     non_header = ["chrom", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"]
#     nTuple = collections.namedtuple('Pandas', non_header)
#
#     rows = []
#     for row in xlsx_df.itertuples(index=False, name=None):
#         unique_id, genome, start, stop, strand = row[1], row[2], int(row[3]), int(row[4]), row[5]
#         peak_height = row[replicate_index]
#         if np.isnan(peak_height) or peak_height == 0:
#             continue
#         # TODO MISSING TIS/TTS SPECIFICATIONS
#         if strand == "+":
#             out_start, out_stop = start-1, start+1
#         else:
#             out_start, out_stop = stop-3, stop-1
#
#         result = [genome, "ORFBounder", "TSS", out_start, out_stop, peak_height, strand, ".", "ID=%s;" % (unique_id)]
#
#         rows.append(nTuple(*result))
#
#     return pd.DataFrame.from_records(rows, columns=non_header)
#
# def retrieve_replicates(xlsx_df):
#     """
#     Return a list of the replicates in the input file
#     """
#
#     return [entry[:-12] for entry in xlsx_df.columns if "_peak_height" in entry]
#
# # def write_dataframe_to_file(cur_df, tmp_folder, replicate):
# #     """
# #     write peak_height dataframe to file
# #     """
# #
# #     Path(os.path.dirname(tmp_folder)).mkdir(parents=True, exist_ok=True)
# #     cur_df.to_csv(os.path.join(tmp_folder, replicate + "_peak_height.gff"), sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
#
# def write_dataframe_to_file(cur_df, tmp_folder, replicate):
#     """
#     write peak_height dataframe to file
#     """
#
#     Path(os.path.dirname(tmp_folder)).mkdir(parents=True, exist_ok=True)
#     cur_df.to_csv(os.path.join(tmp_folder, replicate + "_peaks"), sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
#
# def generate_peak_height_output_files(xlsx_df, tmp_folder):
#     """
#     Create an idr input file for each replicate in the input file
#     """
#
#     for replicate in retrieve_replicates(xlsx_df):
#         cur_df = create_idr_format_data_frame(xlsx_df, replicate)
#
#         write_dataframe_to_file(cur_df, tmp_folder, replicate)

def get_wildcards(xlsx_df):
    """
    retrieve wildcards from _peak_height columns
    """

    return [entry[:-12] for entry in xlsx_df.columns if "_peak_height" in entry]


def check_bamfile_input(bam_file_path, wildcards):
    """
    Check bam input path.
    Ensure that there is:
     - one bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if bam_file_path == "" or not isinstance(bam_file_path, str):
        return -1

    valid_bam = []

    _, _, bam_file_list = next(os.walk(bam_file_path))
    bam_file_list = [ file for file in bam_file_list if file.endswith(".bam")]

    for card in wildcards:

        if "TIS" in card:
            #RNATIS_prefix = "RNATIS-" + "-".join(card.split("-")[1:])
            for file in bam_file_list:
                if card in file and not "RNA" in file:# or RNATIS_prefix in file:
                    valid_bam.append(os.path.join(bam_file_path,file))

        if "TTS" in card:
            #RNATTS_prefix = "RNATTS-" + "-".join(card.split("-")[1:])
            for file in bam_file_list:
                if card in file and not "RNA" in file:# or RNATTS_prefix in file:
                    valid_bam.append(os.path.join(bam_file_path,file))

    if len(valid_bam) == 0:
        return -1

    return valid_bam

def retrieve_original_offset_position(start, stop, strand, offset, method):
    """
    retrieve the original position of the peaks
    """

    original_position = 0
    if method == "TIS":
        if strand == "+":
            original_position = start + offset
        else:
            original_position = stop - offset

    else:
        if strand == "+":
            original_position = stop + offset
        else:
            original_position = start - offset

    return original_position

def create_replicate_peak_dict(xlsx_df, replicate, offset_dict):
    """
    create dictionary for a given replicate with the original positions of the peak.
    """

    replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
    if "TIS" in replicate:
        method = "TIS"
    else:
        method = "TTS"

    replicate_dict = {}
    for row in xlsx_df.itertuples(index=False, name=None):
        unique_id, chrom, start, stop, strand = row[1], row[2], int(row[3]), int(row[4]), row[5]
        peak_height = row[replicate_index]

        if np.isnan(peak_height) or peak_height == 0:
            continue

        try:
            offset = offset_dict[method][replicate]
        except KeyError:
            offset = offset_dict[method]["default"]

        original_position = retrieve_original_offset_position(int(start)-1, int(stop)-1, strand, int(offset), method)

        replicate_dict[(chrom, original_position-2, original_position+2, strand)] = (unique_id, peak_height, 0, 0)

    return replicate_dict


def create_dataframe(read_count_dict, replicate):
    """
    create dataframe with results
    """

    header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Original_peak_pos"] \
           + [replicate + "_peak_height", replicate + "_peak_rpkm", replicate + "_read_count"]

    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for (chrom, interval_start, interval_stop, strand), (unique_id, peak_height, rpkm, read_count) in read_count_dict.items():
        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        result_rows.append(nTuple(unique_id, chrom, start, stop, strand, interval_start+2, peak_height, rpkm, read_count))

    return pd.DataFrame.from_records(result_rows, columns=header)

def write_test_file(out_df, output_path):
    """
    write a test file
    """
    Path(os.path.dirname(output_path)).mkdir(parents=True, exist_ok=True)

    io.excel_writer(output_path, {"CDS" : out_df})

def sharpiro_wilk_test(data):
    """
    Perform shapiro wilk test and print results
    """
    print("Performing shapiro wilk test:")
    stat, p = shapiro(data)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
    	print('Sample looks Gaussian (fail to reject H0)')
    else:
    	print('Sample does not look Gaussian (reject H0)')

def d_agostino_test(data):
    """
    Perform d agostinos K2 test and print results
    """
    print("Performing agostino K2 test:")
    stat, p = normaltest(data)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
    	print('Sample looks Gaussian (fail to reject H0)')
    else:
    	print('Sample does not look Gaussian (reject H0)')

def anderson_darling_test(data, output_path, card):
    """
    """
    print("Performing Anderson-darling test:")
    result = anderson(data)
    print('Statistic: %.3f' % result.statistic)
    p = 0
    for i in range(len(result.critical_values)):
    	sl, cv = result.significance_level[i], result.critical_values[i]
    	if result.statistic < result.critical_values[i]:
    		print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
    	else:
    		print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))

def print_list_file(data, file_path):
    with open(file_path, "w") as f:
        for entry in data:
            f.write("%s\n" % entry)

def analyse_distributions_and_plot(data, output_path, card):
    """
    """

    y = np.array([x for x in data])# if x <= 20000])

    y_std = standardize_data(y)

    df = pd.DataFrame(list(zip(y,y_std)), columns=["y", "y_std"])
    print(df.describe())

    size = len(y)
    # fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9,5))

    # axes[0].hist(data)
    # axes[1].plot(sorted(data), stats.lognorm.pdf(data, 1.4932397758225622, -0.7404025861595065, 252.61459557262953))
    # axes[2].plot(data, stats.invgauss.pdf(data, 4.883726683673663, -14.814193193403817, 164.3724468978952))
    # plt.savefig("/mnt/datavault/SPP2002/analysis/test_ORFBounder/overview_dist.png")
    #

    results = chi_square_method(y_std, size, get_all_dists(), output_path, card)
    plot_best_n_distributions(y, results, 2, output_path, card)
    qq_pp_plot(y, results, 2, size, output_path, card)


def check_bamfile_input_RIBO(bam_file_path, wildcards):
    """
    Check bam input path.
    Ensure that there is:
     - one bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if bam_file_path == "" or not isinstance(bam_file_path, str):
        return -1

    valid_bam = []

    _, _, bam_file_list = next(os.walk(bam_file_path))
    bam_file_list = [file for file in bam_file_list if file.endswith(".bam")]

    for card in wildcards:
        ribo = "RIBO-%s-%s" % tuple(card.split("-")[1:])

        for file in bam_file_list:
            if ribo in file:
                valid_bam.append(os.path.join(bam_file_path,file))

    if len(valid_bam) == 0:
        return -1

    return valid_bam

def create_replicate_orf_dict(xlsx_df, replicate):
    """
    create dictionary for a given replicate with the original positions of the peak.
    """

    replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")

    replicate_dict = {}
    for row in xlsx_df.itertuples(index=False, name=None):
        unique_id, chrom, start, stop, strand = row[1], row[2], int(row[3]), int(row[4]), row[5]
        peak_height = row[replicate_index]

        if np.isnan(peak_height) or peak_height == 0:
            continue

        replicate_dict[(chrom, start-1, stop-1, strand)] = (unique_id, peak_height, 0, 0)

    return replicate_dict

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Filter the final output file.')
    parser.add_argument("-i", "--input_table", action="store", dest="input_table", required=True, help= "Table created by ORFBounder merge or call_ORFBounder.py (Requires atleast 2 replicates).")
    parser.add_argument("--offset_json", action="store", dest="offset_json", required=True, help="A JSON file containing all information about the offsets.")
    parser.add_argument("--mapping_tis", action="store", dest="mapping_tis", default="", help="The mapping method used for the TIS data.")
    parser.add_argument("--mapping_tts", action="store", dest="mapping_tts", default="", help="The mapping method used for the TTS data.")
    parser.add_argument("--genome", action="store", dest="genome", required=True, help="The genome file.")
    parser.add_argument("--bamfiles", action="store", dest="bam_files", required=True)
    # parser.add_argument("-t","--tmp_folder", action="store", dest="tmp_folder", required=True, help="folder for storing temporary files.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    xlsx_df = read_input_table(args.input_table)
    offset_dict = misc.build_offset_dictionary(args.offset_json, args.genome, args.mapping_tis, args.mapping_tts)

    wildcards = get_wildcards(xlsx_df)
    bam_files = list(check_bamfile_input(args.bam_files, wildcards))
    print(bam_files)
    for idx, card in enumerate(wildcards):
        print(idx, card)
        read_count_dict = create_replicate_peak_dict(xlsx_df, card, offset_dict)
        interlap_dict, total_mapped = expr.create_interlap_dict(bam_files[idx])
        for (chrom, start, stop, strand), val in read_count_dict.items():

            read_count = expr.count_reads(chrom, start, stop, strand, interlap_dict)
            rpkm = expr.calculate_rpkm(total_mapped[chrom], read_count, int(stop)-int(start)+1)

            read_count_dict[(chrom,start,stop,strand)] = (val[0], val[1], rpkm, read_count)

        df = create_dataframe(read_count_dict, card)

        analyse_distributions_and_plot([x for x in df.loc[df["Strand"] == "+", "%s_peak_rpkm" % card].to_list()], args.output_path, card)

if __name__ == '__main__':
    main()
