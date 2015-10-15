#!/usr/bin/env python
#
# Copyright 2015 Tamas Nagy, Fatima Ugur, Ruilin Tian
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from __future__ import print_function, division, unicode_literals
import argparse, sys, os, re
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy

prot_groups = None
sns.set_style("white")


class HelpParser(argparse.ArgumentParser):
    """
    Like the normal parser except shows the help message on default
    instead of requiring the user to type -h to see it
    """
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def process_data(directory="ms_data"):
    """
    Runs all data processing steps, returns four dataframes

    Returns:
        prot_groups         - original data
        intensity_cols      - raw intensity values (wcl only)
        normed_intensities  - normalized intensity values (wcl only)
        log_ratios          - log fold change ratios vs WT with only proteins >
                              7.75 fold change in at least one experiment present
    """
    prot_groups = pd.read_table(os.path.join(directory, "proteinGroups.txt"),
                                converters={"Reverse": bool, "Potential contaminant": bool})

    psites = pd.read_table(os.path.join(directory, "Phospho (STY)Sites.txt"),
                           converters={"Reverse": bool, "Diagnostic peak": bool,
                                       "Potential contaminant": bool})
    # Remove contaminants, contaminant column doesn't include all
    # reverse proteins
    prot_groups = prot_groups[~prot_groups["Protein IDs"].str.contains("CON__") & ~prot_groups["Protein IDs"].str.contains("REV__")]

    # grab only the columns with the intensity values for the WCL
    # experiments and update the index
    intensity_cols = prot_groups.loc[:, [colname for colname in prot_groups.columns if colname.startswith("Intensity ")
                                         and colname.endswith("_WCL")]]
    intensity_cols.index = prot_groups["Fasta headers"].str.split(" ").str.get(1)
    # set zero intensities to null
    intensity_cols = intensity_cols[intensity_cols != 0]

    # find experiments that had <250 proteins with non-zero intensities
    low_count_exps = intensity_cols.columns[intensity_cols[intensity_cols != 0].count() < 250]

    # Normalize intensities in each experiment by dividing by the sum of
    # all intensities for that experiment
    normed_intensities = intensity_cols/intensity_cols.sum()
    normed_intensities = normed_intensities[list(set(intensity_cols.columns) - set(low_count_exps))]
    normed_intensities.columns = [colname+"_normed" for colname in normed_intensities.columns]
    num_exps = len(["_normed" in colname for colname in normed_intensities.columns])
    ratios_all_present_prots = normed_intensities[normed_intensities.count(axis=1) == num_exps]
    ratios_all_present_prots = ratios_all_present_prots.divide(ratios_all_present_prots["Intensity Control_WCL_normed"],
                                                               axis=0)
    ratios_all_present_prots.drop("Intensity Control_WCL_normed", axis=1, inplace=True)

    log_ratios = np.log2(ratios_all_present_prots)

    # throw out all proteins that have less than a 7.75 fold change in all
    # experiments. This was done to reduce the size of the "significant" protein
    # set while preserving the clustering of experiments together.
    log_ratios = log_ratios[~pd.isnull(log_ratios[np.abs(log_ratios) > np.log2(7.75)].sum(axis=1))]
    log_ratios.columns = re.findall("Intensity ([A-z_0-9]*)_WCL_normed", "\n".join(log_ratios.columns))

    return prot_groups, intensity_cols, normed_intensities, log_ratios


def plot_normed_intensities(normed_intensities, path=""):
    """
    Plot the normalized intensities to show that the distributions
    overlap post-normalization
    """
    for col in normed_intensities.columns:
        sns.kdeplot(np.log2(normed_intensities[col]), label=col, shade=True)
        plt.legend(bbox_to_anchor=(1.3, 0.5), loc="center right")
    sns.despine()
    sns.axlabel("Log$_2$ Normalized Intensities", "Density")
    plt.tight_layout()
    plt.savefig(os.path.join(path, "normed_intensities_wcl.png"), dpi=300)


def plot_scatter_intensities(intensity_cols, experiment1, experiment2="Control_WCL", path=""):
    """
    Create a scatter plot of the per-protein raw intensities between
    two experiments. Experiment names should be the names of the columns
    minus the "Intensity " part.
    """
    if "Intensity %s"%experiment1 not in intensity_cols.columns:
        raise ValueError("Supplied column name not in dataframe")
    f = plt.figure(figsize=(3, 3))
    plt.loglog()
    tmp_prot_groups = intensity_cols.loc[:, ["Intensity %s"%experiment1, "Intensity %s"%experiment2]]
    sns.regplot(tmp_prot_groups.ix[:, 1], tmp_prot_groups.ix[:, 0], fit_reg=False)
    xs = np.linspace(tmp_prot_groups.min()[0], tmp_prot_groups.max()[0])
    plt.plot(xs, tmp_prot_groups.sum()[0]/tmp_prot_groups.sum()[1]*xs, "r", alpha=0.4)
    sns.axlabel("%s Raw Intensities"%experiment1.replace("_", " "),
                "%s Raw Intensities"%experiment2.replace("_", " "))
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(path, "%sv%s_intensities.png"%(experiment1, experiment2)), dpi=300)


def cluster(log_ratios, cophenetic_cutoff=0.5):
    """
    Computes clusters using hierarchical clustering and the given cophenetic
    distance cutoff

    Returns:
        row_linkage
        col_linkage
        cluster_groups, dict: mapping cluster number (integer) to list of
        gene names in that cluster
    """
    dists = scipy.spatial.distance.pdist(log_ratios.values)
    row_linkage = scipy.cluster.hierarchy.linkage(dists, method='complete')
    col_linkage = scipy.cluster.hierarchy.linkage(
        scipy.spatial.distance.pdist(log_ratios.values.T), method='complete')
    log_ratios["cluster"] = scipy.cluster.hierarchy.fcluster(row_linkage, cophenetic_cutoff*dists.max(), 'distance')
    return row_linkage, col_linkage, {name: group.index.tolist() for name, group in log_ratios.groupby("cluster")}


def plot_heatmap(df, row_linkage=None, col_linkage=None, font_size=10, legend_title="", path=""):
    """
    Plots a given df using the provided row_linkage and col_linkage output from
    scipy's linkage function. The other parameters are for adjusting the plot
    formatting
    """
    cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    # ignore genes that have log ratio changes less than 1 for all experiments
    g = sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage, cmap=cmap, figsize=(20, 20))
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=30, horizontalalignment='right', fontsize=font_size)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=font_size)
    plt.setp(g.cax.yaxis.get_majorticklabels(), fontsize=font_size)

    # fix ugly positioning
    hm_pos, rd_pos, cd_pos = g.ax_heatmap.get_position(), g.ax_row_dendrogram.get_position(), \
                             g.ax_col_dendrogram.get_position()
    g.ax_heatmap.set_position([hm_pos.x0-rd_pos.width*.25, hm_pos.y0, hm_pos.width*0.5, hm_pos.height])
    g.ax_row_dendrogram.set_position([rd_pos.x0, rd_pos.y0, rd_pos.width*0.75, rd_pos.height])
    g.ax_col_dendrogram.set_position([cd_pos.x0-rd_pos.width*.25, cd_pos.y0, cd_pos.width*0.5, cd_pos.height/2])
    cax_pos = g.cax.get_position()
    g.cax.set_position([cax_pos.x0*.75, cax_pos.y0*.25, cax_pos.width*0.5, cax_pos.height*0.5])
    g.cax.set_title(legend_title, fontsize=font_size)

    plt.savefig(os.path.join(path, "sig_heatmap.png"), dpi=300, bbox_inches="tight")


def main(directory):
    # import data as pandas dataframes, the contaminant columns are
    # converted to booleans so that the pandas parser is faster and less
    # grouchy
    prot_groups, intensity_cols, normed_intensities, log_ratios = process_data(directory)

    row_linkage, col_linkage, flat_clusters = cluster(log_ratios)
    # plot
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = "plots"
        try:
            os.makedirs(path)
        except OSError:
            if not os.path.isdir(path):
                raise
        plot_normed_intensities(normed_intensities, path=path)
        plot_scatter_intensities(intensity_cols, experiment1="EtOH_Kin3KO_WCL", path=path)
        plot_scatter_intensities(intensity_cols, experiment1="EtOH_menadione_WCL", path=path)
        plot_heatmap(log_ratios, row_linkage=row_linkage, col_linkage=col_linkage,
                     legend_title="log normed\nintensity\nchange vs WT", path=path)

if __name__ == "__main__":
    parser = HelpParser(description='Analyze MS data for the Et0h group')
    parser.add_argument('-d', '--dir', help='Path to the directory\
            containing the MaxQuant files', required=True)
    args = vars(parser.parse_args())
    
    if os.path.isdir(args["dir"]):
        main(args["dir"])
    else:
        sys.stderr.write('error: Not a valid directory\n')
        sys.exit(2)
