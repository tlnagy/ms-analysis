#!/usr/bin/env python
#
# Copyright 2015 Tamas Nagy, Fatima Ugur, Ruilin Tian
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from __future__ import print_function, division, unicode_literals
import argparse, sys, os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class HelpParser(argparse.ArgumentParser):
    """
    Like the normal parser except shows the help message on default
    instead of requiring the user to type -h to see it
    """
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def plot_normed_intensities(normed_intensities):
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
    plt.savefig("normed_intensities_wcl.png", dpi=300)


def plot_scatter_intensities(intensity_cols, experiment1, experiment2="Control_WCL"):
    """
    Create a scatter plot of the per-protein raw intensities between
    two experiments. Experiment names should be
    """
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
    plt.savefig("%sv%s_intensities.png"%(experiment1, experiment2), dpi=300)


def main(directory):
    # import data as pandas dataframes, the contaminant columns are
    # converted to booleans so that the pandas parser is faster and less
    # grouchy
    psites = pd.read_table(os.path.join(directory, "Phospho (STY)Sites.txt"),
                           converters={"Reverse": bool, "Diagnostic peak": bool,
                                       "Potential contaminant": bool})
    prot_groups = pd.read_table(os.path.join(directory, "proteinGroups.txt"),
                                converters={"Reverse": bool, "Potential contaminant": bool})

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

    # plot
    plot_normed_intensities(normed_intensities)



    # TODO: Actual implement some kind of analysis here
    yield

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
