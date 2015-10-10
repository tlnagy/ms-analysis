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

def main(directory):
    # import data as pandas dataframes, the contaminant columns are
    # converted to booleans so that the pandas parser is faster and less
    # grouchy
    psites = pd.read_table(os.path.join(directory, "Phospho (STY)Sites.txt"), 
            converters={677: bool, 211:bool, 678:bool})
    prot_groups = pd.read_table(os.path.join(directory, "proteinGroups.txt"), 
            converters={421: bool, 422:bool})

    # TODO: Actual implement some kind of analysis here
    yield

if __name__ == "__main__":
    parser = HelpParser(description='Analyze MS data for the Et0h group')
    parser.add_argument('-d','--dir', help='Path to the directory\
            containing the MaxQuant files', required=True)
    args = vars(parser.parse_args())
    
    if os.path.isdir(args["dir"]):
        main(args["dir"])
    else:
        sys.stderr.write('error: Not a valid directory\n')
        sys.exit(2)
