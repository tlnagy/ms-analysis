# ms-analysis
Analysis of MS results for PUBS

## Reproducing results

The code can be run by 

```
./run_analysis.py -d ms_data
```

## Development

To build on top of this code, first import it (if it is in the current
directory) like so

```
import run_analysis as ra
```

then make sure to run the `process_data` function and point it to the
location of the ms_data directory containing the data from the MS
experiment e.g.

```
ms_data/
├── Phospho\ (STY)Sites.txt
├── experimentalDesignTemplate.txt
├── parameters.txt
├── peptides.txt
├── proteinGroups.txt
├── summary.txt
└── tables.pdf

0 directories, 7 files
```

so we have 

```
path = "ms_data"
prot_groups, intensity_cols, normed_intensities, log_ratios = run_analysis.process_data(path)
```
where `path` is the location of the ms_data folder relative to the current
directory. And that's it! Now you access to the processed data. The other
functions within `run_analysis` can be accessed similarly
