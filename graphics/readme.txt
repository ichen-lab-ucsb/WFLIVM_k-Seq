This folder contains scripts for processing k-seq fitting results to produce files containing all merged data and additional analyses like sequence family, catalytic enhancement, and promiscuity index to be used for subsequent analysis. Additional information and scripts can be found at https://github.com/ichen-lab-ucsb/WFLIVM_k-Seq.

— Required inputs:
bfo-results.csv, blo-results.csv, bwo-results.csv, bio-results.csv, bvo-results.csv, bmo-results.csv

— MergeSubs_SplitFams.py: script to merge all bxo-results.csv files into a single file. Also determines the Hamming distance of each sequence to each of five wild type sequences and assigns each sequence to associated sequence families. Assigns a color_label to each family for subsequent plotting. Outputs a merged csv file required for HistogramFits.py and BackgroundNormalizer.py.

— HistogramFits.py: script to plot histograms of all data with each bxo substrate and fit to a bimodal Gaussian distribution of log10-transformed kA_50% values and reports fitting parameters. Outputs a png file of plots and csv files of fitting parameters for each substrate.

— BackgroundNormalizer.py: calculates catalytic enhancement (r) values for each sequence using kA values and background rate estimates, as well as for 2.5% and 97.5% kA estimates. Calculates sums, harmonic means, and ratios of catalytic enhancement values. Outputs a csv file containing these values appended to the input file.

— PromIndex_TableMaker.py: script to make csv files containing sequence and activity data to be used as inputs for calculating promiscuity indices using the online tool available at http://hetaira.herokuapp.com/. Requires the output file from BackgroundNormalizer.py.

— PromIndexResultsMerger.py: script to concatenate and merge the csv files produced as outputs from the online tool available at http://hetaira.herokuapp.com/. The promiscuity index values (I) are merged to the file produced from BackgroundNormalizer.py.