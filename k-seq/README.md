# WFLIVM_k-Seq

This folder contains the script to process count files output from
[EasyDIVER](https://github.com/ichen-lab-ucsb/EasyDIVER) and prepare for k-seq fitting.

### Count data preprocessing

`count-data-preprocessing.py` script contains the data processing pipeline to parse count files,
filter sequences, quantify sequence amount and calculate reacted fraction for *k*-Seq fitting.
'k-seq' package is required for the script.

Use 'count-data-preprocessing.py -h' for more details.

### *k*-Seq fitting
We used [k-seq](https://github.com/ichen-lab-ucsb/k-seq) to fit the preprocessed data to kinetic model with
bootstrapping. Please follow the package instruction for installation.

The fitting was performed using the script 'k-seq-fitting.py' in the k-seq package.
An example on fitting BFO data:

```shell
k-seq-fitting.py \
    --seq_data bfo-preped.pkl \      # output from count-data-preprocessing.py
    --table_name reacted_fraction \  # table to fit
    --bootstrap_num 1000 \           # number of bootstrap resampling
    --bs_method data \               # bootstrap by resampling the data points
    --large_data \                   # save data during fitting
    --core_num 12 \                  # run 12 threads in parallel for fitting
    --output bfo-output              # folder for output files
```
