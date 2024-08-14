#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


for chrom in range(1,38+1):
    # load the output file from imputation_accuracy.py script, which has a header and is tab-separated
    array_imputed_differences = pd.read_csv(
        f"173karray_173kimputed/differences_chr{chrom}_173karray_173kimputed.tsv",
        sep="\t", header="infer", index_col=0
    )
    print(f'chr{chrom}')

    # count the number of rows from the XX sites that do not have a missing value (that is, NA)
    maximum_possible_correct = array_imputed_differences.count(axis=0)
    print(maximum_possible_correct)
    # calculating the proportion correct per dog
    proportion_correctly_imputed_per_dog = (array_imputed_differences == 0).sum(axis=0) / maximum_possible_correct
    # to view this proportion correct per dog in python:
    print(pd.DataFrame(proportion_correctly_imputed_per_dog).to_string(header=False, index=True))
    # to save this matrix to csv
    proportion_correctly_imputed_per_dog.to_csv(f'173karray_173kimputed/chr{chrom}_concordPerDog_173a173i.csv')
    
    # count the number of columns from the XX sites that do not have a missing value (that is, NA)
    maximum_possible_correct_col = array_imputed_differences.count(axis=1)
    print(maximum_possible_correct_col)
    # calculating the proportion correct per dog
    proportion_correctly_imputed_per_site = (array_imputed_differences == 0).sum(axis=1) / maximum_possible_correct_col
    # to view this proportion correct per dog in python:
    print(pd.DataFrame(proportion_correctly_imputed_per_site).to_string(header=False, index=True))
    # to save this matrix to csv
    proportion_correctly_imputed_per_site.to_csv(f'173karray_173kimputed/chr{chrom}_concordPerSite_173a173i.csv')
           
    # plot a histogram for impuatation accuracy per dog
    plt.figure()
    plt.hist(proportion_correctly_imputed_per_dog);
    plt.title(f"Chr{chrom} 173k array v 173k imputed")
    plt.xlabel("Proportion of concordant sites")
    plt.ylabel("Number of dogs");
    plt.savefig(f'173karray_173kimputed/chr{chrom}_concordPerDog_173a173i.pdf')
