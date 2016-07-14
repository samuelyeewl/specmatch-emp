#!/bin/bash
set -x
seq 5100 100 5500 | parallel python tests/library_match.py "./lib/library_reduced.h5" "./results/reduced/{1}_results.csv" {1} 100
seq 5100 100 5500 | awk '{print $1"\t"$1+100}' | parallel --colsep '\t' python tests/match_analysis.py "./lib/library_reduced.h5" "./results/reduced/{1}_results.csv" "/Users/samuel/Dropbox/SpecMatch-Emp/plots/lincomb/" "{1}" "'Wavelength: {1} - {2} A'"
seq 5100 100 5500 | parallel python tests/library_match_lincomb.py "./lib/library_reduced.h5" "./results/reduced/{1}_results.csv" "./results/reduced/lincomb5_{1}_results.csv" 5 {1} 100
seq 5100 100 5500 | awk '{print $1"\t"$1+100}' | parallel --colsep '\t' python tests/match_lincomb_analysis.py "./lib/library_reduced.h5" "./results/reduced/lincomb5_{1}_results.csv" "/Users/samuel/Dropbox/SpecMatch-Emp/plots/lincomb/" "{1}_lincomb5" "'Wavelength: {1} - {2} A, Linear Combination: Best 5'"
seq 2 8 | parallel python tests/library_match_lincomb.py "./lib/library_reduced.h5" "./results/reduced/5300_results.csv" "./results/reduced/lincomb{1}_5300_results.csv" {1} 5300 100
seq 2 8 | parallel python tests/match_lincomb_analysis.py "./lib/library_reduced.h5" "./results/reduced/lincomb{1}_5300_results.csv" "/Users/samuel/Dropbox/SpecMatch-Emp/plots/lincomb/" "5300_lincomb{1}" "'Wavelength: 5300 - 5400 A, Linear Combination: Best {1}'"