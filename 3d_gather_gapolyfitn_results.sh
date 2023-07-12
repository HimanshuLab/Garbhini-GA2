#! /bin/bash

gapolyfitn_results_dir=$1;
outdir=$2;
suffix=$3;

grep -R "^Term" $gapolyfitn_results_dir | sed 's/:/\t/g;s/.txt//g;s/_itr/\t/g' >  $2"/gapolyfitn_terms"$3".txt";
grep -R "RMSE"  $gapolyfitn_results_dir | sed 's/:/\t/g;s/.txt//g;s/_itr/\t/g' | awk '{print $1"\t"$2"\t"$4}' > $2"/gapolyfitn_rmse"$3".txt";
grep -R "Var1"  $gapolyfitn_results_dir | grep -v ">x<" | sed 's/<strong>//g;s/<\/strong>//g;s/.txt//g;s/://g;s/Var1/\t/g;s/_itr/\t/g' > $2"/gapolyfitn_variables"$3".txt";

Rscript ./3e_gather_gapolyfitn_results.R $2"/gapolyfitn_terms"$3".txt" $2"/gapolyfitn_rmse"$3".txt" $2"/gapolyfitn_variables"$3".txt" $2"/gapolyfitn_formulas"$3".txt" $2"/gapolyfitn_formulas_consistent"$3".txt"
