#!/bin/sh

sample_dir=$1 #eg. 'TIGER/BSP-OR_001'
file_pattern=$2 #eg. 'BSP-OR_001'

cd $sample_dir



################################

# run base caller

java -jar ../../TIGER_scripts/base_caller.jar -r $file_pattern'.tiger_input.txt' -o $file_pattern".base_caller.txt" -n bi

################################

# run allele frequency estimator

java -jar ../../TIGER_scripts/allele_freq_estimator.jar -r $file_pattern'.tiger_input.txt' -o $file_pattern".bmm_allele_freqs.txt" -n bi -w 1000

################################

# apply Beta mixture model

Rscript --vanilla ../../TIGER_scripts/beta_mixture_model.R $file_pattern".bmm_allele_freqs.txt" $file_pattern".bmm_intersections.txt"

################################

# prepare files for HMM probability estimation

perl ../../TIGER_scripts/prep_hmm_probs.pl -s $file_pattern -m $file_pattern'.tiger_input.txt' -b $file_pattern".base_caller.txt" -c ../../TIGER_scripts/n2_chrom_lengths.tiger.txt -o $file_pattern'.tiger.hmm_prep_prob.txt'

################################

# calculate transmission and emission probabilities for HMM

perl ../../TIGER_scripts/hmm_probs.pl -s $file_pattern".bmm_allele_freqs.txt" -p $file_pattern'.tiger.hmm_prep_prob.txt' -a $file_pattern".bmm_intersections.txt" -c ../../TIGER_scripts/n2_chrom_lengths.tiger.txt -o $file_pattern'.'

################################

# Run HMMs

java -jar ../../TIGER_scripts/hmm_play.jar -r $file_pattern".base_caller.txt" -z $file_pattern'._hmm_model' -t bi -o $file_pattern'.hmm_output.txt'

################################

# rough breakpoint estimation

perl ../../TIGER_scripts/estimate_breaks.pl -s $file_pattern -m $file_pattern'.tiger_input.txt' -b $file_pattern'.hmm_output.txt' -c ../../TIGER_scripts/n2_chrom_lengths.tiger.txt -o $file_pattern'.CO_estimates.txt'

################################

# refine breakpoints

perl ../../TIGER_scripts/refine_recombination_break.pl $file_pattern'.tiger_input.txt' $file_pattern'.CO_estimates.breaks.txt'

################################

# smooth refined breakpoints

perl ../../TIGER_scripts/breaks_smoother.pl -b $file_pattern'.CO_estimates.refined.breaks.txt' -o $file_pattern'.smoothed.breaks.txt'
