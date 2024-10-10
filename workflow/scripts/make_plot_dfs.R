library(tidyverse)
library(ggthemes)

hmm_states <- read_csv(snakemake@input[['states']])
hmm_intervals <- read_csv(snakemake@input[['intervals']])
interval_states <- list(hmm_intervals, hmm_states) %>% map(
  \(x) x %>% mutate(
    chromosome = case_match(
      chromosome,
      1 ~ 'chr2L',
      2 ~ 'chr2R',
      3 ~ 'chr3L',
      4 ~ 'chr3R',
      5 ~ 'chrX'
      )
    ) %>% group_by(chromosome))
# split on chromosome because we need to apply different 
# breakpoints to each chromosome
intvls <- interval_states[[1]] %>% select(!chrom_id)
sts <- interval_states[[2]] %>% group_split()
# chromosomes aren't listed in the intervals df if they didn't have a CO 
# so we have to grab the hmm state from those to plot them
chromosomes <- tibble(
  'chromosome' = list('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX')
  )
distinct_chrs <- intvls %>% distinct(chromosome) %>% arrange(chromosome)
missing_chrs <- chromosomes %>%
  filter(! chromosome %in% distinct_chrs$chromosome)
smoothed_st_values <- sts %>% map(
  \(x) x %>%
  filter(chromosome %in% missing_chrs$chromosome) %>%
  select(sample, chromosome, hmm_state1) %>% head(n=1)
  )

chr_ends <- read_tsv(
  snakemake@input[['chrom_lengths']], col_names=c('chromosome', 'length')
  )
chr_ends <- chr_ends %>% mutate(
  chromosome = case_match(
    chromosome,
    1 ~ 'chr2L',
    2 ~ 'chr2R',
    3 ~ 'chr3L',
    4 ~ 'chr3R',
    5 ~ 'chrX'
    )
  )

# make a column that has the breakpoints between hmm states
smoothed_st_values <- smoothed_st_values %>%
  map(
    \(x) x %>% left_join(chr_ends) %>% mutate(start = 1) %>%
      mutate(stop = length) %>% rename(hmm_state = hmm_state1)
    )

st_joined <- smoothed_st_values %>% bind_rows()
intvls <- interval_states[[1]] %>%
  bind_rows(smoothed_st_values) %>% group_by(chromosome) %>% group_split()
brks <- map(
  intvls, \(x) c(1, x[['stop']]) %>% as_tibble_col(column_name='breaks')
  )

# use breakpoints to assign smoothed hmm state to each of the 
# variants in the states df
brks <- brks %>% map2(
  intvls, \(x, y) x %>% mutate(hmm_state = c(y[['hmm_state']], NA))
  )
# the labels list needs to be shorter than the breaks by 1 because it's
# labeling each of the intervals between the breaks
# so lob off the last label, which is NA anyway because it's the last
# position on the chromosome
brks <- brks %>% map(\(x) x %>% distinct(breaks, .keep_all=TRUE))
sts_smooth <- sts %>% map2(brks, \(x, y) x %>% mutate(
  smooth_state = cut(
    position,
    y$breaks,
    labels=y$hmm_state[1:length(y$hmm_state)-1]))
) %>% bind_rows() %>% ungroup()
chr_ends <- chr_ends %>% mutate(
  start = case_when(endsWith(chromosome, 'R') ~ lag(length), .default = 0)
  )
bound_states <- sts_smooth %>% bind_rows()
mutate_chroms <- bound_states %>% left_join(chr_ends) %>%
  mutate(
    full_chrom = case_match(
      chromosome,
      'chr2L'  ~ 'chromosome 2',
      'chr2R' ~ 'chromosome 2',
      'chr3L' ~ 'chromosome 3',
      'chr3R' ~ 'chromosome 3',
      'chrX' ~ 'chromosome X'
      )
    ) %>%
  mutate(position = position / 1000000) %>% filter(smooth_state != "transition")

plot1 <- mutate_chroms %>% mutate(
    mod_pos = case_when(
      endsWith(chromosome, 'R') ~ position + length, .default = position
      )
    )
plot2 <- bound_states %>% left_join(chr_ends) %>%
  filter(startsWith(chromosome, 'chr2'))
intvls_joined <- intvls %>%
  bind_rows() %>%
  mutate(begin = start/1000000) %>%
  select(!start) %>%
  rename(lgth_intvl = length) %>%
  mutate(stop = stop/1000000) %>%
  mutate(lgth_intvl = lgth_intvl/1000000) %>%
  left_join(chr_ends) %>%
  mutate(length=length/1000000) %>%
  select(!start) %>%
  rename(start = begin) %>%
  filter(hmm_state != 'transition') %>%
  mutate(
    full_chrom = case_match(
      chromosome,
      'chr2L'  ~ 'chromosome 2',
      'chr2R' ~ 'chromosome 2',
      'chr3L' ~ 'chromosome 3',
      'chr3R' ~ 'chromosome 3',
      'chrX' ~ 'chromosome X'
      )
    ) %>% filter(startsWith(chromosome, 'chr2'))
intvls_joined$hmm_state <- factor(
  intvls_joined$hmm_state, levels=c('oregonr', 'w1118', 'het')
  )
# intvls_joined$chromosome <- factor(
#   intvls_joined$chromosome, levels=c(
#     ' ',
#     'chrX',
#     'chr2L', 'chr2R', 'chr3L', 'chr3R')
#   )
plot1 %>% write_csv(snakemake@output[['plot1']])
plot2 %>% write_csv(snakemake@output[['plot2']])
intvls_joined %>% write_csv(snakemake@output[['invls']])

# library(tidyverse)
# library(ggthemes)

# hmm_intervals <- read_csv(snakemake@input[['intervals']])
# hmm_states <- read_csv(snakemake@input[['states']])
# interval_states <- list(hmm_intervals, hmm_states) %>% map(
#   \(x) x %>% mutate(
#     chromosome = case_match(
#       chromosome,
#       1 ~ 'chr2L',
#       2 ~ 'chr2R',
#       3 ~ 'chr3L',
#       4 ~ 'chr3R',
#       5 ~ 'chrX'
#       )
#     ) %>% group_by(chromosome))
# # split on chromosome because we need to apply different 
# # breakpoints to each chromosome
# intvls <- interval_states[[1]] %>% select(!chrom_id) %>% group_split()
# sts <- interval_states[[2]] %>% group_split()
# # chromosomes aren't listed in the intervals df if they didn't have a CO 
# # so we have to grab the hmm state from those to plot them
# chromosomes <- tibble(
#   'chromosome' = list(snakemake@params[['chromosomes']]) 
#   )
# distinct_chrs <- interval_states[[1]] %>% distinct(chromosome) %>% 
#   arrange(chromosome)
# missing_chrs <- chromosomes %>% 
#   filter(! chromosome %in% distinct_chrs$chromosome)
# smoothed_st_values <- interval_states[[2]] %>%
#   filter(chromosome %in% missing_chrs$chromosome) %>%
#   select(sample, chromosome, hmm_state1) %>%
#   group_split() %>% map(\(x) head(x, n=1))

# chr_ends <- read_tsv(
#   snakemake@input[['chrom_lengths']], col_names=c('chromosome', 'length')
#   )
# chr_ends <- chr_ends %>% mutate(
#   chromosome = case_match(
#     chromosome,
#     1 ~ 'chr2L',
#     2 ~ 'chr2R',
#     3 ~ 'chr3L',
#     4 ~ 'chr3R',
#     5 ~ 'chrX'
#     )
#   )

# # get the chromosome ends so we can make a list of intervals
# smoothed_st_values <- smoothed_st_values %>%
#   map(
#     \(x) x %>% left_join(chr_ends) %>% mutate(start = 1) %>% 
#       mutate(stop = length) %>% rename(hmm_state = hmm_state1)
#     )

# st_joined <- smoothed_st_values %>% bind_rows()
# intvls <- interval_states[[1]] %>% bind_rows() %>%
#   bind_rows(smoothed_st_values) %>% group_by(chromosome) %>% group_split()
# # create a breaks column to show all the breakpoints in a single column as opposed to the start and stop columns
# # should be all of the stop values with 1 as the start
# brks <- map(
#   intvls, \(x) c(1, x[['stop']]) %>% as_tibble_col(column_name='breaks')
#   )
# # need to add an empty state at the end because we added an extra row
# # this state will be ignored anyway because it's at the end of the chromosome
# brks <- brks %>% map2(
#   intvls, \(x, y) x %>% mutate(hmm_state = c(y[['hmm_state']], NA))
#   )

# # labels should be 1 row shorter than breaks, because the labels are assigned to intervals between breaks
# # last label is the one that should be excluded
# sts[[1]] %>% write_csv('sts_chr2L.csv')
# brks[[1]] %>% write_csv('brks_chr2L.csv')
# brks[[1]]$
# brks[[1]]$hmm_state[1:length(brks[[1]]$hmm_state)-1]
# sts_smooth <- sts %>% map2(brks, \(x, y) x %>% mutate(
#   smooth_state = cut(
#     position,
#     y$breaks,
#     labels=y$hmm_state[1:length(y$hmm_state)-1]
#     )
#   )
# ) %>% bind_rows() %>% ungroup()
# chr_ends <- chr_ends %>% mutate(
#   start = case_when(endsWith(chromosome, 'R') ~ lag(length), .default = 0)
#   )
# bound_states <- sts_smooth %>% bind_rows()
# mutate_chroms <- bound_states %>% left_join(chr_ends) %>%
#   mutate(
#     full_chrom = case_match(
#       chromosome,
#       'chr2L'  ~ 'chromosome 2',
#       'chr2R' ~ 'chromosome 2',
#       'chr3L' ~ 'chromosome 3',
#       'chr3R' ~ 'chromosome 3',
#       'chrX' ~ 'chromosome X'
#       )
#     ) %>% 
#   mutate(position = position / 1000000) %>% filter(smooth_state != "transition")

# plot1 <- mutate_chroms %>% mutate(
#     mod_pos = case_when(
#       endsWith(chromosome, 'R') ~ position + length, .default = position
#       )
#     )
# plot2 <- bound_states %>% left_join(chr_ends) %>% 
#   filter(startsWith(chromosome, 'chr2'))
# intvls_joined <- intvls %>%
#   bind_rows() %>%
#   mutate(begin = start/1000000) %>%
#   select(!start) %>%
#   rename(lgth_intvl = length) %>%
#   mutate(stop = stop/1000000) %>%
#   mutate(lgth_intvl = lgth_intvl/1000000) %>%
#   left_join(chr_ends) %>%
#   mutate(length=length/1000000) %>%
#   select(!start) %>%
#   rename(start = begin) %>%
#   filter(hmm_state != 'transition') %>% 
#   mutate(
#     full_chrom = case_match(
#       chromosome,
#       'chr2L'  ~ 'chromosome 2',
#       'chr2R' ~ 'chromosome 2',
#       'chr3L' ~ 'chromosome 3',
#       'chr3R' ~ 'chromosome 3',
#       'chrX' ~ 'chromosome X'
#       )
#     ) %>% filter(startsWith(chromosome, 'chr2'))
# intvls_joined$hmm_state <- factor(
#   intvls_joined$hmm_state, levels=c('oregonr', 'w1118', 'het')
#   )
# # intvls_joined$chromosome <- factor(
# #   intvls_joined$chromosome, levels=c(
# #     ' ',
# #     'chrX',
# #     'chr2L', 'chr2R', 'chr3L', 'chr3R')
# #   )
# plot1 %>% write_csv(snakemake@output[['plot1']])
# plot2 %>% write_csv(snakemake@output[['plot2']])
# intvls_joined %>% write_csv(snakemake@output[['intvls']])
