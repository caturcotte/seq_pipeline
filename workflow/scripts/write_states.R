library(tidyverse)
library(parallel)
library(furrr)

nnodes <- snakemake@threads - 1
cl <- parallel::makeCluster(nnodes)
plan(cluster, workers = cl)
hmm_intervals <- read_csv(snakemake@input[['intervals']])
hmm_states <- read_csv(snakemake@input[['states']])
hmm_intervals <- hmm_intervals %>% mutate(
  chromosome = case_match(
    chromosome,
    1 ~ 'chr2L',
    2 ~ 'chr2R',
    3 ~ 'chr3L',
    4 ~ 'chr3R',
    5 ~ 'chrX'
  )
)

hmm_states <- hmm_states %>% mutate(
  chromosome = case_match(
    chromosome,
    1 ~ 'chr2L',
    2 ~ 'chr2R',
    3 ~ 'chr3L',
    4 ~ 'chr3R',
    5 ~ 'chrX'
  )
)
intervals_2 <- hmm_intervals %>% 
  select(sample, chromosome, start, stop, hmm_state) %>% 
  filter(chromosome == 'chr2L' | chromosome == 'chr2R') %>% 
  mutate(
    start = case_when(
      chromosome == 'chr2R' ~ start + 23513712, .default = start)
    ) %>% 
  mutate(
    stop = case_when(chromosome == 'chr2R' ~ stop + 23513712, .default = stop)
    ) %>% 
  select(!chromosome) %>% 
  group_by(sample)
states_2 <- hmm_states %>% 
  rename(sample = "Sample") %>% 
  select(sample, chromosome, position, base_geno, hmm_state1, hmm_state2) %>% 
  mutate(
    position = case_when(
      chromosome == 'chr2R' ~ position + 23513712, .default = position
      )
    ) %>% 
  select(!chromosome) %>% 
  group_by(sample)
states_2_list <- states_2 %>% group_split()
intervals_2_list <- intervals_2 %>% group_split()
starts <- intervals_2 %>% select(start, hmm_state) %>% group_split()
breaks_list <- 
  intervals_2 %>%
  ungroup() %>%
  ## one tibble with columns 'start' and 'end' per 'category':
  nest_by(sample) %>%
  ## dataframe into list of tibbles, named with category
  pull(data, sample) %>%
  ## tibbles into vector of breaks
  map(~ .x %>% as.matrix %>% c %>% unique %>% sort) 

## get tile by indexing into "break_list" via mapping: 
states_2_hmm <- states_2 %>% 
  mutate(
    tile = future_map2(position, sample,
                ~ cut(.x, breaks_list[[.y]]), .progress=TRUE
    ) %>% unlist
  )
states_2_hmm %>% write_csv(snakemake@outputs[[1]])
parallel::stopCluster(cl)
