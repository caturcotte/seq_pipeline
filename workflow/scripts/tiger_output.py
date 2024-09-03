import os
os.environ['MODIN_ENGINE'] = "unidist"
os.environ['UNIDIST_BACKEND'] = "mpi"
import modin.pandas as pd
import numpy as np


# def create_TIGER_master_df(project_dir=str, file_pattern=str):
def create_TIGER_master_df(files_list: list, ref_parental_genotype, alt_parental_genotype):
    #this dataframe this function makes is relatively large because it gives you genotype and HMM state calls for every SNP marker across each chromosome for each sample
    
    #This is closest to the 'raw' data that TIGER outputs, but this format is mostly useful for plotting and visualization of individual chromosomes for manual inspection of SNPs and crossover calls

    #get list of file names for all TIGER outputs
    # files_list = get_TIGER_files(project_dir=project_dir, file_pattern=file_pattern)

    data = []

    for file in files_list:

        #make dataframe; adjust column names if you end up using different TIGER output files than I did (file_pattern='.CO_estimates.txt')
        df = pd.read_csv(
            file,
            sep="\t",
            header=None,
            names = [
                "Sample",
                "chromosome",
                "position",
                "base_geno",
                "hmm_state1",
                "hmm_state2",
                "reference",
                "ref_reads",
                "variant",
                "var_reads"
                ]
            )

        #This dict will convert parental arabidopsis genotype calls (Col and Ler) to C. elegans genotype calls (N2 is Bristol WT, CB4856 is Hawaiian WT)
        #Replace these values with your preferred genotype labels as needed
        genotype_dict = {
            'CC':ref_parental_genotype,
            "CL":"het",
            "LL":alt_parental_genotype,
            "CU":"u"+ref_parental_genotype,
            "LU":'u'+alt_parental_genotype,
            "UN":"unknown",
            '?':"?"
            }
        df = df.replace(genotype_dict)
        
        #for our specific cross scheme to ID crossovers, heterozygous and homozygous CB4856 calls are should all be considered CB4856
        df['hmm_state1'] = df['hmm_state1'].replace({'het':alt_parental_genotype})
        
        #add a unique identifier for each chromosome
        #example... BSP-OR-001-1 = (project)-(gamete)-(sample number)-(chromosome)
        df['chrom_id'] = df['Sample'] + '-' + str(df['chromosome'])
        # (lambda row: row['Sample']+"-"+str(row['chromosome']), axis=1)

        print(file, 'done')
        data.append(df)

    # make dataframe from all samples in data list
    data = pd.concat(data)
    print(data)

    return data

def make_transition_intervals_df(master_df):
    
    #These new transition intervals are the potential crossover intervals that can be classified later as crossover or 'not crossover'

    new_dfs = []

    for chrom_id in master_df.chrom_id.unique():

        df = master_df[master_df['chrom_id']==chrom_id]

        starts = list(df.start)
        stops = list(df.stop)
        states = list(df.hmm_state)

        new_intervals = []

        for i in range(len(stops)-1):

            if states[i] != states[i+1]:

                row = [ df['sample'].unique()[0], df['chromosome'].unique()[0], stops[i], starts[i+1], 'transition', df['chrom_id'].unique()[0], (starts[i+1]-stops[i]) ]

                new_intervals.append(row)

        new_intervals_df = pd.DataFrame(new_intervals, columns=df.columns)


        full_intervals_df = pd.concat([df, new_intervals_df]).sort_values(by='start').reset_index(drop=True)

        new_dfs.append(full_intervals_df)


    print(new_dfs)
    return pd.concat(new_dfs)




def create_state_intervals_df(files_list: list, ref_parental_genotype, alt_parental_genotype):
    #This will create a much smaller dataframe than the marker dataframe above by collapsing runs of the same HMM state into intervals and creating a new row for a 2-marker transition interval between HMM states

    # files_list = get_TIGER_files(project_dir=project_dir, file_pattern=file_pattern)

    genotype_dict = {'CC':ref_parental_genotype, "CL":alt_parental_genotype, "LL":alt_parental_genotype, "CU":"u"+ref_parental_genotype, "LU":'u'+alt_parental_genotype, "UN":"unknown", '?':"?"}
    data = []

    for file in files_list:

        #create df given columns in files using file_pattern='.CO_estimates.breaks.txt' in TIGER outputs
        df = pd.read_csv(
            file,
            sep="\t",
            header=None,
            names = [
                "sample",
                "chromosome",
                "start",
                "stop",
                "hmm_state"
                ]
            )
        genotype_dict = {
            'CC':ref_parental_genotype,
            "CL":"het",
            "LL":alt_parental_genotype,
            "CU":"u"+ref_parental_genotype,
            "LU":'u'+alt_parental_genotype,
            "UN":"unknown",
            '?':"?"
            }
        print(df)
        df['hmm_state'] = df['hmm_state'].replace(genotype_dict)

        # df['hmm_state'] = df.apply(lambda row: genotype_dict[(row['hmm_state'])], axis=1)
        df['chrom_id'] = df['sample'] + '-' + str(df['chromosome'])
        df['length'] = df['stop'] - df['start'] + 1
#         df['gamete'] = df['sample'].unique()[0].split("-")[1]

        data.append(df)


    master_df = pd.concat(data)

    master_df = make_transition_intervals_df(master_df)

    master_df = master_df.sort_values(by=['chrom_id', 'chromosome', 'start']).reset_index(drop=True)

    return master_df


def get_marker_counts_per_interval(master_bases_df, master_intervals_df):
    
    #This is useful to retain marker density information for each interval, which is needed for later maths and classification of crossover intervals

    new_dfs = []

    for chrom_id in master_bases_df.chrom_id.unique():

        bases_df = master_bases_df[master_bases_df['chrom_id']==chrom_id]
        intervals_df = master_intervals_df[master_intervals_df['chrom_id']==chrom_id]

        starts = list(intervals_df.start)
        stops = list(intervals_df.stop)

        counts = []

        for i in range(len(starts)):

            counts.append(len(bases_df[(bases_df['position']>=starts[i]) & (bases_df['position']<=stops[i]) ]))

        intervals_df['marker_counts'] = counts

        new_dfs.append(intervals_df)

    markers_counted_df = pd.concat(new_dfs)

    return markers_counted_df

#dictionary of chromosome names and corresponding lengths, change these as needed!
# there were no decent SNPs on chr4 so it's omitted
dm6_fasta_names = list(snakemake.config['chromosomes'].keys())
dm6_chrom_lengths = list(snakemake.config['chromosomes'].values())
# dm6_chrom_lengths = [15114068, 15311845, 13819453, 17493838, 20953657, 17739129]

tiger_chrom_len_dict = dict(zip([1,2,3,4,5,6], dm6_chrom_lengths))
tiger_chrom_name_dict = dict(zip(dm6_fasta_names, [1,2,3,4,5,6]))

ref_parental_genotype = 'w1118' #change this to your reference genotype
alt_parental_genotype = 'oregonr' #change this to your alternate parental genotype

# sample_vcfs_path = '/Users/zac/Desktop/VCF_to_TIGER/sample_vcf/'

# cb_ref_vcf = '/Users/zac/Desktop/VCF_to_TIGER/reference_vcf/CB.aligned.to.N2.SNPs.noHet.no-repeats.vcf'

# sample_vcfs = get_files_in_directory(sample_vcfs_path)
# sample_vcfs

#generate TIGER marker and interval dataframes
tiger_marker_df = create_TIGER_master_df(
    snakemake.input['co_estimates'],
    ref_parental_genotype,
    alt_parental_genotype,
    # project_dir='TIGER_test_output',
    # file_pattern='.CO_estimates.txt'
    )
tiger_marker_df.to_parquet(
    snakemake.output['hmm_states']
    # 'TIGER_hmm_states.all_markers.pickle.gzip', compression='gzip'
    )
tiger_pre_intervals = create_state_intervals_df(
    snakemake.input['breaks'],
    ref_parental_genotype,
    alt_parental_genotype
    # project_dir='TIGER_test_output',
    # file_pattern='.CO_estimates.breaks.txt'
    )
tiger_intervals = get_marker_counts_per_interval(
    tiger_marker_df,
    tiger_pre_intervals
    )

tiger_intervals.to_csv(
    snakemake.output['hmm_intervals']
)
#     'TIGER_hmm_intervals.pickle.gzip',
#     compression='gzip'
#     )
# tiger_intervals