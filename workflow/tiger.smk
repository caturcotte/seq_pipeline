rule create_tiger_inputs:
input:
    "{sample}_str"
def get_TIGER_inputs_and_mkdirs(new_tiger_directory: str):

    os.mkdir(new_tiger_directory)

    # Create an empty list to store the filenames
    file_list = []

    # Use the os module to get a list of all the files in the current directory

    for filename in os.listdir(os.getcwd()):

        if os.path.isfile(os.path.join(os.getcwd(), filename)):

            if filename.split('.')[-2] == "tiger_input":

                os.mkdir(new_tiger_directory+'/'+filename.split('.')[0])

                os.replace( os.getcwd()+'/'+filename , new_tiger_directory+'/'+filename.split('.')[0]+'/'+filename )

def create_TIGER_inputs_and_directory(master_dataframe, new_project_dir: str):

    if new_project_dir not in os.listdir():

        for sample_type in master_dataframe.sample_type.unique():

            for sample_num in master_dataframe.sample_num.unique(): #subset master by individual samples...
    
                filename = str(sample_type) + "_" + str(sample_num)+".tiger_input.txt"
    
                df = master_dataframe[(master_dataframe['sample_type']==sample_type) & (master_dataframe['sample_num']==sample_num)]
    
    
                df['tiger_chrom'] = df.apply(lambda row: tiger_chrom_name_dict[row['chromosome']], axis=1)
    
                df.to_csv(filename, sep="\t", header=False, index=False, columns = ['tiger_chrom', 'position', 'reference', 'ref_reads', 'variant', 'variant_reads'])
    
            get_TIGER_inputs_and_mkdirs(new_tiger_directory=new_project_dir)
    
        else:
            print("No TIGER input or folders made. This project folder already exists...")

def run_TIGER_pipeline(master_dataframe, project_dir: str):

    #try to make TIGER inputs


    if project_dir not in os.listdir(): # if project folder not made...

        print('making TIGER inputs...')
        create_TIGER_inputs_and_directory(master_dataframe, new_project_dir=project_dir) #make inputs, project directory, and sample sub-directories

        for sample_dir in os.listdir(project_dir):

            sample_folder = os.getcwd()+'/'+project_dir+'/'+sample_dir+'/'

            if sample_dir != '.DS_Store':

                if sample_dir+'.tiger_input.txt' in os.listdir(sample_folder) and sample_dir+'.CO_estimates.txt' not in os.listdir(sample_folder): #if folder has input but missing this TIGER output, do TIGER

                    subprocess.run(['sh', 'run_TIGER.sh', os.getcwd()+'/'+project_dir+'/'+sample_dir+'/', sample_dir])

                elif sample_dir+'.tiger_input.txt' in os.listdir(sample_folder) and sample_dir+'.CO_estimates.txt' in os.listdir(sample_folder): #if folder has input and has this TIGER output, skip...
                    print(sample_dir, 'already processed...')
                    continue

                elif sample_dir+'.tiger_input.txt' not in os.listdir(sample_folder): #warn if input not in folder for some reason
                    print(sample_dir, 'needs TIGER input...')



    elif project_dir in os.listdir(): #if project directory alredy exists
        for sample_dir in os.listdir(project_dir): #for every sample directory, run tiger

            if sample_dir != '.DS_Store': #skip mac's DS_store file...useless

                sample_folder = os.getcwd()+'/'+project_dir+'/'+sample_dir+'/'

                if sample_dir+'.tiger_input.txt' in os.listdir(sample_folder) and sample_dir+'.CO_estimates.txt' not in os.listdir(sample_folder): #if folder has input but missing this TIGER output, do TIGER

                    subprocess.run(['sh', 'run_TIGER.sh', os.getcwd()+'/'+project_dir+'/'+sample_dir+'/', sample_dir])

                elif sample_dir+'.tiger_input.txt' in os.listdir(sample_folder) and sample_dir+'.CO_estimates.txt' in os.listdir(sample_folder): #if folder has input and has this TIGER output, skip...
                    print(sample_dir, 'already processed...')
                    continue

                elif sample_dir+'.tiger_input.txt' not in os.listdir(sample_folder): #warn if input not in folder for some reason
                    print(sample_dir, 'needs TIGER input...')
                    