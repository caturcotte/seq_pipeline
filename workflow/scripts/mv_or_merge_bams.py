import subprocess

try:
    s = snakemake.input[1]
    subprocess.run(['sambamba', 'merge', '-t', str(snakemake.threads), *[str(i) for i in snakemake.output], *[str(i) for i in snakemake.input]])
except:
    subprocess.run(['mv', str(snakemake.input[0]), str(snakemake.output[0])])