library(RTIGER)
setupJulia(JULIA_HOME="/work/users/c/a/cannecar/co_analysis/seq_pipeline/julia-1.0.5/bin")
sourceJulia()

filePath <- snakemake@input[[1]]
sampleName <- snakemake@params[['sample_name']]
generation <- snakemake@params[['generation']]
outPath <- dirname(filePath)

expDesign <- data.frame(files=filePath, name=sampleName)

# TODO: generalize these at some point
chromLengths <- c(
  23513712,
  25286936,
  28110227,
  32079331,
  23542271,
  3667352,
  )
names(chromLengths) <- c('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chrY')

if (generation == 'F1') {
  myres = RTIGER(
    expDesign = expDesign,
    outputdir = outPath,
    seqlengths = chromLengths,
    rigidity = 5,
    autotune=TRUE,
    nstates=2
    )
} else if (generation == 'G0') {
  myres = RTIGER(
    expDesign = expDesign,
    outputdir = outPath,
    seqlengths = chromLengths,
    autotune=TRUE,
    rigidity = 5
    )
}