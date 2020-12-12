#!/user/bin/env Rscript

## rreggiar@ucsc.edu

## generate 'colData' for tximeta aggregation of data using the naming conventions established 
##+ for organizing RNAseq data

# data/
## experiment_name/
### sample_name -- condition.rep.input_vol* *(if needed)

#[tximeta vignette](https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html#Analysis_starts_with_sample_table)
#[alevin process aware](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

suppressPackageStartupMessages({
library(tidyverse)
library(tximeta)
library(SummarizedExperiment)
library(HDF5Array)
})

parse.input <- function(data.path, output.name){
    ## assembles metadata object from input data

    files <- dir(file.path(data.path))

    # names(files) <- files

    # print(files)

    output.files <- list.files(data.path, full.names = TRUE)
    names(output.files) <- output.files
    output.files.samples <- list.files(output.files, full.names = TRUE)
    names(output.files.samples) <- output.files.samples
    # print(output.files.samples)


    # creates empty df to fill with sample info
    sample.info <- data.frame('files' = file.path(),
    						  'names' = character(),
                              'condition' = character(),
                              'rep' = numeric(), 
                              'input_vol' = numeric())

    # iteratively add to the empyty df the parsed sample names
    # here with format 'condition.rep'
    lapply(output.files.samples, function(x){

    		x.path <- Sys.glob(file.path(x,paste0("*v35.",output.name,"*"), "quant.sf"))[1]
    		x <- basename(x)
            sample.split <- str_split(x, '[.]', n=3)[[1]]
            condition <- sample.split[1]
            rep <- sample.split[2]
            input.vol <- sample.split[3]
            temp.df <- data.frame('files' = x.path,
            					  'names' = x,
                                  'condition' = condition,
                                  'rep' = rep,
                                  'input_vol' = input.vol)
            sample.info <<- rbind(sample.info, temp.df)
    })

    # convert to data frame
    # merge in the rest of the metadata
    duplicates <- c('panc.10.2.5', 'panc.6.2.5', 'panc.7.2.5', 'panc.9.2.5')#,
    #                 'ctrl.3.2.5', 'ctrl.2.3.4', 'ctrl.1.3.0', 'ctrl.1.5.0',
    #                 'ctrl.1.2.0', 'ctrl.1.1.5' , 'ctrl.1.0.5')

    head(sample.info)

    sample.info.df <<- 
        sample.info %>%
        bind_rows() %>%
        filter(! names %in% duplicates,
    		   condition != 'luad') %>%
    		   # ! grepl('intra', input_vol)) %>%
        filter(! is.na(files)) %>%
        mutate(condition = relevel(as.factor(condition), ref = 'ctrl'))

    print(sample.info.df)
}


paths.in <- scan(file=file("stdin", "r"), what="character", n=4)
data.path <- paths.in[1]
txome.path <- paths.in[2]
output.name <- paths.in[3]
txome.tsv <- paths.in[4]

print('parse input dir')
parse.input(data.path, output.name)

print('load tximeta')
tximeta::loadLinkedTxome(txome.path)

tx.to.gene <- read_csv('/public/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv', col_names=F)


if (output.name == 'ucsc.rmsk.salmon'){
	metaSkip = TRUE
	# print('import data with tximeta')
	# txi <- tximeta::tximeta(coldata = sample.info.df, type = 'salmon', skipMeta=metaSkip)

	outpath <- file.path("/public/groups/kimlab/exoRNA-biomarkers-panc/output.data/")

	# print('save tx h5')
	# saveHDF5SummarizedExperiment(txi, dir=paste0(outpath,paste0(output.name, '_h5_se')) , replace=TRUE)

	gxi <- tximeta::tximeta(coldata = sample.info.df, 
		type = 'salmon', 
		skipMeta=metaSkip, 
		txOut=FALSE, 
		tx2gene=tx.to.gene)

	print('save gene h5')
	saveHDF5SummarizedExperiment(gxi, dir=paste0(outpath,paste0(output.name, '_gene_h5_se')), replace=TRUE)


} else{
	metaSkip = FALSE
	print('import data with tximeta')
	txi <- tximeta::tximeta(coldata = sample.info.df, type = 'salmon', skipMeta=metaSkip)

	outpath <- file.path("/public/groups/kimlab/exoRNA-biomarkers-panc/output.data/")

	print('save tx h5')
	saveHDF5SummarizedExperiment(txi, dir=paste0(outpath,paste0(output.name, '_h5_se')) , replace=TRUE)

	gxi <- summarizeToGene(txi)

	print('save gene h5')
	saveHDF5SummarizedExperiment(gxi, dir=paste0(outpath,paste0(output.name, '_gene_h5_se')), replace=TRUE)

	if (!is.null(txome.tsv)){
		cg <- read.delim(txome.tsv, header = TRUE, as.is = TRUE)

		colnames(cg)[colnames(cg) == 'intron'] <- 'unspliced'
		split.gxi <- tximeta::splitSE(gxi, cg, assayName = 'counts')

		saveHDF5SummarizedExperiment(split.gxi, 
			dir=paste0(outpath,paste0(output.name, '_gene_split_h5_se')), 
			replace=TRUE)
	}

}




