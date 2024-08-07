suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("glmGamPoi"))

#' @param sc_path Path to scRNA counts file (standard CytoSPACE input file format)
#' @param ct_path Path to cell type labels file (standard CytoSPACE input file format)
#' @param st_path Path to ST data (standard CytoSPACE input file format)
#' @param outdir Path to output directory
#' @param prefix Prefix of output file
#' @param disable_downsampling If TRUE, disables downsampling of scRNA-seq dataset for SCTransform() regardless of dataset size
#' @import data.table dplyr Seurat
get_cellfracs_seuratv3 <- function(sc_path, ct_path, st_path, outdir, prefix, disable_downsampling=FALSE){
    message(Sys.time(), " Load ST data")
    if(endsWith(st_path,'.csv')){
        st_delim <- ','
    } else{
        st_delim <- '\t'
    }
    st <- fread(st_path, sep = st_delim, header = TRUE, data.table = FALSE)
    rownames(st) = st[,1]; st = st[,-1]
    st[is.na(st)] <- 0
    st <- as.matrix(st) 
    st <- CreateSeuratObject(st) %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() # edit by Sahil Sahni
    message(Sys.time(), " Load scRNA data")
    if(endsWith(sc_path,'.csv')){
        sc_delim <- ','
    } else{
        sc_delim <- '\t'
    }
    scrna <- fread(sc_path, sep = sc_delim, header = TRUE, data.table = FALSE)
    rownames(scrna) = scrna[,1]; scrna = scrna[,-1]
    scrna[is.na(scrna)] <- 0
    num_cells <- dim(scrna)[2]
    scrna <- as.matrix(scrna)

    if(endsWith(ct_path,'.csv')){
        ct_delim <- ','
    } else{
        ct_delim <- '\t'
    }

    celltypes <- fread(ct_path, sep = ct_delim, header = TRUE, data.table = FALSE)
    rownames(celltypes) = celltypes[,1]
    colnames(celltypes)[2] <- "CellType"
    celltypes <- celltypes[colnames(scrna),]

    scrna <- CreateSeuratObject(counts = scrna, meta.data = celltypes) 
    Idents(scrna) <- scrna$CellType
    cell_types <- unique(scrna$CellType)
    num_cell_types <- length(cell_types)
    names(cell_types) <- make.names(cell_types)
    Idents(scrna) <- make.names(Idents(scrna))
    if((num_cells > 10000) && (!disable_downsampling)) {
        warning("Downsampling to 10000 cells with equal number of cells for each cell type, for the purposes of cell type fraction estimation only. ",
                "This will not affect the dataset in the main CytoSPACE run. ",
                "To disable downsampling, please run get_cellfracs_seuratv3.R separately with the --disable-fraction-downsampling flag specified. ",
                "The output file [prefix]Seurat_weights.txt can then be passed to CytoSPACE with the -ctfep flag.")
        scrna <- subset(scrna, downsample = round(10000/num_cell_types))
    }
    if((dim(scrna)[2] < 300)) {
        warning("Please note that there may be an error in SCTransform() if there are too few cells available in the scRNA-seq dataset.")
    }
    scrna <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() # edit by Sahil Sahni
    
    cell_index_vec <- Idents(scrna)
    message(Sys.time(), " Integration ")
    anchors <- FindTransferAnchors(reference = scrna, query = st, normalization.method = "LogNormalize")
    predictions.assay <- TransferData(anchorset = anchors, refdata = cell_index_vec, k.weight=46)
    colnames(predictions.assay) <- gsub('prediction.score.', '', colnames(predictions.assay))
    predictions.assay <- predictions.assay[, 2:(ncol(predictions.assay)-1), drop=FALSE]
    col_sums <- colSums(predictions.assay)
    cellfrac <- col_sums / sum(col_sums)
    cellfrac <- data.frame(Index = names(cellfrac), Fraction = cellfrac)
    cellfrac$Index <- cell_types[cellfrac$Index]

    write.table(predictions.assay, file.path(outdir, paste0(prefix, "Seurat_cellfracs.txt")), quote = FALSE, sep = "\t")
    write.table(t(cellfrac), file.path(outdir, paste0(prefix, "Seurat_weights.txt")), quote = FALSE, sep = "\t", col.names = FALSE)
} 


parser <- ArgumentParser()
# specify options 
parser$add_argument("--scrna-path", type="character", default="",
                    help="Path to single cell RNA-seq data, with rows as genes and columns as single cells")
parser$add_argument("--ct-path", type="character", default="",
                    help="Path to cell type annotation file, with first column showing single cell names, and second column showing cell type annotations")
parser$add_argument("--st-path", type="character", default="",
                    help="Path to Spatially-Resolved Transcriptomics (SRT) data, with rows as genes and columns as spots")
parser$add_argument("--outdir", type="character", default="./", help="Output directory")
parser$add_argument("--prefix", type="character", default="CytoSPACE_input.", help="Prefix for output files")
parser$add_argument("--disable-fraction-downsampling", action="store_true", help="Does not downsample scRNA-seq dataset if specified")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
get_cellfracs_seuratv3(sc_path = args$scrna_path, ct_path = args$ct_path, st_path = args$st_path,
                        outdir = args$outdir, prefix = args$prefix, disable_downsampling=args$disable_fraction_downsampling)
