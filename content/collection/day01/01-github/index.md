---
date: "2024-07-25"
draft: false
excerpt: You can share information about yourself with the community on GitHub by
  creating a profile README. GitHub shows your profile README at the top of your profile
  page.
links:
- icon: door-open
  icon_pack: fas
  name: website
  url: https://bakeoff.netlify.com/
- icon: github
  icon_pack: fab
  name: code
  url: https://github.com/apreshill/bakeoff
subtitle: Put your best foot forward, first.
title: A GitHub profile
weight: 1
---

# Tutorial: Summary data-based Mendelian Randomisation

Mendelian randomization (MR) has emerged as a powerful tool in epidemiology, enabling researchers to estimate causal relationships between modifiable risk factors and health outcomes using genetic variants. Often referred to as “Nature’s Randomised Control Trial”, MR is based on the principle that genetic variants associated with exposures are randomly distributed at conception. Consequently, using genetically-determined exposures can address confounding and reverse causation that plagues observational associations.

In particular, [summary-data Mendelian randomization (SMR)](https://yanglab.westlake.edu.cn/software/smr/#Overview) has gained traction due to its ability to utilize summary statistics from large genome-wide association studies (GWAS), allowing for robust causal inference without the need for individual-level data. SMR is also designed for efficient handling of large molecular datasets, for instance methylation and gene expression data. This makes it convenient to simultaneously analyse causal relationships between molecular risk factors and a single disease outcome, as well as between multiple molecular risk factors and various molecular outcomes. 

For this tutorial, I have attempted to create an end-to-end pipeline for SMR analysis, starting from the downloading of publicly available genetic association data all the way to visualisation of the results. This tutorial borrows from two excellent resources that have helped me grasp SMR. The gene-outcome pair and datasets used are as per [training material](https://cnsgenomics.com/data/teaching/GNGWS22/module2/SMR.pdf) from the 2022 edition of the Genetics and Genomics Winter School hosted by the University of Queensland. The format for the tutorial and code for preparing the GWAS data are from the [blog of a medical student](http://www.alearnerlin.top/index.php/创建您的网站与区块/博客/) from China which I enjoyed for its intuitive flow.

## Biological background

In this practical, we will be using SMR to assess the evidence for a causal relationship between *CACNA2D4* and schizophrenia. *CACNA2D4* encodes the auxiliary subunit α2δ-4 of voltage-gated calcium channels and belongs to the [*CACN* family of genes](https://pubmed.ncbi.nlm.nih.gov/22488967/) that have been associated with brain circuits, cognitive function and neuropsychiatric disorders. [Knockout mice for *CACNA2D4](https://pubmed.ncbi.nlm.nih.gov/35353835/)* exhibit significant behavioral impairments relevant to cognitive and motor functions, suggesting a role in neuropsychiatric conditions.

![theory-1.png](Tutorial%20Summary%20data-based%20Mendelian%20Randomisatio%2015ecdcf09b31807093d2cb638e8ec681/theory-1.png)

SMR can effectively investigate causal links and biological pathways between the *CACNA2D4* gene and schizophrenia by integrating data from genome-wide association studies (GWAS) and expression quantitative trait locus (eQTL) studies. GWAS identifies single nucleotide polymorphisms (SNPs) associated with schizophrenia, while eQTL studies reveal how these SNPs affect *CACNA2D4* expression across different genotypes. By leveraging genetic variants linked to both *CACNA2D4* and schizophrenia, researchers can assess whether variations in *CACNA2D4* expression might causally influence the risk of developing schizophrenia.

![theory-2.png](Tutorial%20Summary%20data-based%20Mendelian%20Randomisatio%2015ecdcf09b31807093d2cb638e8ec681/theory-2.png)

## **Step 1: Installing the SMR software**

To download and install the SMR software on a Mac, follow these steps:

1. **Download the Software**:
    - Click on the [download link](https://yanglab.westlake.edu.cn/software/smr/#Download) for the latest version of the software. For Mac users, select `smr_Mac_v1.03.zip`.
2. **Extract the Downloaded File**:
    - Once the download is complete, locate the `smr_Mac_v1.03.zip` file
    - Double-click on the zip file to extract its contents. This will create a folder containing the executable files necessary for running SMR.
3. **Move the Executable to /usr/local/bin**:
    - Open a Terminal window (Applications > Utilities > Terminal).
    - Move the executable to **`/usr/local/bin`** using the command (Replace `/path/to/extracted/smr` with the actual path of the executable file you extracted from the zip folder)
        
        ```bash
        sudo mv /path/to/extracted/smr /usr/local/bin
        ```
        
    
    This directory is standard for user-installed binaries, allowing you to run the SMR command from any location, rather than just from the directory containing the executable.
    
4. **Verify Installation**:
    - To ensure that SMR is installed correctly, type `smr-1.3.1-macos-arm64` in your Terminal and press Enter. If everything is set up properly, you should see a list of commands or options available for SMR.

## Step 2: Create a working directory

Create a working directory for your SMR analysis, for instance, **`smr_practical`**. The file paths in this tutorial will be based on the name **`smr_practical`**, so adjust accordingly if you choose a different name for your directory.

## Step 3: Download the required data for SMR

This step covers the preparation of data for [SMR.](http://SMR.In) Considering that we would like to investigate causal relationships between *CACNA2D4*  and schizophrenia, we would need three datasets:

- **eQTL data for *CACNA2D4*** (SNP-gene expression; SMR BESD format)
- **GWAS data for schizophrenia** (SNP-schizophrenia; SMR .ma format)
- **Genetic for the reference population** (plink binary format, i.e. `.bed` , `.bim` and `.bam` )

### Downloading the eQTL data

For this practical, we will use [blood eQTL data](https://yanglab.westlake.edu.cn/software/smr/#DataResource) (hg19; Westra et al. 2013, Nat Genet) that can be downloaded from the SMR website and has already been formatted for SMR analyses. 

Upon extracting the downloaded .zip file, you will find three files with the extensions `.besd`, `.epi` and `.esi.` These files represent different components of the eQTL summary statistics. Storing eQTL summary statistics in three separate files is essential for efficiently managing large molecular datasets. For detailed information on the contents and column names of each file, refer to the [Data Management](https://yanglab.westlake.edu.cn/software/smr/#DataManagement) section of the SMR website.

Since we are focusing on the *CACNA2D4* gene ([ENSG00000151062](https://www.ensembl.org/Homo_sapiens/geneview?gene=ENSG00000151062)) we should extract eQTL results for *CACNA2D4* from the eQTL dataset using the following command:

```r
smr-1.3.1-macos-arm64 --beqtl-summary westra_eqtl_hg19 --extract-probe myprobe.list  --make-besd --out  CACNA2D4
```

The **`myprobe.list`** file contains the probe ID for *CACNA2D4*, ILMN_1696317, which is listed in the second column of the .epi file. Running the command extracts the eQTL data for *CACNA2D4* in BESD format.

![eqtl_data.png](Tutorial%20Summary%20data-based%20Mendelian%20Randomisatio%2015ecdcf09b31807093d2cb638e8ec681/eqtl_data.png)

⚠️ If you receive an error message or an empty output (zero bytes), ensure you are running the code in the same directory as your eQTL files.

**If you're using your own eQTL data:** you'll need to generate the three required file types from the eQTL summary statistics. This is often necessary for tissue-specific eQTL data, e.g. from the GTEx database .To create these files, follow the steps in the [‘Make a BESD file’ section](https://yanglab.westlake.edu.cn/software/smr/#MakeaBESDfile). You’ll need to generate a `.esd` file for each exposure (e.g., gene) and an `eqtl.flist` file listing the paths to the .esd files (typically created using R). Running the command below will then generate the `.besd`, `.epi`, and `.esi` files in your working directory:

```r
smr --eqtl-flist my.flist --make-besd --out mybesd 
```

### Downloading the GWAS data

GWAS data for various traits is publicly available for download from sources like the [IEU OpenGWAS project](https://gwas.mrcieu.ac.uk/) and the [GWAS catalog.](https://www.ebi.ac.uk/gwas/) 

In this tutorial, we will download [GWAS data for schizophrenia](https://gwas.mrcieu.ac.uk/datasets/ieu-b-5102/) (Trubetskoy et al. 2022, Nature). from the IEU OpenGWAS project. Select **`Download VCF file.`**The download will take approximately 25 minutes.

Once downloaded, use the following code to extract the relevant columns and format them for SMR:

```r
setwd("/Users/konstanzetan/Desktop/smr_practical")
# load required libraries
library(VariantAnnotation)
library(gwasglue)

# read GWAS data for desired trait
gwasdata<-readVcf("ieu-b-5102.vcf.gz")
gwasdata<-gwasvcf_to_TwoSampleMR(gwasdata,"outcome")
gwasdata<-data.frame(gwasdata$SNP,
                 gwasdata$effect_allele.outcome,
                 gwasdata$other_allele.outcome,
                 gwasdata$eaf.outcome,
                 gwasdata$beta.outcome,
                 gwasdata$se.outcome,
                 gwasdata$pval.outcome,
                 gwasdata$samplesize.outcome)
colnames(gwasdata)<-c("snp","A1","A2","freq","b","se","p","n")
write.table(gwasdata,"gwas_schizophrenia.ma",quote = FALSE,row.names = FALSE)
```

### Downloading the genetic data for the reference population

We will download the genetic data from the [MAGMA software](https://cncr.nl/research/magma/) homepage, which provides genetic data in a population-specific format and in the PLINK binary format required by SMR. Ensure that the reference population matches the one used for obtaining GWAS and eQTL summary statistics, which in this case is European.

SMR uses genetic data for the reference population during quality control, to filter out SNPs that have abnormal allele frequencies, as well as to perform the [**HEterogeneity In Dependent Instruments (HEIDI)](https://mr-dictionary.mrcieu.ac.uk/term/heidi/)** test that assesses the probability of a single genetic variant underlying the putative causal (SMR-significant) association between a given exposure and trait. Briefly, this is based on an assumption that if there is a single causal variant, there should be homogeneity in the causal estimates calculated for all SNPs in LD with the SNP being used in SMR as the instrumental variable. Heterogeneity in causal estimates indicate the likelihood of more than one genetic variant contributing to the observed association.

⚠️ Genetic data can be downloaded from the 1000Genomes FTP site, but these files lack RSIDs. Therefore, you will need to manually map the RSIDs to the genetic data to ensure that the variant identifiers match across GWAS, QTL, and genetic data files, allowing the SMR analysis to run successfully.

## Step 4: Running the SMR

Here comes the easiest and sweetest part! Once all the files have been properly formatted, this should run smoothly

```r
smr-1.3.1-macos-arm64 --bfile ./g1000_eur --gwas-summary ./gwas_schizophrenia.ma --beqtl-summary ./westra_eqtl_hg19/CACNA2D4 --out CACNA2D4schizhophrenia --thread-num 10
```

⚠️  The ‘file could not be found/opened’ error message is a common one at this step, and typically results from directory mis-specification. In case it helps with troubleshooting, the directory structure for this SMR analysis is as follows: 

![smr_final_dir.png](Tutorial%20Summary%20data-based%20Mendelian%20Randomisatio%2015ecdcf09b31807093d2cb638e8ec681/smr_final_dir.png)

## Step 5: Visualising results

SMR results are typically visualised as [locus plots](https://yanglab.westlake.edu.cn/software/smr/#SMRlocusplot19), which show the overlap between genetic associations for the exposure and the outcome.  Code for generating such plots from SMR-formatted files is available on the SMR website under the section [SMR locus plot](https://yanglab.westlake.edu.cn/software/smr/#SMRlocusplot19). However, I prefer using the R package [LocusZoom](https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html) because of its customisation options. In this section, I share my code fro creating LocusZoom plots.

First, retrieve the eQTL data for the gene of interest. Since our eQTL data is in **`.besd`** format, we need to obtain the summary statistics in a human-readable format as input for LocusZoom, using the following command: 

```r
smr-1.3.1-macos-arm64 --beqtl-summary ./westra_eqtl_hg19/CACNA2D4 --query 1 --probe ILMN_1696317 --out CACNA2D4_eqtl 
```

The output file, `CACNA2D4_eqtl` , contains association statistics for 117 SNPs tested for their relationship with *CACNA2D4* in a human-readable format. Note that the p-value threshold is set to 1, such that all SNP-*CACNA2D4* tests in the specified region are included, not just significant results. 

With the GWAS and eQTL data in human-readable format, we can use the following R code to generate the locus plots:

```r
# load required libraries
library(dplyr)
library(locuszoomr)
library(tidyr)
library(EnsDb.Hsapiens.v75) # used for hg19, for hg38 onward use AnnotationHub

# define window size for plot as +/-500kb around the SMR instrumental variable
## cis-eQTL analyses are often restricted to 1Mb between SNP and eQTL; 500kb seems to be a nice window for viewing
window_start <- 1901485-500000 # coordinate of SMR IV -5000
window_end <- 1901485+500000 # coordinate of SMR IV + 5000

# read eqtl data
cacna2d4_eqtl <- read.delim("CACNA2D4_eqtl.txt") # read in eqtl results
cacna2d4_eqtl <- cacna2d4_eqtl %>%
	rename(snp = SNP) 

# create locus object for plot
eqtl_loc <- locus(data = cacna2d4_eqtl, chrom = 'Chr', labs = 'snp', pos = 'BP', p = 'p',
seqname = 12, xrange = c(window_start, window_end), ens_db = EnsDb.Hsapiens.v75) # change seqname to the chromosome where the gene is located
 
# read gwas data
schizophrenia_gwas <- read.delim("gwas_schizophrenia.ma", sep = " ")
# format gwas data
schizophrenia_gwas  <- schizophrenia_gwas %>%
  semi_join(cacna2d4_eqtl, by = c("snp")) %>%  # Keep only rows with matching SNPs in eqtl dataset 
  left_join(cacna2d4_eqtl %>% dplyr::select(snp, BP, Chr), by = c("snp"))  # Left join to get position data 

# create locus object for the plot
gwas_loc <- locus(data = schizophrenia_gwas, chrom = 'Chr', labs = 'snp', pos = 'BP', p = 'p',
seqname = 12, xrange = c(window_start, window_end), ens_db = EnsDb.Hsapiens.v75) # make sure the column names for the variables match between the association datasets! 

# label the instrumental variable on the plots
## here I edited the default scatterplot function to annotate SMR IV on the genetic association plots

pl <- quote({
  # add custom text label for instrumental variable
  lx <- loc$data$BP[loc$data$snp == "rs1044825"]
  ly <- loc$data$logP[loc$data$snp == "rs1044825"]
  #text(lx, ly, "16:30095236", pos = 4, cex = 0.8)
  # modify instrumental variable shape and color
  points(lx, ly, pch = 23, col = "purple", bg = "purple")
})

# define plot function 
scatter_plot_rmcol_indexsnp <- function (loc, index_snp = loc$index_snp, sentinel = loc$sentinel, pcutoff = 5e-08, chromCol = "royalblue", 
    sigCol = "red", cex = 1, cex.axis = 0.9, cex.lab = 1, xlab = NULL, 
    ylab = NULL, yzero = (loc$yvar == "logP"), xticks = TRUE, 
    border = FALSE, showLD = TRUE, LD_scheme = c("grey", "royalblue", 
        "cyan2", "green3", "orange", "red", "purple"), legend_pos = "topleft", 
    labels = NULL, label_x = 4, label_y = 4, add = FALSE, align = TRUE, 
    ...) 
{
    if (!inherits(loc, "locus")) 
        stop("Object of class 'locus' required")
    data <- loc$data
    if (is.null(xlab)) 
        xlab <- paste("Chromosome", loc$seqname, "(Mb)")
    if (is.null(ylab)) {
        ylab <- if (loc$yvar == "logP") 
            expression("-log"[10] ~ "P")
        else loc$yvar
    }
    hasLD <- "ld" %in% colnames(data)
    if (!"bg" %in% colnames(data)) {
        if (showLD & hasLD) {
            data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
            data$bg[is.na(data$bg)] <- 1L
            data$bg[data[, loc$labs] == index_snp] <- 7L
            data <- data[order(data$bg), ]
            LD_scheme <- rep_len(LD_scheme, 7)
            data$bg <- LD_scheme[data$bg]
        }
        else {
            data$bg <- chromCol
            if (loc$yvar == "logP") 
                data$bg[data[, loc$p] < pcutoff] <- sigCol
            data$bg[data[, loc$labs] == sentinel] <- "green"
        }
    }
    if (align) {
        op <- par(mar = c(ifelse(xticks, 3, 0.1), 3.5, 2, 1.5))
        on.exit(par(op))
    }
    if (loc$yvar == "logP" & !is.null(pcutoff)) {
        abl <- quote(abline(h = -log10(pcutoff), col = "darkgrey", 
            lty = 2))
    }
    else abl <- NULL
    pch <- rep(21L, nrow(data))
    if ("pch" %in% colnames(data)) 
        pch <- data$pch
    col <- "black"
    if ("col" %in% colnames(data)) 
        col <- data$col
    new.args <- list(...)
    if (add) {
        plot.args <- list(x = data[, loc$pos], y = data[, loc$yvar], 
            pch = pch, bg = data$bg, cex = cex)
        if (length(new.args)) 
            plot.args[names(new.args)] <- new.args
        return(do.call("points", plot.args))
    }
    ylim <- range(data[, loc$yvar], na.rm = TRUE)
    if (yzero) 
        ylim[1] <- min(c(0, ylim[1]))
    plot.args <- list(x = data[, loc$pos], y = data[, loc$yvar], 
        pch = pch, bg = data$bg, col = col, las = 1, font.main = 1, 
        cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, xlim = loc$xrange, 
        ylim = ylim, xlab = if (xticks) xlab else "", ylab = ylab, 
        bty = if (border) "o" else "l", xaxt = "n", tcl = -0.3, 
        mgp = c(1.7, 0.5, 0), panel.first = abl)
    if (length(new.args)) 
        plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    if (!is.null(labels)) {
        i <- grep("index", labels, ignore.case = TRUE)
        if (i) 
            labels[i] <- index_snp
        ind <- match(labels, data[, loc$labs])
        if (any(is.na(ind))) {
            message("label ", paste(labels[is.na(ind)], collapse = ", "), 
                " not found")
        }
        lx <- data[ind, loc$pos]
        ly <- data[ind, loc$yvar]
        labs <- data[ind, loc$labs]
        add_labels(lx, ly, labs, label_x, label_y, cex = cex.axis * 
            0.95)
    }
    if (xticks) {
        axis(1, at = axTicks(1), labels = axTicks(1)/1e+06, cex.axis = cex.axis, 
            mgp = c(1.7, 0.4, 0), tcl = -0.3)
    }
    else if (!border) {
        axis(1, at = axTicks(1), labels = FALSE, tcl = -0.3)
    }
    if (!is.null(legend_pos)) {
        if (showLD & hasLD) {
            legend(legend_pos, legend = c("0.8 - 1.0", "0.6 - 0.8", 
                "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"), title = expression({
                r^2
            }), y.intersp = 0.96, pch = 21, col = "black", pt.bg = rev(LD_scheme[-c(1, 
                7)]), bty = "n", cex = 0.8)
        }
    }
}
# generate locus plots
oldpar <- set_layers(2)
scatter_plot_rmcol_indexsnp(eqtl_loc, chromCol = "grey", col = NA, xticks = FALSE, pcutoff = 0.05, panel.last = pl) # adjust colors for significant points, shape of lead snp
scatter_plot_rmcol_indexsnp(gwas_loc, chromCol = "grey", col = NA, xticks = FALSE, pcutoff = 0.05, panel.last = pl) # adjust colors for significant points, shape of lead snp
genetracks(eqtl_loc, filter_gene_biotype = "protein_coding", cex.text = 0.5)
par(oldpar) 
```

The resulting plot visualizes genetic associations within a 1 Mb region centered on the SNP used as the instrumental variable (rs1044825) for the SMR analysis between *CACNA2D4* and schizophrenia. This SNP is highlighted in purple.

![plot_annnotated.png](Tutorial%20Summary%20data-based%20Mendelian%20Randomisatio%2015ecdcf09b31807093d2cb638e8ec681/plot_annnotated.png)

## Next steps

In practice, SMR is often supplemented with [genetic colocalization](https://chr1swallace.github.io/coloc/articles/a02_data.html) analysis to investigate the assumption that a single causal variant underlies the observed association. Colocalization analysis provides statistical evidence to determine whether the patterns displayed on the locus plot are statistically similar, beyond mere visual resemblance. I might do a follow-up post on colocalization analysis in the future. 🙂

## Further learning

- CNSGenomics provides an [excellent practical guide](https://cnsgenomics.com/data/teaching/GNGWS24/module6/Practicals/Module6_SMR_practical.pdf) on SMR, featuring questions that encourage critical thinking about the implications of SMR outputs and the design of follow-up experiments.

## Help me make this better

If you’ve tried this tutorial and encountered issues, or if you have ideas on how I can code things better, I’d be grateful if you could leave a comment below, or email me at [konstanz001@e.ntu.edu.sg](mailto:konstanz001@e.ntu.edu.sg)
