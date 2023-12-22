#' Effect of Body Mass Index (BMI) on Coronary Artery Disease (CAD)
#'
#' It contains independent datasets from three genome-wide association studies (GWASs):
#' \enumerate{
#' \item Exposure dataset: A GWAS for BMI in round 2 of the UK BioBank (sample size: 336,107), \url{http://www.nealelab.is/uk-biobank}.
#' \item Outcome dataset: A GWAS A GWAS for CAD from the CARDIoGRAMplusC4D consortium  (sample size: ~185,000), with genotype imputation using the 1000 Genome Project, (PubMed 26343387).
#' \item Selection dataset: A GWAS for BMI in the Japanese population (sample size: 173,430), (PubMed 28892062).
#' }
#'
#' @docType data
#' @usage data(bmi.cad)
#'
#' @format A \code{data.frame} with 1119 rows and 42 variables.
#'
#' @keywords datasets
#'
#' @source \url{https://github.com/qingyuanzhao/mr.raps}
#'
"bmi.cad"


#' Effects of HDL subfractions on Coronary Artery Disease (CAD)
#'
#' It contains SNP-exposure and SNP-outcome association summary statistics from the following genome-wide association studies (GWASs):
#' \enumerate{
#' \item Exposure dataset: A GWAS for traditional lipids \url{http://csg.sph.umich.edu/willer/public/lipids2013/} (PubMed 24097068); A GWAS for subfractions \url{http://www.computationalmedicine.fi/data\#NMR_GWAS} (PubMed 27005778).
#' \item Outcome dataset: A GWAS for CAD from the CARDIoGRAMplusC4D consortium \url{http://www.cardiogramplusc4d.org/data-downloads/} (PubMed 28714975).
#' \item Selection dataset: A GWAS for traditional lipids \url{https://www.ebi.ac.uk/gwas/studies/GCST007141} (PubMed 29507422); A GWAS for subfractions \url{http://csg.sph.umich.edu/boehnke/public/metsim-2017-lipoproteins/} (PubMed 29084231).
#' }
#'
#' @docType data
#' @usage data(hdl_subfractions)
#'
#' @format A list with 3 elements with the first element being a data.frame with 24 columns (see below for column descriptions), the second element being an empty data.frame, the third element being an estimated correlation matrix.
#' \describe{
#'   \item{SNP}{SNP rsid}
#'   \item{effect_allele}{effect allele}
#'   \item{other_allele}{other allele}
#'   \item{gamma_exp1,...,gamma_exp9}{Estimated associations of each SNP with respectively traditional lipids HDL, LDL, TG, and HDL subfractions S-HDL-P, S-HDL-TG, M-HDL-P, M-HDL-C, L-HDL-P, and L-HDL-C.}
#'   \item{se_exp1,...,se_exp9}{Standard error estiamtes for gamma_exp1,...,gamma_exp9}
#'   \item{gamma_out1}{Estimated SNP-outcome association}
#'   \item{se_out1}{Standard error estimate for gamma_out1}
#'   \item{selection_pvals}{p-values for SNP-exposure associations in the selection dataset}
#' }
#' @keywords datasets
#'
#' @source \url{https://github.com/tye27/mr.divw}
#'
"hdl_subfractions"




