#' Extracts a set of markers from PLINK, VCF, or BCF files
#' @description
#' The 'extractmarkers' function extracts a set of markers, either using their rsID or chromosome position, from PLINK, VCF, or BCF files. 
#' In using rsIDs for extraction, a text file containing one marker per line is accepted. When using GRCh38/GRCh37 positions, the text file should contain 4 columns:
#' the first column should contain the chromosome number only (e.g. "1" and not "chr1"); the second column should contain the starting bp position; 
#' the third column the final bp position; and the fourth column should specify the desired output name with no file extension (e.g. "snpsfromchr1").
#' @param input.file is either a VCF or BCF file present in the working directory. File path is also accepted, but not suggested.
#' @param plink.files (required) asks whether the input files are in PLINK format. Default is FALSE.
#' @param bed.file should be provided if plink.files = TRUE
#' @param bim.file should be provided if plink.files = TRUE
#' @param fam.file should be provided if plink.files = TRUE
#' @param markerID (required) asks whether the markers will be extracted by rsID. Set to FALSE if using positions.
#' @param snps.list (required) is a text file containing the list of markers to be extracted. 
#' If using rsIDs, each line should correspond to one marker. If using POS (positions), there should be 4 columns: [1] chromosome numner (integer), [2] starting base position, [3] final base position, and [4] output file name. 
#' @param output is the user's desired output name. This will only be applied when extracting markers via rsID. 
#' @param plink.exec (required) takes in the file path for the plink executable. Suggestion is to have plink.exe in the working directory.
#' @examples
#' extractmarkers(input.file = "mygenotype.vcf", plink.files = FALSE, markerID = TRUE, snps.list = "listsnps.txt", output = "success", plink.exec = "./plink.exe")
#' @examples extractmarkers(plink.files = TRUE, bed.file = "data.bed", bim.file = "data.bim", fam.file = "data.fam", markerID = FALSE, snps.list = "listsnps.txt", plink.exec = "./softwares/plink.exe")
#' @export

extractmarkers <- function(
    input.file = input, 
    plink.files = FALSE,
    bed.file = bed,
    bim.file = bim,
    fam.file = fam,
    markerID = FALSE,
    snps.list = snps.file, 
    output = output.name,
    plink.exec = plink.path){
  
  plink_exec <- function(PLINKoptions = "") system(paste(plink.exec,PLINKoptions)) #specify path to plink
  
  if(plink.files == TRUE){
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file) # read in plink files
    
    if(markerID == TRUE){
      plink_exec(stringr::str_c("--bfile ", input, " --const-fid 0 --cow --extract ", snps.list, " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", output))
    } else if(markerID == FALSE){
      
      positions <- read.csv(snps.list, sep="")
      pos_list <- split(positions, 1:nrow(positions))
      
      for (i in pos_list){ # require file to be [1] chr number, [2] starting bp, [3] last bp, [4] output names
        plink_exec(stringr::str_c(
          "--bfile ", 
          input, 
          " --cow --chr ", i[,1], " --from-bp ", i[,2], 
          " --to-bp ", i[,3],
          " --recode vcf --keep-allele-order --out ", i[4])) }
    } else{
      stop("markerID required. Specify TRUE if using SNPs ID for extraction, FALSE if using POS.")
    }
    
  } else if(plink.files == FALSE){
    
    if(markerID == TRUE){
      
      if(tools::file_ext(snps.list) == "txt"){
        # for vcf input files
        if(tools::file_ext(input.file) == "vcf"){
          plink_exec(stringr::str_c(
            "--vcf ", input.file, 
            " --cow --extract ", snps.list, 
            " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", output))
          
          # for bcf input files
        } else if(tools::file_ext(input.file) == "bcf"){
          plink_exec(stringr::str_c(
            "--bcf ", input.file, 
            " --const-fid 0 --cow --extract ", snps.list, 
            " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", output))
        }
      } else {
        stop("SNP file should be a txt file with markers on each line.")
      }
      
    } else if(markerID == FALSE){
      
      if(tools::file_ext(snps.list) == "txt"){
        
        positions <- read.csv(snps.list, sep="")
        pos_list <- split(positions, 1:nrow(positions))
        
        # for vcf input files
        if(tools::file_ext(input.file) == "vcf"){
          
          for (i in pos_list){ 
            plink_exec(stringr::str_c(
              "--vcf ", 
              input.file, 
              " --cow --chr ", i[,1], " --from-bp ", i[,2], 
              " --to-bp ", i[,3],
              " --recode vcf --keep-allele-order --out ", i[4]))
          }
          
          # for bcf input files
        } else if(tools::file_ext(input.file) == "bcf"){
       
             for (i in pos_list){ 
            plink_exec(stringr::str_c(
              "--bcf ", 
              input.file, 
              " --cow --chr ", i[,1], " --from-bp ", i[,2], 
              " --to-bp ", i[,3],
              " --recode vcf --keep-allele-order --out ", i[4]))

        }
      } 
      } else {
        stop("SNP position file should be a txt file indicating columns [1] the chr number, 
             [2] the starting bp pos, [3] the last bp pos, [4] output name for the file.")
      }
    } else {
    stop("markerID required. Specify TRUE if using SNPs ID for extraction, FALSE if using POS.")
  }
  }
}