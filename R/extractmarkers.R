#' Extracts a set of markers from PLINK, VCF, or BCF files
#' @param input.file is either a VCF or BCF file present in the working directory. File path is also accepted, but not suggested.
#' @param plink.files (required) asks whether the input files are in PLINK format. Default is FALSE.
#' @param bed.file should be provided if plink.files = TRUE
#' @param bim.file should be provided if plink.files = TRUE
#' @param fam.file should be provided if plink.files = TRUE
#' @param snps (required) is a text file containing the list of markers to be extracted. Each line should correspond to one marker.
#' @param output (required) is the user's desired output name.
#' @param plink.exec (required) takes in the file path for the plink executable. Suggestion is to have plink.exe in the working directory.
#' @examples
#' extractmarkers(input.file = "mygenotype.vcf", plink.files = FALSE, snps = "listsnps.txt", output = "success", plink.exec = "./plink.exe")
#' @examples extractmarkers(plink.files = TRUE, bed.file = "data.bed", bim.file = "data.bim", fam.file = "data.fam", snps = "listsnps.txt", output = "success", plink.exec = "./softwares/plink.exe")
#' @export

extractmarkers <- function(
      input.file = input, 
      plink.files = FALSE,
      bed.file = bed,
      bim.file = bim,
      fam.file = fam,
      snps = snps.file, 
      output = output.name,
      plink.exec = plink.path){
   
   plink_exec <- function(PLINKoptions = "") system(paste(plink.exec,PLINKoptions))
   
   # for plink input files
   if(plink.files == TRUE){
      input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
      plink_exec(stringr::str_c("--bfile ", input, " --const-fid 0 --cow --extract ", snps, " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", output))
      } else if(plink.files == FALSE){
         if(tools::file_ext(snps) == "txt"){
            # for vcf input files
            if(tools::file_ext(input.file) == "vcf"){
               plink_exec(stringr::str_c("--vcf ", input.file, " --cow --extract ", snps, " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", output))
               
               # for bcf input files
            } else if(tools::file_ext(input.file) == "bcf"){
               plink_exec(stringr::str_c("--bcf ", input.file, " --const-fid 0 --cow --extract ", snps, " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", output))
            }
         } else {
            stop("SNP file should be a txt file with markers on each line.")
         }
   }
}
