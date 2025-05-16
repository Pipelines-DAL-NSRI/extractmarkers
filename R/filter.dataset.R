#' @title Basic quality control of variants.
#' @description
#' This function performs QC of the variant dataset based on commonly used filters.
#' @param input.file is either a VCF or BCF file present in the working directory. File path is also accepted, but not suggested.
#' @param plink.files (required) asks whether the input files are in PLINK format. Default is FALSE.
#' @param bed.file should be provided if plink.files = TRUE
#' @param bim.file should be provided if plink.files = TRUE
#' @param fam.file should be provided if plink.files = TRUE
#' @param remove.related (required) asks whether kinship coefficients will be calculated and related samples will be removed.
#' @param kinship.coefficient should be provided if remove.related = TRUE.  Based on kingrelatedness.com/manual, a kinship coefficient range of >0.354 corresponds to duplicate/MZ twins, [0.177, 0.354] to 1st-degree relationships, [0.0844, 0.177] to 2nd-degree relationships, and [0.0442, 0.0884] to 3rd-degree relationships.
#' @param geno.value is the maximum threshold of missingness per variant.
#' @param maf.value is the maximum threshold of missingness per variant.
#' @param mind.value is the percentage threshold of individuals with missing data
#' @param limit.LD asks whether the function should prune the markers that are in linkage equilibrium. 
#' @param indep.pairwise.kb should be set if limit.LD is TRUE. This is the window size in variant count or kilobase units. 
#' @param indep.pairwise.ct should be set if limit.LD is TRUE. This is the variant count to shift the window at each end of step.
#' @param indep.pairwise.r2 should be set if limit.LD is TRUE. This is the pairwise threshold. At each step, pairs of variants in the current window with squared correlation greater than the threshold are noted and the variants are pruned from the window until no such pairs remain.
#' @examples
#' filter.dataset("myvcf.vcf", plink.files = FALSE, remove.related = TRUE, kinship.coefficient = 0.0844, geno.value = 0.5, maf.value = 0.05, mind.value = 0.1, limit.LD = FALSE)
#' filter.dataset("myvcf.vcf", plink.files = FALSE, remove.related = FALSE, geno.value = 0.5, maf.value = 0.05, mind.value = 0.1, limit.LD = TRUE, indep.pairwise.kb = 200, indep.pairwise.ct = 30, indep.pairwise.r2 = 0.5)
#' filter.dataset(plink.files = TRUE, bed.file = "mybed.bed", bim.file = "mybim.bim", fam.file = "myfam.fam", remove.related = FALSE, geno.value = 0.5, maf.value = 0.05, mind.value = 0.1, limit.LD = TRUE, indep.pairwise.kb = 300, indep.pairwise.ct = 25, indep.pairwise.r2 = 0.4)
#' @export

filter.dataset <- function(
    input.file = input, 
    plink.files = FALSE,
    bed.file = bed,
    bim.file = bim,
    fam.file = fam,
    remove.related = TRUE,
    kinship.coefficient = value,
    geno.value = geno,
    maf.value = maf,
    mind.value = mind,
    limit.LD = FALSE,
    indep.pairwise.kb = window.size,
    indep.pairwise.ct = step.size,
    indep.pairwise.r2 = r2.threshold,
    plink.exec = plink.path,
    plink2.exec = plink2.path){
  
  plink_exec <- function(PLINKoptions = "") system(paste(plink.exec,PLINKoptions)) #specify path to plink1.9
  plink_exec2 <- function(PLINKoptions = "") system(paste(plink2.exec,PLINKoptions)) #specify path to plink2.0
  
  #ensure parameters are read as integers
  kinship.coefficient <- as.integer(kinship.coefficient)
  geno.value <- as.integer(geno.value)
  maf.value <- as.integer(maf.value)
  mind.value <- as.integer(mind.value)
  indep.pairwise.kb <- as.integer(indep.pairwise.kb)
  indep.pairwise.ct <- as.integer(indep.pairwise.ct)
  indep.pairwise.r2 <- as.integer(indep.pairwise.r2)
  
  if(plink.files == TRUE){
  input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file) # read in plink files
  
  if(remove.related == TRUE){
    # checking related samples
    plink_exec2(stringr::str_c("--bfile ", 
                               input, 
                               " --make-king triangle bin --out related"))
    
    plink_exec2(stringr::str_c("--bfile ",
                               input, 
                               " --king-cutoff related ", 
                               kinship.coefficient, 
                               " --make-bed --out unrelated"))
    
    plink_exec(stringr::str_c(
      "--bfile unrelated --geno ",
      geno.value,
      " --mind ", 
      mind.value,
      " --maf ",
      maf.value,
      " --make-bed --out filtered"))
    
    if(limit.LD == TRUE){
      plink_exec(stringr::str_c(
        "--bfile filtered --indep-pairwise ",
        indep.pairwise.kb, " ",
        indep.pairwise.ct, " ",
        indep.pairwise.r2,
        " --recode vcf --out filtered.pruned"))
    }
    
  } else if(remove.related == FALSE){

    plink_exec(stringr::str_c(
      "--bfile ",
      input,
      " --geno ",
      geno.value,
      " --mind ", 
      mind.value,
      " --maf ",
      maf.value,
      " --make-bed --out filtered"))
    
    if(limit.LD == TRUE){
      plink_exec(stringr::str_c(
        "--bfile filtered --indep-pairwise ",
        indep.pairwise.kb, " ",
        indep.pairwise.ct, " ",
        indep.pairwise.r2, 
        " --recode vcf --out filtered.pruned"))
      
    }
    
  } else {
    stop("Indicate if remove.related is TRUE or FALSE.")
  }

  
  } else if(plink.files == FALSE){
    
      if(remove.related == TRUE){
        if(tools::file_ext(input.file) == "vcf"){
          plink_exec2(stringr::str_c(
            "--vcf ", input.file, 
            " --make-king triangle bin --out related"))
          
          plink_exec2(stringr::str_c("--vcf ", input.file, 
                                     " --king-cutoff related ", 
                                     kinship.coefficient, 
                                     " --make-bed --out unrelated"))
          
          # for bcf input files
        } else if(tools::file_ext(input.file) == "bcf"){
          plink_exec2(stringr::str_c(
            "--bcf ", input.file, 
            " --make-king triangle bin --out related"))
          
          plink_exec2(stringr::str_c("--bcf ", input.file, 
                                     " --king-cutoff related ", 
                                     kinship.coefficient, 
                                     " --make-bed --out unrelated"))
        } else {
          stop("Non-plink files are neither VCF or BCF files. Input file is in an unacceptable format.")
        }
        
        plink_exec(stringr::str_c(
          "--bfile unrelated --geno ",
          geno.value,
          " --mind ", 
          mind.value,
          " --maf ",
          maf.value,
          " --make-bed --out filtered"))
        
        if(limit.LD == TRUE){
          plink_exec(stringr::str_c(
            "--bfile filtered --indep-pairwise ",
            indep.pairwise.kb, " ",
            indep.pairwise.ct, " ",
            indep.pairwise.r2, 
            " --recode vcf --out filtered.pruned"))
        }
        
      } else if(remove.related == FALSE){
        if(tools::file_ext(input.file) == "vcf"){
        
        plink_exec(stringr::str_c(
          "--vcf ",
          input.file,
          " --geno ",
          geno.value,
          " --mind ", 
          mind.value,
          " --maf ",
          maf.value,
          " --make-bed --out filtered"))
        
        if(limit.LD == TRUE){
          plink_exec(stringr::str_c(
            "--bfile filtered --indep-pairwise ",
            indep.pairwise.kb, " ",
            indep.pairwise.ct, " ",
            indep.pairwise.r2,
            " --recode vcf --out filtered.pruned"))
        }
        } else if(tools::file_ext(input.file) == "bcf"){
          
          plink_exec(stringr::str_c(
            "--bcf ",
            input.file,
            " --geno ",
            geno.value,
            " --mind ", 
            mind.value,
            " --maf ",
            maf.value,
            " --make-bed --out filtered"))
          
          if(limit.LD == TRUE){
            plink_exec(stringr::str_c(
              "--bfile filtered --indep-pairwise ",
              indep.pairwise.kb, " ",
              indep.pairwise.ct, " ",
              indep.pairwise.r2, " ",
              " --recode vcf --out filtered.pruned"))
          }
      }
      }  else {
        stop("Indicate if remove.related is TRUE or FALSE.")

  } 
  } else {
    stop("Indicate if plink.files is TRUE or FALSE")
  }
  }