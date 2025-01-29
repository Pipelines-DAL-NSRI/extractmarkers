# extractmarkers

## Description
Function for extracting markers from PLINK, VCF, or BCF files. 

## Parameters
@param input.file is either a VCF or BCF file present in the working directory. File path is also accepted, but not suggested. Can skip if using PLINK files.
@param plink.files (required) asks whether the input files are in PLINK format. Default is FALSE.
@param bed.file should be provided if plink.files = TRUE
@param bim.file should be provided if plink.files = TRUE
@param fam.file should be provided if plink.files = TRUE
@param snps (required) is a text file containing the list of markers to be extracted. Each line should correspond to one marker.
@param output (required) is the user's desired output name.
@param plink.exec (required) takes in the file path for the plink executable. Suggestion is to have plink.exe in the working directory.

## Use
```
extractmarkers(input.file = "mygenotype.vcf", plink.files = FALSE, snps = "listsnps.txt", output = "success", plink.exec = "./plink.exe")

extractmarkers(plink.files = TRUE, bed.file = "data.bed", bim.file = "data.bim", fam.file = "data.fam", snps = "listsnps.txt", output = "success", plink.exec = "./softwares/plink.exe")
```
