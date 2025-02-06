# extractmarkers

## Description
The 'extractmarkers' function extracts a set of markers, either using their rsID or chromosome position, from PLINK, VCF, or BCF files. In using rsIDs for extraction, a text file containing one marker per line is accepted. When using GRCh38/GRCh37 positions, the text file should contain 4 columns: the first column should contain the chromosome number only (e.g. "1" and not "chr1"); the second column should contain the starting bp position; the third column the final bp position; and the fourth column should specify the desired output name with no file extension (e.g. "snpsfromchr1").

## Function
```
extractmarkers(
    input.file = input, 
    plink.files = FALSE,
    bed.file = bed,
    bim.file = bim,
    fam.file = fam,
    markerID = FALSE,
    snps.list = snps.file, 
    output = output.name,
    plink.exec = plink.path)
```

## Parameters
@param *input.file* is either a VCF or BCF file present in the working directory. File path is also accepted, but not suggested. Can skip if using PLINK files.

@param *plink.files* (required) asks whether the input files are in PLINK format. Default is FALSE.

@param *bed.file* should be provided if plink.files = TRUE

@param *bim.file* should be provided if plink.files = TRUE

@param *fam.file* should be provided if plink.files = TRUE

@param *markerID* (required) asks whether the markers will be extracted by rsID. Set to FALSE if using positions.

@param *snps.list* (required) is a text file containing the list of markers to be extracted. 
#' If using rsIDs, each line should correspond to one marker. If using POS (positions), there should be 4 columns: [1] chromosome numner (integer), [2] starting base position, [3] final base position, and [4] output file name. 

@param *output* (required) is the user's desired output name.

@param *plink.exec* (required) takes in the file path for the plink executable. Suggestion is to have plink.exe in the working directory.

## Use
```
extractmarkers(input.file = "mygenotype.vcf", plink.files = FALSE, markerID = TRUE, snps.list = "listsnps.txt", output = "success", plink.exec = "./plink.exe")

extractmarkers(plink.files = TRUE, bed.file = "data.bed", bim.file = "data.bim", fam.file = "data.fam", markerID = FALSE, snps.list = "listsnps.txt", plink.exec = "./softwares/plink.exe")
```
