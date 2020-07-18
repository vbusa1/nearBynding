# nearBynding

`nearBynding` is an R package that discerns RNA structure proximal to protein 
binding within regions of the transcriptome defined by the user. Input CLIP 
protein-binding data can either be in aligned BAM or peak-called BED/bedGraph 
formats. RNA structure can either be internally calculated via CapR or can be 
provided by the user as a BED/bedGraph. RNA structure binding profiles can be 
visually and mathematically compared between proteins across multiple formats.

## Installation

The most up-to-date version can be loaded to R via

```
devtools::install_git("vbusa1/nearBynding"")
library(nearBynding)
```

## External Software Dependencies

`nearBynding` requires three external softwares. Add all dependency directories 
to your PATH after installation.

#### bedtools

bedtools is available for installation 
[here.](https://bedtools.readthedocs.io/en/latest/content/installation.html)

Installation instructions will vary by operating system.

#### CapR

Download the zip file from the [github repository](https://github.com/fukunagatsu/CapR), 
unzip the file, and move it to a directory where you want to permanently store 
the function.

In the command line, access the folder where the unzipped file is stored.

```
cd CapR-master
make
./CapR
```

If installation is successful, the final line will output

`Error: The number of argument is invalid.`

#### StereoGene

Download the zip file from the [github repository](https://github.com/favorov/stereogene), 
unzip the file, and move it to a directory where you want to permanently store 
the function.

In the command line, access the folder where the unzipped file is stored.

```
cd stereogene-master
cd src
make
./stereogene -h
```

If installation is successful, the final line will output a menu of argument options.
