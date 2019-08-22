# genome-assembly-pipelines
Pipelines for making genome assemblies

# About the pipelines
These pipelines are written in pure bash, that has advantages and disadvanteges.
It might have been a better idea to write them in a workflow management system
like Snakemake or Nextflow.

# Installation
Installation is simple: Clone this folder `git clone https://github.com/hwalinga/genome-assembly-pipelines`. Or `wget` the specific file you want with for example `https://raw.githubusercontent.com/hwalinga/genome-assembly-pipelines/master/ngs.bash`. Running the program is simple, just use `bash path/to/pipeline.bash`.

Almost all Linux computers will have a sufficient modern bash running to have no problems. bash on a Mac is probably too old to run and needs to updated.

You still have to make sure that the you have installed the dependencies for the assemblying itself. You can install this with Anaconda or Homebrew. (SOAPnuke has to be installed manually)

### Core (You most likely have this installed already installed)
- free
- ps
- printf
- awk
- perl

### Common
- GNU Parallel

### Next Generationg Sequencing (NGS)
- SOAPnuke
- fastp
- spades
- wkhtmltopdf (optional)
- seqtk (for phage assembly)

### Oxford Nanopore Technologies (ONT)
- qcat
- filtlong
- flye
- medaka
- racon

### Coverage plots
* GNU Plot (Version >= 5.0)
* samtools (Verion >= 1.0)
* minimap2

# Usage

For both pipelines applies: Make a workspace where you output will end up with,
and run `bash path/to/pipeline.bash with the correct argments`, when in this
workspace. (`bash path/to/pipeline.bash --help` can help you with how to use.)

## NGS

## ONT
