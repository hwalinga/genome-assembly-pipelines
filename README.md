# genome-assembly-pipelines
Pipelines for making genome assemblies

# About the pipelines
These pipelines are written in pure bash, that has advantages and disadvanteges.
It might have been a better idea to write them in a workflow management system
like Snakemake or Nextflow.

# Installation
Installation is simple: Clone this folder, you can just have this in your homefolder.
```
cd  # Go to the home folder.
git clone https://github.com/hwalinga/genome-assembly-pipelines  # Clone this directory.
```

Or `wget` the specific file you want with for example `https://raw.githubusercontent.com/hwalinga/genome-assembly-pipelines/master/ngs.bash`.

Running the program is simple, just use `bash path/to/pipeline.bash`.

(For advanced usage you can make the scripts available with the `$PATH` variable
and `chmod +x` them.)

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

### Misc.
* zenity (For the interactive use)

# Usage

For both pipelines applies: Make a workspace where you output will end up with,
and run `bash path/to/pipeline.bash with the correct argments`, when in this
workspace. (`bash path/to/pipeline.bash --help` can help you with how to use,
but the help is also printed below.)

## NGS

```
Usage:
ngs.bash [--options] "SOURCE1_glob" "SOURCE2_glob"
Provide the source1 and source2 interleaved glob patterns quoted, so that
they will not be expanded before the program can read them.
--test:
    Show commands to run, but not execute them.
--nomove:
    If source files are on a different parition than the current directory,
    they will first will be moved to this partition.
    If you have a reliable connection between the two drives, you can disable
    that default behavior with this flag.
    (Results still end up on the current directory regardless of this option.)
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly.
    (NB. However, at the moment there are no checks in place that check
    if a part of the pipeline has already been completed in a previous run.)
--bacteria,--phage,--both:
    Set one of this option if you have bacteria or phage DNA.
    In this case the difference with a phage assembly is that it samples the
    data, so that the assembly is faster. There isn't any difference further
    at the moment.
--nocov:
    Do not plot the coverage plots.
--help,-h
    Plot this help and exit.
```

For NGS you will have paired-end reads. The paths to this data must be quoted.
(Or at least the glob pattern must be quoted.)
For example:

```
bash ~/genome-assembly-pipelines/ngs.bash "path/to/files/*_1.fq.gz" "path/to/files/*_2.fq.gz"
```

For easy use, just let the program prompt you for the directory with the fastq files.

```
bash ~/genome-assembly-pipelines/ont.bash -d
```

(This will assume the files are matching the pattern '*1.fq.gz' and '2.fq.gz'.)

NB. Known bug. You can have the wildcard (\*) within a folder name, but in that case
your folders, cannot contain any spaces. Also, having spaces in your work
directory and source files will error out the program.

## ONT

```
Usage:
ont.bash [--options] "FolderPath" OR/AND "FastqFiles"
--test:
    Show commands to run, but not execute them.
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly.
    (Currently not implemented.)
--nocov:
    Do not plot the coverage plots.
-i [FILE]
    Instead of supplying as an argument you can pass the fastq file
    with the -i option. Note that this way the file will not be copied.
    (Copying might be desirable if it is on an unstable filesystem.)
--help,-h
    Plot this help and exit.
```

The ONT pipeline starts with the fastq files.
(Basecalling the fast5 files should be done on a powerfull machine with GPU acceleration.)
You can provide a list of these files or the folder in which to find the
fastq files as arguments.
For example:

```
bash ~/genome-assembly-pipelines/ont.bash /path/to/fastq/files/fastq_pass/
```

For easy use, just let the program prompt you for the directory with the fastq files.

```
bash ~/genome-assembly-pipelines/ont.bash -d
```

## sb-ont users
If you are on the sb-ont machine, there are specific aliases that will help
you making your life easier. Run showhelp and see under the section
"assembly aliases". You can also find these and the pdf of showhelp at
https://github.com/hwalinga/ont-linux-cluster-setup/
