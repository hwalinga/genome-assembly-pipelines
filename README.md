Pipelines for making genome assemblies

# Updates

#### Update 30 Nov

You can now restart the ont pipeline and it will continue where it left off.
It will however continue after the last created folder. So if you know
there was a crash in the last folder. Just `rm -rf` that folder.

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

### Core (You most likely have this already installed)
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
-d [DIR]
    You can also provide the directory with all the fastq files with this option.
    If you leave this without any argument you will be prompted (zenity required).
    The files in this directoy must have the "*1.fq.gz" and "2.fq.gz" suffix.
-p
    This option will prompt you automatically for all the options
    (zenity required)
    (Currently not very well implemented)
--help,-h
    Print this help and exit.
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

(This will assume the files are matching the pattern '\*1.fq.gz' and '\*2.fq.gz'.)

NB. Known bug 1. You can have the wildcard (\*) within a folder name, but in that case
your folders cannot contain any spaces.

NB. Known bug 2. Also, having spaces in your work
directory and source files will error out the program.

NB. Known bug 3. You can also not combine spaces in the basename of the file with globs.

(As you realise, spaces are annoying. You can replace spaces with underscores with for example:)

```
for i in *; do mv "$i" "${i/ /_}"; done
```

## ONT

#### Important

The ONT pipeline makes some assumptions which you might want to change
in the script.

The first and probably most important one is that you must give the correct
kit to use for the demultiplexing. Change the line `KIT=".."` for the correct
kit you used. Run `qcat --list-kits`, to see what are avaible.

Next, the demultiplexing is only keeps on going with barcodes that have
minimum amount of data available. You can change this in the line that states:
`find demultiplex ! -name none.fastq -type f -size +150k ...` where the size
the minimum size is the barcode data must contain.

Last but not least, flye is optimized for a certain contig lenght. You must
change that if you assemble something with a different length.
The *-g* option in flye indicates the predicted contig length.

Usage:
ont.bash [--options] "FolderPath" OR/AND "FastqFiles"
--test:
    Show commands to run, but not execute them.
--keep:
    This flag will make sure the inbetween results will be kept.
    This way you can resume the pipeline when it crashed unexpectedly.
    (NB. Currently --keep is not implemented.)
    But you can rerun the script and it will run from the last result.
    However, it will also run if that result is incomplete. So if you know
    that there was a crash in the last result, it might be best to remove this
    folder so that is the first result.
--nocov:
    Do not plot the coverage plots.
-i [FILE]
    Instead of supplying as an argument you can pass the fastq file
    with the -i option. Note that this way the file will not be copied.
    (Copying might be desirable if it is on an unstable filesystem.)
    If you leave this without any argument you will be prompted (zenity required).
-d [DIR]
    You can also provide the directory with all the fastq files with this option.
    If you leave this without any argument you will be prompted (zenity required).
-p
    This option will prompt you automatically for all the options
    (zenity required)
    (Currently not very well implemented)
--help,-h
    Print this help and exit.
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
you making your life easier. Run `showhelp` and see under the section
"aliases for assembly / bulk drive interaction". You can also find these and the pdf of showhelp at
https://github.com/hwalinga/ont-linux-cluster-setup/

## NB
For long jobs please note that you don't have to leave your terminal open,
but have to make use of the `screen` program, that is also explained in the
section "running long jobs" of the `showhelp`.
