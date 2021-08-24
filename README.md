
# Msuite2: Multi-mode DNA methylation data analysis suite
Msuite2 is the successor of Msuite.<br />
Version 2.0.0, Aug 2021<br />
Authors: Lishi Li, Yunyun An, Pengxiang Yuan, Li Ma, Xin Jin, Xin Hong, Kun Sun<br />
Software implemented by Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please read the license file under `license` directory.

---

## Installation
`Msuite2` is written in `Perl` and `R` for Linux/Unix platform. To run `Msuite2` you need a Linux/Unix
machine with `Bash 4 (or higher)`, `Perl 5.10 (or higher)` and `R 2.10 (or higher)` installed.

This source package contains pre-compiled executable files using `G++ v4.8.5` for Linux x86_64 system.
If you could not run the analysis normally (which is usually caused by low version of `libc++` library),
or you want to build a different version optimized for your system, you can re-compile the programs:
```
user@linux$ make clean && make
```

Note that `Msuite2` depends on the following software:

* [bowtie2](https://github.com/BenLangmead/bowtie2 "bowtie2")
* [samtools](http://samtools.sourceforge.net/ "samtools")

Please install them properly and make sure that they are included in your `PATH`.
In addition, please make sure that the version of your `g++` compiler is higher than 4.8
(you can use `g++ -v` to check it).

Before running `Msuite2`, genome indices must be built. To this end, we have prepared a utility named
`build.index.sh` under the `build.index` directory. To use it, you need to prepare the genome sequence
(either in one multi-fasta file or a directory containing the sequences for individual chromosomes) and
RefSeq annotation for your genome (we have included files for mm10 and hg38 in this package).
You can download the annotations for other species/genome versions from the
[UCSC genome browser](http://genome.ucsc.edu/ "UCSC Genome Browser").
(The RefSeq annotation is used to profile the methylation level around Transcription Start Sites as a
quick quality control of your data.)

Then, you can build the genome indices using the following command:
```
user@linux$ build.index/build.index.sh GENOME.FA(or GENOME.DIR) REFSEQ.txt Genome.ID
```
Note that this utility will automatically incorporate the Lambda genome to build the genome indices.
`Gzip` or `Bzip2` compression of `GENOME.FA` and `REFSEQ.txt` are also supported, but the files must
contain the corresponding suffix (i.e., `REFSEQ.txt.gz` for `Gzip` compressed and `REFSEQ.txt.bz2` for
`Bzip2` compressed file). The `Genome.ID` is an identifier that you specified to name your genome and
the indices will be written to the `index` directory under the root of `Msuite2`. You can add as many
genomes to `Msuite2` as you need.

## Run Msuite2
The main program is `msuite2`. You can add its path to your `.bashrc` file under the `PATH` variable
to call it from anywhere, or you can run the following command to add it to your current session:
```
user@linux$ PATH=$PATH:$PWD
```

Call `Msuite2` without any parameters to see the usage (or use '-h' option):
```
########## Msuite2: Multi-mode DNA methylation data analysis suite ##########

Author : Kun Sun (sunkun@szbl.ac.cn)
Version: v2.0.0 (Aug 2021)

Usage: msuite [options] -x index -1/-U Read1.fq [ -2 Read2.fq ] -o out.dir

Compulsory parameters:

  -1/-U Read1.fq   Specify the path to the files containing read 1
                   If your data is Paired-end, specify read 2 files using '-2' option
                   Note that if -U is used, '-2' will be ignored no matter it's set or not

                   If you have multiple files, use ',' to separate them or use '*' syntax
                   (Note that single quotation marks are required for '*' syntax)

  -x index         Specify the genome index
                   Please refer to README file on how to build index for Msuite

  -o out.dir       Specify the output directory
                   Note that your specified directory will be created if it does not exist
                   otherwise the files under that directory could get over-written


Optional parameters:

  -2 Read2.fq      Specify the path to the file containing read 2
                   Use this parameter if your data is generated in paired-end mode

                   If you have multiple files, use ',' to separate them or use '*' syntax
                   Note that all files must be correctly paired in '-1' and '-2'

  -3               Use 3-letter alignment (default)
  -4               Use 4-letter alignment
                   Note that the above two options are mutually exclusive

  -m BS/TAPS       Specify the library protocol (default: BS)
                   Note that only 'TAPS' and 'BS' are acceptable

  -c cycle         Specify the seqeuencing cycles of the data (default: auto-detect)

  -k kit           Specify the library preparation kit (default: illumina)
                   Note that the current version supports 'illumina', 'nextera' and 'bgi'

  --phred33        Read cycle quality scores are in Phred33 format (default)
  --phred64        Read cycle quality scores are in Phred64 format
                   Note that the above two options are mutually exclusive

  -q score         The minimum quality score to keep the cycle (default: 20)
                   Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred

                   Sometimes quality scores start from 35 ('#') in the FASTQ files,
                   in this case you could adjust '-q' option, e.g., '--phred33 -q 22'

  --minsize size   Minimum read size to be kept for alignment (default: 36)

  --cut-r1-head N  Cut the head N cycles in read1 (default: 0)
  --cut-r1-tail N  Cut the tail N cycles in read1 (default: 0)
  --cut-r2-head N  Cut the head N cycles in read2 (default: 0)
  --cut-r2-tail N  Cut the tail N cycles in read2 (default: 0)

                   Note that the total cut basepairs in read 1 and 2 must be the same
                   (i.e., cut-r1-head + cut-r1-tail must equal to cut-r2-head + cut-r2-tail)

                   In BS-seq, read2 starts from the 3'-end which frequently suffer from
                   DNA damage issues, and the DNA repair step in library prepraration usually
                   uses un-methylated Cytosines which lead to bias in DNA methylation calling
                   as commonly seen in the M-bias plot.

  --minins MIN     Minimum insert size (default: 0)
  --maxins MAX     Maximum insert size (default: 1000)
                   Note that the above two options will be ignored for Single-End data

  --align-only     Stop after alignment (i.e., do not perform DNA methylation call and
                   visualization around TSS; default: not set)

  --CpH            Set this flag to call methylation status of CpH sites (default: not set)

  -p threads       Specify how many threads should be used (default: use all threads)

  -h/--help        Show this help information and quit
  -v/--version     Show the software version and quit

Please refer to README file for more information.

```

**IMPORTANT NOTE**: If your data is generated using BS-seq protocol, you MUST use the 3-letter mode and set
`-m BS`. 4-letter mode ONLY supports processing of TAPS/5hmC-CATCH data where the non-CpG methylation is
very low (e.g., most somatic tissues in human). In addition, Msuite2 could directly analyze the data generated
by ATAC-me or similar protocols via setting `-k nextera`. You can use '-c cycle' option to control the cycles
that you want to analyze (if you do not want to analyze all cycles for some reason).


### Example 1
Your data is generated using TAPS protocol in 75 bp * 2 (paired-end) mode, and you want to align your data to
the hg19 reference genome in 4-letter mode, and you want to use 16 threads to speed-up the analysis, then you
can run:
```
user@linux$ Msuite2 -1 /path/to/read1.fq -2 /path/to/read2.fq -x hg19 \
                   -4 -m TAPS -p 16 -o /path/to/output/dir
```

### Example 2
Your data is generated using BS-seq protocol in 100 bp * 1 (single-end) mode while you only want to analyze the
first 75 bp of your reads (e.g., due to sequencing quality considerations), and you want to align your data to
the mm10 reference genome (note that you MUST use 3-letter mode here), and use 32 threads to speed-up the analysis,
	then you can run:
```
user@linux$ Msuite2 -1 /path/to/read1.fq.gz -x mm10 -c 75 \
                   -3 -m BS -p 32 -o /path/to/output/dir
```

### Example 3
If your data is generated using BS protocol in paired-end mode, and you have 3 lanes of data, you want to align
your data to the hg19 reference genome, you want to skip the head/tail 5/10 cycles in both reads to suppress the
issues by DNA overhang, and you want to use 48 threads to speed-up the analysis, then you can run:
```
user@linux$ Msuite2 -1 /path/to/lane1.read1.fq.gz,/path/to/lane2.read1.fq.gz,/path/to/lane3.read1.fq.gz \
                   -2 /path/to/lane1.read2.fq.gz,/path/to/lane2.read2.fq.gz,/path/to/lane3.read2.fq.gz \
                   --cut-r1-head 5 --cut-r1-tail 10 --cut-r2-head 5 --cut-r2-tail 10 \
				   -x hg19 -p 48 -o /path/to/output/dir
```

If you want to use add all the '.fq' files in your path, you can use the `*` syntax:
```
user@linux$ Msuite2 -1 '/path/to/lane*.read1.fq.gz' \
                   -2 '/path/to/lane*.read2.fq.gz' \
				   --cut-r1-head 5 --cut-r1-tail 10 --cut-r2-head 5 --cut-r2-tail 10 \
                   -x hg19 -p 48 -o /path/to/output/dir
```
Note that the single quotation mark is essential to protect the '\*' syntax from been extracted by your shell.

<br />
`Msuite2` will check the data and dependent programs then generate a `makefile` under `/path/to/output/dir`
('-o' option). Then you can go to `/path/to/output/dir` and run `make` to perform the analysis:
```
user@linux$ cd /path/to/output/dir; make
```

We have prepared a testing dataset under the `testing_dataset` directory. It contains *in silico* generated
reads following the TAPS (TET-assisted pyridine borane sequencing) protocol using
[SHERMAN](http://www.bioinformatics.babraham.ac.uk/projects/sherman/) software (key parameters: C-&gt;T
conversion rate: 20% for CpG sites, C-&gt;T conversion rate in CpH sites: 0.5%, error rate: 0.1%).
Note that the reads are restricted to CT/GA-rich regions (CT proportion >=80% or GA proportion >=80%) to
illustrate the advantage of 4- over 3-letter alignment. We also have prepared a work shell to run `Msuite2`
on this dataset using both 3- and 4-letter modes:
```
user@linux$ ./run_testing_dataset.sh
```
Note that this script will automatically build the indices for hg19 genome if you have not done this before.

You can compare the performance of 3-letter and 4-letter alignments by inspecting the outputs, which will be
written to `testing_dataset/Msuite2.Mode3/` and `testing_dataset/Msuite2.Mode4/`.


## Outputs explanation
`Msuite2` outputs all the results for the given region in the directory specified by `-o OUTDIR` option.
`Msuite2` will write the analysis report into a HTML file named `Msuite2.report/index.html`, which records the
essential statistics and visualizations for quality control, mappability, overall methylation level on CpG
sites, M-bias plot, and conversion rate (estimated using reads mapped to the Lamda genome).

The alignment results are recorded in the file `Msuite2.final.bam` (in standard BAM format) and "Msuite2.rmdup.sam"
(in standard SAM format). The methylation calls are recorded in the file `Msuite2.CpG.meth.call`,
`Msuite2.CpH.meth.call` and `Msuite2.CpG.meth.bedgraph`.

You can run `make clean` in the OUTDIR to delete the intermediate files to save storage space.


## Utilities
### Mviewer
`Msuite2` contains a visualization tool named `Mviewer`, adapted from the authors' previous
[BSviewer](http://sunlab.cpy.cuhk.edu.hk/BSviewer/) software. It is specially optimized to be compatiable
with `Msuite2` alignment results and provides nucleotide-level, genotype-preserved DNA methylation data
visualization. For more information, please refer to README file in `Mviewer` directory.


### Others
`Msuite2` also provides other utilities under the `util` directory.<br />

The `profile.meth.pl` program is designed to summarize the methylome into bins. You can use it to prepare data
for `Circos` plots.<br />

The `extract.meth.in.region` program is desgined to extract the covered CpG sites, C-count and T-count in the
given regions (e.g., CpG islands, promoters).<br />

The `pe_bam2bed.pl` and `se_bam2bed.pl` are designed to translate the aligned BAM file into BED format file, and
`bed2wig` is designed to translate BED file into WIG files (e.g., for coverage profiles).<br />

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
Msuite2 package is freely available at
[https://github.com/hellosunking/Msuite2/](https://github.com/hellosunking/Msuite2/ "Msuite2 @ Github").


