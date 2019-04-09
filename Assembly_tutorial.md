Transcriptome Assembly Practical
--

#### The following details the steps involved in:

- Generating a _de novo_ RNA-Seq assembly using the Oyster River Protocol
- Evaluating the quality of the assembly
  - BUSCO and TransRate
- Quantifying transcript expression levels
  - Using Salmon
- Functionally annotating transcripts (won't get here, but see code)
  - dammit
- Predicting coding regions (won't get here, but see code)
  - TransDecoder

#### Oyster River Protocol

One of the nice things about the ORP is that it does, automatically, the assembly, evaluation, and quantification of the assembly. I will more fully introduce these procedures in the lecture part of the day.

All required software and data are provided pre-installed on a Docker VM. For when you get home, and want to install this on your own machines, see http://oyster-river-protocol.readthedocs.io/en/latest/aws_setup.html and https://oyster-river-protocol.readthedocs.io/en/latest/docker_install.html

The workshop materials here expect that you have some familiarity with UNIX. If you need a refresher, the Command Line Bootcamp is good (http://rik.smith-unna.com/command_line_bootcamp/). For extended prep, see DataCarpentry and SoftwareCarpentry organizations.

#### About the Tutorial

This tutorial uses a very small dataset, so that everyone can have the opportunity to complete an assembly without using a truly massive amount of compute. **BUT**, if you can assembly the small dataset, you can assembly a bigger one, too, assuming you have a big enough computer.

#### About the data

The data used in this tutorial are from SRR1221220, which is a RNAseq dataset from the Japanese fire belly newt (_Cynops pyrrhogaster_). This is what I did, before the class.

- Assembled the full dataset using the ORP
- Quantitated gene expression using Salmon
- Extracted out 50 of the highest expression transcripts (numbers 451-500 to be exact).
- Mapped the reads to the fasta file of those 50 contigs using bwa mem, producing a SAM file.
- Used Picard to move from SAM to fastq
- Extracted a random 80,000 read pairs from the larger fastq using seqtk.

#### Run Docker image

What you'll notice is that we'll just launch the image in our terminal - no using the web-browser interface. I'm sure Bastian will at this point tell us why this might not be a good choice, but it works for me, and should for you, too.

```
docker run -it -v /media/penelopeprime/RNA-Sequence\ \
Analysis:/home/training/share:ro \
macmaneslab/orp:2.2.5_ebi2019 bash
```

At this point, your terminal prompt should look something like this `training@1f518d44b503:~$`



#### Tutorial Begins

The 1st thing you should always do when logging in to a new machine is explore the directory structure. Where are the data, where are programs installed, etc.

Make a directory that will contain all of your assembly materials. In general, it's smart practice to have individual folders for each step in your bioinformatics pipeline.

```
mkdir $HOME/assembly_practical && cd $HOME/assembly_practical
```

What is this `&&` thing?? It basically serves to link two commands together, **IF** the 1st one succeeds. So, `do command 1, do command 2 if 1 succeeds, do command 3 if command 2 succeeds`

#### Assemble using the ORP

Activate the conda environment:

```
conda activate orp
```

At this point, your command prompt should look something like  `(orp) training@1312b3d53c18:~$`. Note the prefix `(orp)`.


```
$HOME/Oyster_River_Protocol/oyster.mk \
STRAND=RF \
TPM_FILT=1 \
MEM=10 \
CPU=8 \
READ1=$HOME/Oyster_River_Protocol/sampledata/test.1.fq.gz \
READ2=$HOME/Oyster_River_Protocol/sampledata/test.2.fq.gz \
RUNOUT=sampledata
```

the final assembly will be at `assemblies/sampledata.ORP.fasta`

Let's undestand this command:

- `$HOME/Oyster_River_Protocol/oyster.mk` the ORP is written as a Makefile, which is a nice way to organize computational pipelines. With few exceptions, if your run fails mid-way through the process, restarting is using the same command will pick up at the point at which it failed.
- `STRAND` Was the library prepared using a strand-specific approach. In 2019, most libraries are, and `RF` is the most common.
- `TPM_FILT` Do you want to remove* lowly expressed, mostly non-biological transcripts. (you probably do)
- `MEM` How much RAM do you require (in Gb)? I usually set this to about 10% less than what the computer has.
- `CPU` use all that you have.
- `READ1` and `READ2` The ORP requires PE reads, sorry. Assembling with SE reads is not really worth it.. The reads can be gzipped. Always safe(r) to use the full path to the reads.
- `RUNOUT` Name the output. This name will be included in the final assembly, so choose wisely.  No special characters (e.g., `\|*?/`)
- appending `--dry-run` to the end of the command will print out to your screen the commands that will be run on your dataset, but they won't actually be run.   


What is the ORP doing. We'll talk about this in the lecture part of the class.

1. Trim adapters and low quality bases (Phred<2)
2. Correct the reads
3. Asssembly using Trinity, SPAdes k=55, SPAdes k=75, TransABySS
4. For iso-groups from the 4 assemblies, and pick the best member of each isogroup.
5. BLAST assemblies *using Diamond to Swiss-Prot* and make sure that we have all the unique stuff from each assembly.
6. cd-hit the assembly
7. Run Salmon on the assembly
8. Optionally, filter very lowly expressed trasnscripts, and make sure we are not throwing out "real" stuff.
9. Run BUSCO
10. Run TransRate.
11. Print a report


How many resources do you need?

- **RAM** The amount of RAM you need scales with the number of unique kmers. The number of unique kmers is positively correlated with number of reads, which is a far more accessible number. In general, you need about 1GB RAM per million reads.
- **CPUs** More is better, but machines with between 12 and 64 cores are common. Note the assembly process is not able to leverage MPI.

##### Evaluating content using BUSCO and TransRate

So you have an assembly, now how good is it? One way is to look at the assembly content. Are all the expected genes present? Fo this part of the exercise, We're going to use a 'real' assembly.


```
mkdir $HOME/assembly_practical/assembly_eval && cd $HOME/assembly_practical/assembly_eval
```

There are some assemblies in `$HOME/share/Day3/assemblies/` and some reads in `$HOME/share/Day3/reads/`... Each classroom row is going to evaluate a different assembly.

Row1=gar
Raw2=lamprey
Row3=mouse
Row4=turtle
Row5=wallaby

Here is the command for another organism, the deer. Your command will be different!!! The reads are a subset of size 1M read pairs.

```
$HOME/Oyster_River_Protocol/report.mk \
ASSEMBLY=$HOME/share/Day3/assemblies/deer.ORP.fasta \
CPU=8 \
MEM=10 \
LINEAGE=eukaryota_odb9 \
READ1=$HOME/share/Day3/reads/deer.1.fq \
READ2=$HOME/share/Day3/reads/deer.2.fq \
RUNOUT=deer
```


Understanding this command:

- `$HOME/Oyster_River_Protocol/report.mk` the ORP is written as a Makefile, which is a nice way to organize computational pipelines. With few exceptions, if your run fails mid-way through the process, restarting is using the same command will pick up at the point at which it failed.
- `MEM` How much RAM do you require (in Gb)? I usually set this to about 10% less than what the computer has.
- `CPU` use all that you have.
- `LINEAGE` Which BUSCO database are you going to use. Here, we will use the eukaryote database.
- `READ1` and `READ2` The ORP requires PE reads, sorry. Assembling with SE reads is not really worth it.. The reads can be gzipped. Always safe(r) to use the full path to the reads.
- `RUNOUT` Name the output. This name will be included in the final assembly, so choose wisely.  No special characters (e.g., `\|*?/`)
- appending `--dry-run` to the end of the command will print out to your screen the commands that will be run on your dataset, but they won't actually be run.   



This will take about 15 minutes to run. One member for your group - write results on the board. Which assembly is the best?

**Note:** The ORP runs BUSCO automatically, and has done so for your dummy assembly. See `$HOME/assembly_practical/reports/run*ORP/short_summary*txt`. A real assembly would score **much** better, hopefully...

**Note:** The ORP runs TransRate automatically, and has done so for your dummy assembly. See `$HOME/assembly_practical/reports/transrate*/assemblies.csv`. This little code snippet will make it easier to view this file.

```
paste <(sed -n 1p $HOME/assembly_practical/reports/transrate*/assemblies.csv | tr , '\n') \
<(sed -n 2p $HOME/assembly_practical/reports/transrate*/assemblies.csv | tr , '\n')
```

```
n_seqs	104
smallest	206
largest	3036
n_bases	90282
mean_len	868.09615
n_under_200	0
n_over_1k	33
n_over_10k	0
n_with_orf	53
mean_orf_percent	75.34392
n90	357
n70	875
n50	1326
n30	1751
n10	2218
gc	0.48795
bases_n	0
proportion_n	0.0
fragments	79999
fragments_mapped	79400
p_fragments_mapped	0.99251
good_mappings	77022
p_good_mapping	0.96279
bad_mappings	2378
potential_bridges	20
bases_uncovered	1605
p_bases_uncovered	0.01778
contigs_uncovbase	32
p_contigs_uncovbase	0.30769
contigs_uncovered	2
p_contigs_uncovered	0.01923
contigs_lowcovered	15
p_contigs_lowcovered	0.14423
contigs_segmented	8
p_contigs_segmented	0.07692
score	0.57688
optimal_score	0.66115
cutoff	0.23415
weighted	5638.68556
```

##### Quantification

The ORP does this.. using Salmon..

#### Annotation using dammit (time permitting)
See http://dib-lab.github.io/dammit/install/

Installing (let's go rogue!!)

```
conda deactivate #to exit your current conda environment.  
conda create --name dammit # make a new conda environent
conda activate dammit #run the new envoronment
conda install -y -c bioconda dammit # install dammit, this will take a few minutes
```

##### I've installed the databases for you, using this command. *You don't have to do this*

```
dammit databases --database-dir $HOME/share/Day3/dammit_dbs \
--install --busco-group eukaryota
```

##### Running the Annotation (untested and potentially broken)

```
mkdir -p $HOME/assembly_practical/dammit/ && cd $HOME/assembly_practical/dammit

dammit annotate $HOME/assembly_practical/assemblies/sampledata.ORP.fasta \
--busco-group eukaryota \
--n_threads 8 \
--database-dir $HOME/share/Day3/dammit_dbs
```


#### Bibliography

1. Oyster River Protocol: doi 10.7717/peerj.5428
2. BUSCO: doi: 10.1093/molbev/msx319
3. TrasRate: doi: 10.1101/gr.196469.115
4. Salmon: doi: 10.1038/nmeth.4197
5. dammit: http://dib-lab.github.io/dammit/
