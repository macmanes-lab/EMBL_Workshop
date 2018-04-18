Transcriptome Assembly Practical
--

#### The following details the steps involved in:

- Generating a _de novo_ RNA-Seq assembly using the Oyster River Protocol
- Evaluating the quality of the assembly
  - BUSCO and TransRate
- Quantifying transcript expression levels
  - Using Salmon
- Functionally annotating transcripts (time permitting)
  - dammit
- Predicting coding regions (won't get here, but see code)
  - TransDecoder

#### Oyster River Protocol

One of the nice things about the ORP is that it does, automatically, the assembly, evaluation, and quantification of the assembly. I will more fully introduce these procedures in the lecture part of the day.

All required software and data are provided pre-installed on an Amazon EC2 AMI. For when you get home, and want to install this on your own machines, see http://oyster-river-protocol.readthedocs.io/en/latest/aws_setup.html

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

#### Tutorial Begins

The 1st thing you should always do when logging in to a new machine is explore the directory structure. Where are the data, where are programs installed, etc.

Make a directory that will contain all of your assembly materials. In general, it's smart practice to have individual folders for each step in your bioinformatics pipeline.

```
mkdir $HOME/assembly_practical && cd $HOME/assembly_practical
```

What is this `&&` thing?? It basically serves to link two commands together, **IF** the 1st one succeeds. So, `do command 1, do command 2 if 1 succeeds, do command 3 if command 2 succeeds`

##### Assemble using the ORP

There is one small bug in the Shannon assembler, which I fixed after Nico made the Docker container. Easy fix (which you should never have to do, again).


```
sed -i 's_--prefix_-q 33 --prefix_g' $HOME/Oyster_River_Protocol/software/Shannon/run_quorum.py
```

```
$HOME/Oyster_River_Protocol/oyster.mk main \
MEM=7 \
CPU=2 \
READ1=$HOME/share/Day3/read.1.fastq \
READ2=$HOME/share/Day3/read.2.fastq \
RUNOUT=ORPtest_YOURNAME
```

Let's unpack this:

- `$HOME/Oyster_River_Protocol/oyster.mk` the ORP is written as a Makefile, which is a nice way to organize computational pipelines. With few exceptions, if your run fails mid-way through the process, restarting is using the same command will pick up at the point at which it failed.
- `main` Do the whole protocol. If you specified just `trim`
 just the trimming step would be done, and nothing else.
- `MEM` How much RAM do you require (in Gb)? I usually set this to about 10% less than what the computer has.
- 'CPU' use all that you have.
- `READ1` and `READ2` The ORP requires PE reads, sorry. Assembling with SE reads is not really worth it.. The reads can be gzipped. Always safe(r) to use the full path to the reads.
- `RUNOUT` Name the output. This name will be included in the final assembly, so choose wisely.
- appending `--dry-run` to the end of the command will print out to your screen the commands that will be run on your dataset, but they won't actually be run.   


How many resources do you need?

- **RAM** The amount of RAM you need scales with the number of unique kmers. The number of unique kmers is positively correlated with number of reads, which is a far more accessible number. In general, you need about 1GB RAM per million reads.
- **CPUs** More is better, but machines with between 24 and 64 cores are common. Note the assembly process is not able to leverage MPI.

##### Evaluating content using BUSCO

So you have an assembly, now how good is it? One way is to look at the assembly content. Are all the expected genes present? Fo this part of the exercise, We're going to use a 'real' assembly.

**Note:** The ORP runs BUSCO automatically, and has done so for your dummy assembly. See `$HOME/assembly_practical/reports/run*orthomerged/short_summary*txt`. A real assembly would score **much** better, hopefully...

```
mkdir $HOME/assembly_practical/assembly_eval && cd $HOME/assembly_practical/assembly_eval

python $(which run_BUSCO.py) -c 2 \
-m transcriptome \
-i $HOME/share/Day04/SRR1221220.orthomerged.fasta \
--out SRR1221220_ORPtest_YOURNAME
```

Unpacking
- `$(which run_BUSCO.py)` This is a trick to call a script using it's full path, without having to type the full path out yourself. Type `which run_BUSCO.py` and see what it returns..
- `-m transcriptome` you want to run BUSCO in transcriptome mode. If you were trying to analyze a genome, then your type `-m genome`
- `-i $HOME/share/SRR1221220.orthomerged.fasta` this is the input
- `--out SRR1221220` the output prefix.
- `-c 2` the number of CPUs to use to the analysis. Use as many as you have!


```
C:0.3%[S:0.0%,D:0.3%],F:0.3%,M:99.4%,n:303

1	Complete BUSCOs (C)
0	Complete and single-copy BUSCOs (S)
1	Complete and duplicated BUSCOs (D)
1	Fragmented BUSCOs (F)
301	Missing BUSCOs (M)
303	Total BUSCO groups searched
```

##### Evaluating assembly structure using TransRate

```
cd $HOME/assembly_practical/assembly_eval

transrate -o transrate_ORPtest_YOURNAME -t 2 \
-a $HOME/assembly_practical/assemblies/ORPtest_YOURNAME.orthomerged.fasta \
--left $HOME/assembly_practical/rcorr/ORPtest_YOURNAME.TRIM_1P.cor.fq \
--right $HOME/assembly_practical/rcorr/ORPtest_YOURNAME.TRIM_2P.cor.fq
```
Unpacking

- `-t 2` Two threads. Use more for a real assembly.
- `-a` The input assembly
- `--left` and `--right` The reads used to make the assembly.

Note: The ORP runs TransRate automatically, and has done so for your dummy assembly. See `$HOME/assembly_practical/reports/transrate*/assemblies.csv`

This little code snippet will make it easier to view this file.

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


```

cd $HOME/assembly_practical/assembly_eval/

salmon index --no-version-check \
--type quasi \
-k 31 \
-i ORPtest_YOURNAME.ortho.idx \
-t $HOME/assembly_practical/assemblies/ORPtest_YOURNAME.orthomerged.fasta

salmon quant --no-version-check -p 2 \
-i ORPtest_YOURNAME.ortho.idx \
--seqBias --gcBias -l a \
-1 $HOME/assembly_practical/rcorr/ORPtest_YOURNAME.TRIM_1P.cor.fq \
-2 $HOME/assembly_practical/rcorr/ORPtest_YOURNAME.TRIM_2P.cor.fq \
-o $HOME/assembly_practical/assembly_eval/salmon_orthomerged_ORPtest_YOURNAME
```

Note: The ORP runs BUSCO automatically, and has done so for your dummy assembly. See `$HOME/assembly_practical/assembly_eval/salmon_orthomerged_ORPtest_YOURNAME/quant.sf` A real assembly would score much better, hopefully...


```
Name    Length  EffectiveLength TPM     NumReads
Single_29       800     665.243 15070.675944    1080.000000
Single_23       519     381.486 15963.011041    656.000000
Single_25       1141    872.521 9075.326798     853.000000
Single_103      237     82.000  1970.671194     17.407591
Single_101      350     188.824 2310.625863     47.000000
Single_55       366     224.992 13368.014888    324.000000
Single_56       443     268.606 13506.878854    390.823744
Single_50       936     786.806 6177.772999     523.612910
Single_53       1355    1235.517        7926.712663     1055.000000
Shannon_ORPtest.shannon_cremaining1_62_0        3036    3096.499        7700.570713     2568.647157
```

#### Annotation using dammit (time permitting)
See http://www.camillescott.org/dammit/installing.html

##### installing (this has been done for you)

```
sudo apt-get update
sudo apt-get install python-pip python-dev python-numpy
pip install --upgrade pip
sudo pip install -U setuptools
sudo pip install dammit
```

##### Running the Annotation (untested and potentially broken)

```
mkdir $HOME/busco_dbs

dammit databases --database-dir $HOME/dammit_dbs \
--install --busco-group eukaryota

dammit annotate $HOME/assembly_practical/assemblies/ORPtest_YOURNAME.orthomerged.fa \
--busco-group eukaryota \
--n_threads 4 \
--database-dir $HOME/dammit_dbs/
--full

```



#### Bibliography

1. Oyster River Protocol: https://www.biorxiv.org/content/early/2017/11/22/177253
2. BUSCO: doi: 10.1093/molbev/msx319
3. TrasRate: doi: 10.1101/gr.196469.115
4. Salmon: doi: 10.1038/nmeth.4197
5. dammit: http://www.camillescott.org/dammit/
