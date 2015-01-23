Quantification of Representative Sequences (QRS) pipeline v. 1.1.0
------------------------------------------------------------------

AUTHOR: Enrique Gonzalez-Tortuero (gonzalez@igb-berlin.de)

NAME: Quantification of Representative Sequence (QRS.py)

DESCRIPTION: QRS is a script (originally written in Perl 5.14.2 but now in Python 3.4 with access to R using RPy2 module) that allows analysing NGS datasets automatically to study population structure. Currently, this pipeline is optimized for pyrosequencing platform sequences but it is possible to use it for another NGS technologies like Ion Torrent. This program may operate either in batch processing mode where all sequences are automatically analysed in an unsupervised way or it may interact with the user at various checkpoints if no parameters have been specified prior to execution.

USAGE:

The program QRS can be executed as batch command or as an interactive program if no parameters are defined. To understand how QRS works in batch mode, we provide several use-case scenarios.

	Creating a reference HMM file (-RM)
	-----------------------------------

# Basic parameters

When this option is used, the batch command can be executed using only basic parameters:
		QRS.py -RM -folder=test -reffiles=test.fasta

In this case, we call QRS to create a reference HMM file of test.fasta in the folder test using default parameters (the aligner is PRANK and QRS will call ReformAl to fine-tune the alignment).
However, QRS can create more than a single HMM file if you specify different reference files joined with a single dash (-) in the -reffiles parameter, like in the following example, where the program will create two HMM files using default parameters:
		QRS.py -RM -folder=test -reffiles=test1.fasta-test2.fasta

# Alignment options

Another option that you can modify is the multiple sequence alignment program. You can choose the program with the -aligner parameter. PRANK is the default aligner program, like in this example where QRS will create a single reference HMM file:
		QRS.py -RM -folder=test -aligner=prank -reffiles=test.fasta

To employ a different aligner for the sequence alignment task simply modify the -aligner parameter. In the following example, QRS will call MUSCLE to align a reference file to create a single reference HMM file:
		QRS.py -RM -folder=test -aligner=muscle -reffiles=test.fasta

Here, we use MAFFT with an iterative global alignment algorithm called G-INS-i to perform an accurate alignment:
		QRS.py -RM -folder=test -aligner=mafft-ginsi -reffiles=test.fasta

The last option you may modify is whether to use ReformAl to fine-tune the alignment or not. By default, this option is enabled (-reformal=yes), like in the following example, where QRS will create two reference HMM files after calling GramAlign to align your reference sequences:
		QRS.py -RM -folder=test -aligner=gramalign -reformal=yes -reffiles=test1.fasta-test2.fasta

If you want to disable this option, you have to type -reformal=no. In this example, QRS will align your three reference files using the default aligner (i.e., PRANK) and will not curate the alignments before creating three reference HMM files:
		QRS.py -RM -folder=test -reformal=no -reffiles=test1.fasta-test2.fasta-test3.fasta

In the following example, we call QRS to align four reference files using MAFFT with their iterative local alignment algorithm (L-INS-i) and we do not want to fine-tune the alignments:
		QRS.py -RM -folder=test -aligner=mafft-linsi -reformal=no -reffiles=test1.fasta-test2.fasta-test3.fasta-test4.fasta

Finally, if you are working with standard alignments (like SILVA [REFERENCE] or GreenGenes [REFERENCE]), perhaps you want to avoid the alignment step. To do it, you can use the parameter "-noalign" like in the following example, where we want to create a HMM file from SILVA alignment (in fasta format):
		QRS.py -RM -folder=test -noalign -reffiles=silva.fasta

	Describing the 454 data set (-AM1)
	----------------------------------

# Basic parameters

When this option is used, the batch command cannot be executed using only basic parameters. Instead, you have to specify if you have a 454 FASTA file or a 454 FASTQ file. In this example the input file is a FASTA file: 
	QRS.py -AM1 -folder=test -informat=fasta -fasta=file.fna

Here, the input file is a FASTQ file:
	QRS.py -AM1 -folder=test -informat=fastq -fastq=file.fastq

# Another input files
If you have a FASTA file as input, you can add a quality file to do all basic statistics on the base quality using the parameter -quality:
	QRS.py -AM1 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual

If you have a paired FASTQ file, you have to add the parameter -paired=yes and the name of the paired FASTQ:
	QRS.py -AM1 -folder=test -informat=fastq -fastq=file.fastq -paired=yes -fastq2=file2.fastq

	Processing a 454 data set (-AM2)
	--------------------------------

# Basic parameters

Before executing this part, you have to evaluate the characteristics of your data set according to the previous step (see AM1 section). As all 454 pyrosequencing data sets are different depending on the case study, there are no default parameters to filter and trim sequences according to length, quality and complexity. 
However, in a hypothetical case that you do not need to filter your data set, the input parameters are the input files added as describe before (see AM1 section). Moreover, you have to add the reference HMM file (created in -RM step, see RM section for more details), oligos and design files (read below to know more about these files) and specify if you have reverse barcoded primers or not. In the following example, we use a 454 FASTA file with quality file and our data set has reverse barcoded primers (-bry): 
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -bry

The oligos file is a plain text file that has two columns splitted by tabs (except if the label is barcode, where there are a third column). The first column contains always a label that indicates forward or reverse adapter, forward or reverse primer and barcode. The second contains the DNA sequence for each element and, in the case of barcode label, the third one indicates the barcode ID (MID1, MID2, MID3...). An example of this file is as follows:
		   seqadapfor	SEQUENCE
		   seqadaprev	SEQUENCE
		   forward	SEQUENCE
		   reverse	SEQUENCE
		   barcode	SEQUENCE	BARID

The design file is also a plain text file that has two or three columns (depending of the use of reverse barcoded primers) splitted by tabs. An example of these file is as follows:
		   BARID1	(BARID2)	SAMPLEID

Here, BARID1 is the Barcode ID for the forward primer, BARID2 is the Barcode ID for the reverse primer (if exists) and SAMPLEID is the name of the sample.

# Filtering by length

As it is said in AM2.1 section, there are no default parameters to filter and trim sequences according to length, quality and complexity because it depends of your data set. In the following example calls QRS to accept sequences that are greater than 300 bp in a 454 FASTQ file and have no reverse barcoded primers:
	QRS.py -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -minlen=300 -oligos=oligos.csv -design=design.csv -brn

In this example, QRS is executed to accept sequences that have between 355 and 500 bp (the data is given as 454 FASTA file without quality file):
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -hmmfile=test.hmm -minlen=355 -maxlen=500 -oligos=oligos.csv -design=design.csv -bry

# Filtering by GC content

Another possible filter is based on the GC content. In the following example QRS is called to accept sequences that have more than 60% GC content in a 454 FASTQ file and have no reverse barcoded primers:
	QRS.py -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -mingc=60 -oligos=oligos.csv -design=design.csv -brn

In this example, QRS is executed to accept sequences that have between 40 and 50% GC content (the data is given as 454 FASTA file without quality file):
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -hmmfile=test.hmm -mingc=40 -maxgc=50 -oligos=oligos.csv -design=design.csv -bry

# Filtering and trimming according to quality

We offer another filter based on base quality. In this example, QRS accepts only sequences that have at least a mean sequence quality value of 25:
	QRS.py -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -minqual=25 -oligos=oligos.csv -design=design.csv -bry

If you have a FASTA file, it is convenient to have a quality file to filter by base quality. In the following example, QRS will filter an input 454 FASTA file by length (accepting only sequences that are between 375 and 480 bp long) and quality (rejecting sequences that have a mean sequence quality score smaller than 28). This data set has no reverse barcoded primers:
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -minlen=375 -maxlen=480 -minqual=28 -brn

The following option is used to trim bases that have an insufficient quality score. In the following example, QRS will trim all nucleotides with a quality value less than 25 in the beginning and at the end of the sequences in a non-reverse barcoded 454 data set:
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -trimqualleft=25 -trimqualright=25 -brn

Finally, another way to filter your dataset according to quality is allowing ambiguous characters (Ns). The number of allowed Ns in your dataset can be modified with the parameter -allowns. In the following example, QRS accepts all sequences that have at most one ambiguous character in the sequence:
	QRS.py -AM2 -folder=folder -informat=fastq -fastq=file.fastq -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowns=1 -bry

However, it is not recommended to allow ambiguous characters as they might indicate bad quality nucleotides (Huse et al. 2007). In the following example, QRS removes all sequences that have more than one ambiguous character:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowns=0 -bry

# Filtering by complexity

Finally, if you want to remove all low complexity sequences like homopolymers, you can do this with-filmet and -filthr parameters. The first argument defines the method to remove this kind of sequences (DUST (Morgulis et al. 2006) or Entropy-based filter) and the second argument defines the threshold for the filtering method. For more details on these filtering methods, see PrinSeq manual (http://prinseq.sourceforge.net/manual.html#QCCOMPLEXITY). In the following example, we use DUST to remove all sequences with a complexity greater than 7:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -filmet=dust -filthr=7 -bry

In this example, we use entropy to remove all sequences with a complexity smaller than 70:
	QRS.py -AM2 -folder=folder -informat=FASTA -FASTA=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -filmet=entropy -filthr=70 -bry

In this example, a 454 data set (that consists in a FASTQ file) is filtered by length (accepting only sequences that are between 390 and 500 bp long), base quality (removing all sequences with mean sequence quality score less than 28), and low-complexity based on the DUST filter (considering all sequences with values greater than 5 as low complexity sequences):
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minlen=390 -maxlen=500 -minqual=28 -filmet=dust -filthr=5 -bry

Finally, it is feasible to trim poly-A/Ts tails. In this example, QRS trims the regions that have more than three followed A/T at the beginning or the end of the sequences:
	QRS.py -AM2 -folder=test -informat=fastq -FASTQ=file.fastq -hmmfile=test.hmm -trimtails=3 -oligos=oligos.csv -design=design.csv -bry

# Classifying sequences as a specific marker

You can also modify the maximum allowed e-value to classify a sequence as a specific marker using a HMM profile (-hmmthr). By default, this argument is set to 10-10 as then similar results as with the usage of BLAST are achieved (Altschul et al. 1990). We do not recommend changing this value. However, you can change the value like in the following example where QRS classifies a sequence as a specific marker if the e-value is smaller than 10-25:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -hmmthr=1E-25 -bry

# Assigning sequences to samples

If you want to modify the maximum number of allowed mismatches to detect the barcodes in your data set, you can use the parameter -allowmis. In this example, QRS classifies a non reverse barcoded data set in different samples with the maximum percentage of mismatches from 0.00 to 1.00. By defaulf, it is considered 0.01 and, although we do not recommend modifying this value, it is possible to change it like in this example, where QRS considers a maximum percentage of 1% to detect the barcodes in a reverse barcoded data set:
	QRS.py -AM2 -folder=folder -informat=FASTQ -FASTQ=file.FASTQ -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowmis=0.01 -bry

# Clustering step

Another parameter you can modify is -cutoff. This parameter is used to cluster sequences in order to de-noise your data set, i. e., to remove all spurious nucleotides across all sequences. By default, QRS calculates this value according to the CD-HIT-OTU algorithm (Li et al. 2012) but you can define the value by a number between 0.00 and 1.00. For example, if you want to cluster all sequences that have 99.5% of similarity, you can type:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -cutoff=0.995 -brn

The parameter -minclustersize specifies the minimum size of cluster to consider that analysed sequences are not sequencing artifacts. By default, this filter is enabled and all clusters that have at least three sequences are considered as good clusters. If you want to disable this argument, you have to type -minclustersize=0, like in this example:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minclustersize=0 -bry

# Alignment step

The last options you may modify concern the use of the alignment program and ReformAl to fine-tune the alignment. For more information, see the RM examples.

# A complex example

In the following example, QRS retrieves all sequences that have more than 350 bp, a mean sequence quality score of 30, no ambiguous characters, and a sequence complexity greater than 79 according to the entropy-based filter. This data set has reverse paired barcodes and the script discards all clusters that have only one sequence. Finally, QRS uses KAlign to align all accepted sequences and this alignment is not fine tuned:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minlen=350 -minqual=30 -allowns=0 -filmet=entropy -filthr=79 -bry -minclusterzise=2 -aligner=kalign -reformal=no

	Classifying and quantifying all 454 sequences into representative sequence variants according to Statistical Parsimony (-AM3)
	-----------------------------------------------------------------------------------------------------------------------------

# Basic parameters

When this option is used, the batch command can be executed with only basic parameters, like in the following example:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -derep=derep.list -names=names.list -sample=sample.csv -outfile=allsamples

In the previous example, QRS will use myaligneddata alignment to run TCS for the data and then retrieve a TCS-types FASTA file called allsamples.fasta and a frequencies matrix file called allsamples.abs.csv. The program makes use of some information from the samples.csv file. This file which is a plain text file that contains the following information always splitted by tabs:

                SAMPLEID	DATA1	(...)	DATAN

Here, SAMPLEID is the name of the sample and DATA are data you want to put here (Place, year, species...).

derep.list and names.list are files produced in the -AM2 step during the clustering step.

# Statistical parsimony parameters

You can modify the limit threshold that you want to use to determine representative sequences variants with the -limitid parameter. By default, QRS considers that two sequences belong to the same TCS-type if they have more similarity than 0.95, like in the following example:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -limitid=0.95 -outfile=allsamples

In this example, the limit threshold to determine if two sequences belong from the same TCS-type is 0.995:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -limitid=0.995 -outfile=allsamples

Another parameter you can modify is the fact that you want to consider gaps as another character or not (-considergaps). By default, it is enabled (-considergaps=1):
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -limitid=0.95 -considergaps=1 -outfile=allsamples

However, if you don't want to consider gaps, disable such parameter using -considergaps=0:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -limitid=0.95 -considergaps=0 -outfile=allsamples

Finally, the last parameter you want to modify about the Statistical parsimony algorithm is related with the transition/transversion relation (-transition=). As Templeton et al (1992) indicated in their paper, this parameter has only three possible values: 1 (extreme bias to transitions in the DNA - Kimura-2-Parameter model), 2 or 3 (no differences between transitions and transversion - Jukes-Cantor model). By default, this parameter is set as 1:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -limitid=0.95 -considergaps=1 -transition=1 -outfile=allsamples

Although we don't recommend to modify this value, if you consider than the probability of transition is the same than probability of transversion (i.e., Jukes-Cantor model), you should modify this value as follows:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -limitid=0.95 -considergaps=1 -transition=3 -outfile=allsamples

REQUIREMENTS:

Before using this pipeline, the following Perl modules and programs should be installed:

* Python modules:
	- Biopython (Bio module; Cock et al. 2009)
	- Numerical Python (numpy module; Oliphant 2007)
	- Processes and System information Utility (psutil module)
	- R from Python (rpy2 module)
	- Scientific Python (scipy module; Oliphant 2007)

* R packages:
	- ape (Paradis et al. 2004)
	- seqinr (Charif & Lobry 2007)

* Programs:
	- PrinSeq (Schmieder & Edwards 2011): it is used to generate summary statistics of sequence and quality data and to filter, reformat and trim 454 data. This program is publicly available at http://prinseq.sourceforge.net under the GPLv3 licence.
	- PEAR (Zhang et al. 2014): it is used to merge all paired fastq files after filtering all your sequences according to length, GC content, quality and complexity with PrinSeq. This program is publicly available at http://www.exelixis-lab.org/web/software/pear under the GPLv3 licence.
	- CutAdapt (Martin 2011): a Python script that removes primers and adapters. This program is publicly available at http://code.google.com/p/cutadapt under the MIT licence.
	- HMMER 3.1b (Finn et al 2011): it is used to search sequences using probabilistic models called profile hidden Markov models (profile HMM). This program is publicly available at http://hmmer.janelia.org/ under GPLv3 licence.
	- USEARCH (Edgar 2010): it is used to cluster sequences and detect chimaeras according to the UCHIME algorithm (Edgar et al. 2011). This program is available after requiring it according to the authors' page (http://www.drive5.com/usearch/).
	- At least one of the following aligners:
		* Clustal Omega (Sievers et al. 2011): a hybrid aligner that mixes progressive multiple sequence alignments with the use of Markov models. This program is publicly available at http://www.clustal.org/omega/ under GPLv2 licence.
		* FSA (Bradley et al. 2009): a probabilistic aligner that uses the sequence annealing technique (Schwartz & Pachter 2007) for constructing a multiple alignment from pairwise homology estimation. FSA is publicly available at http://fsa.sourceforge.net/ under GPL licence.
		* GramAlign (Russell et al. 2008): a time-efficient progressive aligner which estimates distances according to the natural grammar present in nucleotide sequences. This program is freely available at http://bioinfo.unl.edu/gramalign.php
		* KAlign (Lassmann & Sonnhammer 2005): a progressive aligner based on Wu-Manber string-matching algorithm. KAlign is publicly available at http://msa.sbc.su.se under GPLv2 licence.
		* MAFFT (Katoh & Standley 2013): an iterative/progressive alignment program that is publicly available at http://mafft.cbrc.jp/alignment/software/ under BSD licence.
		* MUSCLE (Edgar 2004): an iterative aligner that it is freely available at http://www.drive5.com/muscle/.
		* PicXAA (Sahraeian & Yoon 2010): a probabilistic non-progressive alignment algorithm that finds multiple sequence alignments with maximum expected accuracy that it is publicly available at http://gsp.tamu.edu/picxaa under GPLv3 licence.
		* PRANK (Loytynoja & Goldman 2005): a probabilistic multiple alignment program based on maximum likelihood methods used in phylogenetics that considers evolutionary distances between sequences. This program is publicly available at http://code.google.com/p/prank-msa under GPLv3 licence.

Although the aforementioned aligners have already been tested with the QRS, other aligners may be added to use in this pipeline. If you want to work with another alignment program, please feel free to contact with the author (qrspipeline@gmail.com) to include it in the source code of QRS. 
Finally, an optional recommended requirement might be the installation of ReformAlign (Lyras and Metzler submitted). ReformAlign is a recently proposed profile-based meta-alignment approach that aims to fine-tune existing alignments via the employment of standard profiles.

HISTORY: 

v 1.1.0 - Major issues:
          * Script translated to Python 3 with access to R (thanks to RPy2 project)
          * Script infers representative sequence variants according to Statistical Parsimony (Templeton et al. 1992) without executing TCS. This important modification generates the following changes:
                      * Removed the "Declustering" step at the end of -AM2 step: In the previous version, the QRS pipeline declustered all sequences after the alignment step because TCS requires a non-clusterized dataset in order to run Statistical Parsimony. As we implemented such algorithm in our program without executing TCS, this step is unnecessary and the output from this step is, directly, a clustered alignment file in order to apply such algorithm. For this reason, the FASTA header must have the following structure: ">(Sample_Name)_(NumberID);size=(Total number of sequences that are identical)".
                      * Added "-limitid", "-considergaps" and "-transition" parameters. All of them are parameter related to Statistical Parsimony according to the Templeton et al. (1992) algorithm. For more info about how to use them, read all -AM3 step examples.
                      * Requiral of DEREP and NAMES files in -AM3 step: such files are produced after clustering sequences in USEARCH and created in QRS pipeline as plain text files as "derep.list" and "names.list". These files contains two columns separated by a tabulation. The first column have the name for the sequence and the second column have all sequences that are considered the same (in the case of derep.list are 100%% identical, while in names.list are synonims due to the denoising step). In batch mode, the parameters are "-derep=(Namefile)" and "-names=(Namefile)", respectively.
          * Added a new process in the -AM2 step between the filtering sequences by length, quality, GC content and complexity using PrinSeq and the filtering sequences by a specific marker using HMMER processes. This new process is called "Merging sequences using PEAR" and it is only useful if you are running paired fastq files. PEAR is a required program in order to run this pipeline (see REQUIREMENTS).

	  Minor issues:
          * Added "-noalign" parameter in the -RM step. It is only useful if you are working with standard alignments like SILVA or GreenGenes.
          * Added a new aligner PicXAA (Sahraeian & Yoon 2010). In order to use it, you can choose two options: "-aligner=picxaa-pf" (if you use partition function) or "-aligner=picxaa-phmm" (if you use pair HMM in order to evaluate the alignment).
          * Divided "-trimqual" parameter (to trim nucleotides according to quality) as "-trimqualleft=(Number)" and "-trimqualright=(Number)".
          * Fixed minor error about colon in FASTA files after processing in PrinSeq.
          * "-design=" and "-br[y|n]" are only mandatory if you want to split your dataset into samples. If not, you MUST to write "-nosplit" parameter in the batch mode.

v 1.0.0 - Original version as a Perl script.

REFERENCES:

	* Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) Basic local alignment search tool. Journal of Molecular Biology 215: 403-10.
	* Bradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. PLoS Computational Biology 5: e1000392.
	* Charif D, Lobry JR (2007) SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. IN Bastolla U, Porto M, Roman HE, Vendruscolo M. Structural approaches to sequence evolution: Molecules networks and analysis. Springer-Verlag, New York (USA), 207-32.
	* Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11): 1422-23.
	* Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5): 1792-7.
	* Edgar RC (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19): 2460-1.
	* Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27(16): 2194-200.
	* Finn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39: W29-W37.
	* Huse SM, Huber JA, Morrison HG, Sogin ML, Welch DM (2007) Accuracy and quality of massively parallel DNA pyrosequencing. Genome Biology 8: R143.
	* Katoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33: 511-18.
	* Katoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research 30: 3059-66.
	* Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4): 772-80.
	* Lassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. BMC Bioinformatics 6: 298.
	* Li W, Fu L, Niu B, Wu S, Wooley J (2012) Ultrafast clustering algorithms for metagenomic sequence analysis. Briefings in Bioinformatics 13(6): 656-68.
	* Loytynoja A, Goldman N (2005) An algorithm for progressive multiple alignment of sequences with insertions. Proceedings of the National Academy of Sciences of USA 102: 10557-62.
	* Lyras DP, Metzler (Submitted) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach". BMC Bioinformatics. 
	* Martin M (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17.
	* Morgulis A, Gertz EM, Schaffer AA, Agarwala R (2006) A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13: 1028-40.
	* Oliphant TE (2007) Python for scientific computing. Computing in Science & Engineering 9: 10-20.
	* Paradis E, Claude J, Strimmer K (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290. 
	* Russell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. BMC Bioinformatics 9: 306.
	* Sahraeian SME, Yoon BJ (2010) PicXAA: greedy probabilistic construction of maximum expected accuracy alignment of multiple sequences. Nucleic Acids Research 38(15): 4917-28.
	* Schmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. Bioinformatics 27: 863-4.
	* Schwartz AS, Pachter L (2007) Multiple alignment by sequence annealing. Bioinformatics 23: e24-9.
	* Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Soeding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7: 539.
	* Templeton AR, Crandall KA, Sing CF (1992) A cladistic analysis of phenotypic associations with haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram estimation. Genetics 132: 619-33.
	* Zhang J, Kobert K, Flouri T, Stamatakis A (2014) PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30(5): 614-20.
