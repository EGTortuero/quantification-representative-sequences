# quantification-representative-sequences
Automatically exported from code.google.com/p/quantification-representative-sequences

QRS is a script (originally written in Perl 5.14.2, now in Python 3.4 with access to R using RPy2 module) that allows analyzing amplicon sequencing datasets automatically to study population structure. Currently, this pipeline is optimized for pyrosequencing platform sequences but it is possible to use it for another NGS technologies like Ion Torrent. This program may operate either in batch processing mode where all sequences are automatically analyzed in an unsupervised way or it may interact with the user at various checkpoints if no parameters have been specified prior to execution.

When using this pipeline, you must to cite their use:

González-Tortuero, E., Rusek, J., Petrusek, A., Gießler, S., Lyras, D., Grath, S., Castro-Monzón, F. and Wolinska, J. (2015), The Quantification of Representative Sequences pipeline for amplicon sequencing: case study on within-population ITS1 sequence variation in a microparasite infecting Daphnia. Molecular Ecology Resources. doi: 10.1111/1755-0998.12396

Source code (for all versions of the program) is now under GitHub repositories.

HISTORY OF THE SOURCE CODE:

v 1.4.0 - Major issues:

    Modified all parameters to be used with USEARCH 8+. In this program, several parameters were changed in order to cluster sequences.
    Added the parameter "-noalign" in "-AM2" step. It is only useful if you want to launch the "-AM3" step using NJ clustering.
    Added the new parameter "-namesample" that it is mandatory when you have FASTA/FASTQ files for different samples. In this case, this name is used in order to define sequences in each sample and promoting a correct analysis in the "Inferring Representative Sequences" step. The name of the sample should be as is in the first column of the SAMPLES file.
    Changed the way of indicating the number of allowed mismatches in order to detect the barcoded tags into the sequences. Now it is referred as the number of base pairs that could be indicate a mismatch in the barcode detection instead of a percentage. Additionally, this step is optional only if you want to assign sequences to samples.
        When the obtained dataset in "-RM" and "-AM2" steps has more than a thousand sequences, QRS can modify the behaviour of the aligners in order to deal properly with this dataset according to the aligners manuals. 
    Added optional MAFFT parameters in order to deal with huge datasets: "-aligner=mafft-dpparttree", "-aligner=mafft-fastaparttree" and "-aligner=mafft-parttree". All of them are described in Katoh & Toh (2007) and they are recommended for more than 10,000 sequences according to the authors.
    Added a new optional "aligner" that it is the most accurate way to align huge datasets at this moment ("-aligner=muscle-profile"). This method divides the dataset in two: a small subset of ten sequences that belongs to the most abundant sequences ("core") and the rest of the sequences ("rest"). The core subset is aligned using muscle with default parameters and, then, every sequence in the rest subset is aligned against the core using a profile alignment. This method is well recommended when you have several thousands of sequences due to the accuracy of this method (Sievers et al. 2013). However, this method is really slow and it can take more than a day.
    Simplified code in order to detect chimeras and de-noising sequences according to USEARCH v. 7.0.1003+ manual.
    Added new option in the parameter "-filtermethod=" (in "-AM2" step) that allows to execute HECTOR (Wirawan et al. 2014) in order to correct homopolymers in 454 reads. 

    Minor issues: 

    Added Opal (Wheeler & Kececioglu 2007) as a new aligner ("-aligner=opal").
    Added the number of threads parameter in PEAR program.
    In the assigning sequences to samples, you may know in this step hoy many sequences are in each sample after removing primers.
    Modified CLUSTAL Omega parameters in order to deal properly fastafiles in "-RM" and "-AM2" step.
    Fixed small bug about IDs that are no present in DEREP in NAMES files but are in the fasta file when "-AM3" step is executed.
    Fixed small bug about IDs after clustering sequences when a cluster is marked as a chimera.
    Fixed small bug about alignment options in "-AM2" step (it was disabled by error). 

v 1.3.0 - Current stable version.

    Major issues: 

    Fixed a bug in calculation of absolute frequencies using NJ. Now, the log file shows the total number sequences for each network while the absolute matrix should show the number of sequences for each network in every sample.
    Fixed a bug in the number of allowed mismatches to detect primers and number of allowed ambiguities (Ns). Now, the parameters "-allowedmis" and "-allowns" should be work properly.
    Fixed an error about displaying the help in the program using batch mode. Now, the option "--help" should be work properly.
    Fixed an error in Statistical Parsimony method on the "non-parsimony" probability correction. As this step could underestimate the number of networks, we employed the proposed solution of Templeton et al. (1992) to collapse all non-corrected networks only in a single step because it is compute-expensive (and not in all potential steps as they originally proposed). This fix generates the most similar results than TCS (Clement et al. 2000) 

v 1.2.0 - Major issues:

    Added a new parameter "-method" in order to pool your sequences according to Neighbour Joining ("-method=NJ") or Statistical Parsimony ("-method=SP"). 

    Fixed a bug in calculation of absolute frequencies matrix using SP in the last version. Now, the log file shows the total number sequences for each network while the absolute matrix should show the number of sequences for each network in every sample.
    Fixed an error in Statistical Parsimony method on the "non-parsimony" probability correction. This step is made after determining the most representative sequences for each network considering only parsimony probability. 

    Minor issues:

        Added in the output of interactive batch mode the equivalence to the command line options. 

    Added the verbose output in ReformAlign?
    Fixed bugs in command line options in -RM options.
    Fixed bugs in number of CPUs to run HMMER in -RM options.
    Updated references about ReformAlign? 

v 1.1.0 - Major issues:

    Script translated to Python 3 with access to R (thanks to RPy2 project)
    Script infers representative sequence variants according to Statistical Parsimony (Templeton et al. 1992) without executing TCS. This important modification generates the following changes:
        Removed "Declustering" step at the end of -AM2 step: In the previous version, the QRS pipeline declustered all sequences after the alignment step because TCS requires a non-clusterized dataset in order to run Statistical Parsimony. As we implemented such algorithm in our program without executing TCS, this step is unnecessary and the output from this step is, directly, a clustered alignment file in order to apply such algorithm. For this reason, the FASTA header must have the following structure: ">(Sample_Name?)(NumberID);size=(Total number of sequences that are identical)".
        Added "-limitid", "-considergaps" and "-transition" parameters. All of them are parameter related to Statistical Parsimony according to the Templeton et al. (1992) algorithm. For more info about how to use them, read all -AM3 step examples.
        Requiral of DEREP and NAMES files in -AM3 step: such files are produced after clustering sequences in USEARCH and created in QRS pipeline as plain text files as "derep.list" and "names.list". These files contains two columns separated by a tabulation. The first column have the name for the sequence and the second column have all sequences that are considered the same (in the case of derep.list are 100%% identical, while in names.list are synonims due to the denoising step). In batch mode, the parameters are "-derep=(Namefile)" and "-names=(Namefile)", respectively. 
    Added a new process in the -AM2 step between the filtering sequences by length, quality, GC content and complexity using PrinSeq? and the filtering sequences by a specific marker using HMMER processes. This new process is called "Merging sequences using PEAR" and it is only useful if you are running paired fastq files. PEAR is a required program in order to run this pipeline (see REQUIREMENTS). 

    Minor issues:

        Added "-noalign" parameter in the -RM step. It is only useful if you are working with standard alignments like SILVA or GreenGenes?.
        Added a new aligner PicXAA (Sahraeian & Yoon 2010). In order to use it, you can choose two options: "-aligner=picxaa-pf" (if you use partition function) or "-aligner=picxaa-phmm" (if you use pair HMM in order to evaluate the alignment).
        Divided "-trimqual" parameter (to trim nucleotides according to quality) as "-trimqualleft=(Number)" and "-trimqualright=(Number)".
        Fixed minor error about colon in FASTA files after processing in PrinSeq?.
        "-design=" and "-br[y|n]" are only mandatory if you want to split your dataset into samples. If not, you MUST to write "-nosplit" parameter in the batch mode. 

v 1.0.0 - Original version as a Perl script. Originally used in Gonzalez-Tortuero et al (2015) Molecular Ecology Resources to perform all comparisons between 454 amplicon sequencing and cloning/Sanger sequencing datasets 
