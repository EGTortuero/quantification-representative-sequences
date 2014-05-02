#! /usr/bin/perl

## Loading PERL libraries
# Core modules
use Cwd;
use strict;
use warnings;

# CPAN modules
use Bio::AlignIO;
use Bio::SeqIO;
use File::Find::Rule;
use File::Which;
use String::Approx 'amatch';
use Sys::CPU;
use Sys::MemInfo qw (totalmem freemem);

## Starting computer-related variables
my $SO = $^O;
my $ncpu = Sys::CPU::cpu_count();
my $totalmemory = int(&totalmem / 1024**2);
my $freememory = int(&freemem / 1024**2);

# Starting global variables of the script
my $option = "";
my $directory = "";
my $nrparameters = "";

## Print header of the program
print "		#########################################		\n";
print "		#	       QRS v. 1.0       	#		\n";
print "		#---------------------------------------#		\n";
print "		#   Author: Enrique Gonzalez-Tortuero   #		\n";
print "		#########################################		\n\n";
print "QRS v.1.0  Copyright (C) 2014 Enrique Gonzalez-Tortuero\n";
print "This program comes with ABSOLUTELY NO WARRANTY; for details\ntype 'y' in the interface or '--GNU3' in batch mode.\n";
print "This is free software, and you are welcome to redistribute it\nunder certain conditions; see above for details.\n\n";
print "Running in $SO with $totalmemory Mb RAM ($freememory Mb free) and $ncpu processors.\n";

## Implementing several programs to the pipeline
# Processing NGS data files programs
my $prinseqprogram1 = which('prinseq-lite.pl');
my $prinseqprogram2 = which('prinseq-graphs-noPCA.pl');

# Removing primers program
my $cutadaptprogram = which('cutadapt');

# HMMER program
my $hmmerprogram = which('hmmbuild');
my $hmmerprogram3 = which('nhmmer');

# Clustering program (important to denoising and removing chimaeras step).
my $usearchprogram = which('usearch');

# Aligners programs
my $clustaloprogram = which('clustalo');
my $fsaprogram = which('fsa');
my $gramalignprogram = which('GramAlign');
my $kalignprogram = which('kalign');
my $mafftprogram = which('mafft');
my $muscleprogram = which('muscle');
my $prankprogram = which('prank');

# Meta-aligner program
my $reformalprogram = which('ReformAlign');

# Haplotyping programs
## TCS program
my $TCSfile;
my $TCSrule = File::Find::Rule->new;
$TCSrule->file;
$TCSrule->name("TCS1.21.jar");
my @file = $TCSrule->in($ENV{HOME});
my $listfiles = scalar (@file);
if ($listfiles >= 1) {
	$TCSfile = $file[0];
} else {
	HelpingU();
}
undef @file;

# Another programs
my $java = which('java');
my $perldoc = which('perldoc');

## The first parameter is, always, the step you want to run. If there isn't any argument, run interactive mode.
my $parameter = $ARGV[0];
if (!defined($parameter) or $parameter eq '') {
	print "\nIn which folder do you wan to run all analyses?\n";
	$directory = <STDIN>;
	chomp $directory;
	print "\n\nWhat step do you want to run?\n";
	print "a) Aligning and creating HMM for your reference file/s.\n";
	print "b) Calculating the basic statistics for your NGS dataset.\n";
	print "c) Processing all your NGS dataset to obtain an alignment.\n";
	print "d) Processing your alignment to establish TCS-types\n";
	print "y) GNU license version 3\n";
	print "z) Help\n";
	print "Please, enter a letter: ";
	my $letteroption = <STDIN>;
	chomp $letteroption;
	if ($letteroption eq 'a' or $letteroption eq 'A' or $letteroption eq '1') {
		$option = "-RM";
	} elsif ($letteroption eq 'b' or $letteroption eq 'B' or $letteroption eq '2') {
		$option = "-AM1";
	} elsif ($letteroption eq 'c' or $letteroption eq 'C' or $letteroption eq '3') {
		$option = "-AM2";
	} elsif ($letteroption eq 'd' or $letteroption eq 'D' or $letteroption eq '4') {
		$option = "-AM3";
	} elsif ($letteroption eq 'y' or $letteroption eq 'Y') {
		$option = "--GNU3";
	} elsif ($letteroption eq 'z' or $letteroption eq 'Z') {
		$option = "--help";
	} else {
		die "The requested option isn't exists.\n";
	}
} else {
	$nrparameters = scalar (@ARGV);
	for (my $i = 0; $i < $nrparameters; $i++) {
		if ($ARGV[$i] =~ /^-folder\=(\S+)/) {
			$directory = $1;
		}
		if ($ARGV[$i] =~ /^(-[A|R]M\d?)/ or $ARGV[$i] =~ /^-h$/ or $ARGV[$i] =~ /^--help$/ or $ARGV[$i] =~ /^--GNU3$/) {
			$option = $1;
		}
	}
}

# Checking variables
if ($directory eq "" or $option eq "") {
	HelpingU(); 
}

# Displaying help
if ($option eq '-h' or $option eq '--help') {
	HelpingU();
}

# Displaying GNUv3 license
if ($option eq '--GNU3') {
	GNUv3();
}

# Going to the specified directory
my $direcname = getcwd;
if ($directory eq ".") {
	chdir("$direcname") or die "It impossible to go to $direcname folder!\n";
} else {
	chdir("$direcname\/$directory") or die "It impossible to go to $direcname/$directory folder!\n";
}

# Making all analyses
if ($option eq "-RM") {

	my @references = ();
	my @alignedreferences = ();
	my $alignmentfile = "";
	my $alnreffile = "";
	my $aligner = "";
	my $nralignments = "";
	my $reformalswitch = "";
	
	# Batch mode or interactive mode?
	my $parameter2 = $ARGV[3];
	if (!defined($parameter2) or $parameter2 eq '') {
 		print "\nWhat aligner do you want to use (clustalo, fsa, gramalign, kalign, muscle, mafft-fast, mafft-ginsi, mafft-linsi or prank)? By default, it is prank.\n";
		$aligner = <STDIN>;
		chomp $aligner;
		print "";
		print "\nDo you want to use ReformAl to fine-tune the alignment (i.e., as a meta-aligner)?\nBy default, it is yes.\n";
		my $metalign = <STDIN>;
		chomp $metalign;
		if ($metalign eq "" or $metalign eq 'yes' or $metalign eq 'y' or $metalign eq 'Yes' or $metalign eq 'Y') {
			$reformalswitch = 1;
		} elsif ($metalign eq 'no' or $metalign eq 'n' or $metalign eq 'No' or $metalign eq 'N') {
			$reformalswitch = 0;
		} else {
			die "I don't recognize the parameter.\n";
		}
		print "\nWhat is/are the name/s of the reference fasta file/s?\n";
		print "If there are many reference fasta files, please, put the names of all files separated by a dash (-)\n";
		my $refnam = <STDIN>;
		chomp $refnam;
		@references = split(/-/,$refnam);
	} else {
		for (my $i = 3; $i < $nrparameters; $i++) {
			if ($ARGV[$i] =~ /^-aligner\=(\S+)/) {
				$aligner = $1;
			}
			if ($ARGV[$i] =~ /^-reformal\=(\S+)/) {
				my $metalign = $1;
				if ($metalign eq 'yes' or $metalign eq 'y' or $metalign eq 'Yes' or $metalign eq 'Y') {
					$reformalswitch = 1;
				} elsif ($metalign eq 'no' or $metalign eq 'n' or $metalign eq 'No' or $metalign eq 'N') {
					$reformalswitch = 0;
				}
			}
			if ($ARGV[$i] =~ /^-reffiles\=(\S+)/) {
				@references = split(/-/,$1);
			}

		}
	}
	
	$nralignments = scalar(@references);
	if ($reformalswitch eq "") {
		$reformalswitch = 1;
	}
	
	print "\n###############################";
	print "\n# Creating all reference HMM. #";
	print "\n###############################\n";

	#Running the step
	foreach my $referencefile (@references) {
		my $alignedreferencefile = $referencefile;
		my $metalnreffile;
		my $hmmreferencefile;
		if ($referencefile =~ /\.fas$/) {	
			$alignedreferencefile =~ s/\.fas$/\.aligned/;
			$alignmentfile = $alignedreferencefile . ".best.fas";
			$alnreffile = $alignmentfile;
			$alnreffile =~ s/\.best\.fas$/\.fas/;
			if ($reformalswitch == 1) {
				$metalnreffile = $alnreffile;
				$metalnreffile =~ s/\.fas$/\.reformal\.fas/;
				$hmmreferencefile = $metalnreffile;
				$hmmreferencefile =~ s/\.aligned\.reformal\.fas$/\.hmm/;
			} else {
				$hmmreferencefile =~ s/\.aligned\.fas$/\.hmm/;
			}
		} elsif ($referencefile =~ /\.fna$/) {
			$alignedreferencefile =~ s/\.fna$/\.aligned/;
			$alignmentfile = $alignedreferencefile . ".best.fas";
			$alnreffile = $alignmentfile;
			$alnreffile =~ s/\.best\.fas$/\.fna/;
			if ($reformalswitch == 1) {
				$metalnreffile = $alnreffile;
				$metalnreffile =~ s/\.fna$/\.reformal\.fna/;
				$hmmreferencefile = $metalnreffile;
				$hmmreferencefile =~ s/\.aligned\.reformal\.fna$/\.hmm/;
			} else {
				$hmmreferencefile =~ s/\.aligned\.fna$/\.hmm/;
			}
		} elsif ($referencefile =~ /\.fasta$/) {
			$alignedreferencefile =~ s/\.fasta$/\.aligned/;
			$alignmentfile = $alignedreferencefile . ".best.fas";
			$alnreffile = $alignmentfile;
			$alnreffile =~ s/\.best\.fas$/\.fasta/;
			if ($reformalswitch == 1) {
				$metalnreffile = $alnreffile;
				$metalnreffile =~ s/\.fasta$/\.reformal\.fasta/;
				$hmmreferencefile = $metalnreffile;
				$hmmreferencefile =~ s/\.aligned\.reformal\.fasta$/\.hmm/;
			} else {
				$hmmreferencefile =~ s/\.aligned\.fasta$/\.hmm/;
			}
		} else {
			HelpingU();
		}

		print "\nAligning reference file: $referencefile\n";
		if ($aligner eq "clustalo") {
			system("$clustaloprogram -i $referencefile -o $alnreffile -v --threads=$ncpu");
		} elsif ($aligner eq "fsa") {
			system("$fsaprogram $referencefile > $alnreffile");
		} elsif ($aligner eq "gramalign") {
			system("$gramalignprogram -i $referencefile -o $alnreffile -f 2");
		} elsif ($aligner eq "kalign") {
			system("$kalignprogram -i $referencefile -o $alnreffile -f fasta");
		} elsif ($aligner =~ /^mafft/) {
			if ($aligner eq "mafft-ginsi") {
				system("$mafftprogram --maxiterate 1000 --globalpair --ep 0.123 --thread $ncpu $referencefile > $alnreffile");  
			} elsif ($aligner eq "mafft-linsi") {
				system("$mafftprogram --maxiterate 1000 --localpair --ep 0.123 --thread $ncpu $referencefile > $alnreffile");  
			} elsif ($aligner eq "mafft-fast") {
				system("$mafftprogram --ep 0.123 --thread $ncpu $referencefile > $alnreffile");
			}
		} elsif ($aligner eq "muscle") {
			system("$muscleprogram -in $referencefile -out $alnreffile");
		} elsif ($aligner eq "prank" or $aligner eq "") {
			system("$prankprogram -d=$referencefile -o=$alignedreferencefile");
			Cat($alnreffile,$alignmentfile);
			unlink $alignmentfile;
		}
		
		if ($reformalswitch == 1) {
			print "\nImproving the original alignment using ReformAl.\n";
			system("$reformalprogram -i $referencefile -o $metalnreffile -a $alnreffile");
			print "\nCreating the HMM for your alignment/s.\n";
			system("$hmmerprogram $hmmreferencefile $metalnreffile");
		} else {
			print "\nCreating the HMM for your alignment/s.\n";
			system("$hmmerprogram $hmmreferencefile $alnreffile");
		}

		# Cleaning all intermediate files
		unlink $alnreffile;
		if ($reformalswitch == 1) {
			unlink $metalnreffile;
		}
	}

	print "\nIt is done!\nAll alignments were made in ";
	if ($aligner eq "clustalo") {
		print "Clustal Omega (Sievers et al 2011).";
	} elsif ($aligner eq "fsa") {
		print "FSA (Bradley et al 2009).";
	} elsif ($aligner eq "gramalign") {
		print "GramAlign (Russell et al 2008).";
	} elsif ($aligner eq "kalign") {
		print "KAlign (Lassmann & Sonnhammer 2005).";
	} elsif ($aligner =~ /^mafft/) {
		print "MAFFT (Katoh & Standley 2013).";
	} elsif ($aligner eq "muscle") {
		print "MUSCLE (Edgar 2004).";
	} elsif ($aligner eq "prank" or $aligner eq "") {
		print "PRANK (Löytynoja & Goldman 2008).";
	}
	if ($reformalswitch == 1) {
		print " Then, all alignments were fine-tuned using ReformAl (Lyras and Metzler 2014).";
	}
	print " Finally, HMM was made using HMMER (Finn et al 2011).\nCitation:\n";
	if ($aligner eq "fsa") {
		print "Bradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. PLoS Computational Biology. 5:e1000392.\n";
	} elsif ($aligner eq "muscle") {
		print "Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5):1792-7.\n";
	}
	print "Finn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39:W29-W37.\n";
	if ($aligner =~ /^mafft/) {
		print "Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n";
	} elsif ($aligner eq "kalign") {
		print "Lassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. BMC Bioinformatics 6:298\n";
	} elsif ($aligner eq "prank" or $aligner eq "") {
		print "Löytynoja A, Goldman N (2008) A model of evolution and structure for multiple sequence alignment. Philosophical Transactions of the Royal Society B 363(1512):3913-9.\n";
	}
	if ($reformalswitch == 1) {
		print "Lyras DP, Metzler D (2014) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach. Submitted\n";
	}
	if ($aligner eq "gramalign") {
		print "Russell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. BMC Bioinformatics 9:306.\n";
	} elsif ($aligner eq "clustalo") {
		print "Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539\n";
	}
		
} elsif ($option eq "-AM1") {

	# Starting variables
	my $informatfile = "";
	my $NGSfastafile = "";
	my $NGSqualityfile = "";
	my $NGSfastqfile = "";
	my $pairedswitch = "";
	my $NGSfastq2file = "";
	my $phred64switch = "";
	my $NGSgraphdatafile = "";
	my $NGSrootnamefile = "";
	
	my $parameter2 = $ARGV[3];
	if (!defined($parameter2) or $parameter2 eq '') {
		print "\nWhat kind of file(/s) do you have (fastq or fasta file)?\n";
		$informatfile = <STDIN>;
		chomp $informatfile;
		if ($informatfile eq 'fasta') {
			print "\nWhat is its name?\n";
			$NGSfastafile = <STDIN>;
			chomp $NGSfastafile;
			print "\nDo you have a quality file (yes/no)?\n";
			my $qualexists = <STDIN>;
			chomp $qualexists;
			if ($qualexists eq 'yes' or $qualexists eq 'y' or $qualexists eq 'Yes' or $qualexists eq 'Y') {
				print "\nWhat is its name?\n";
				$NGSqualityfile = <STDIN>;
				chomp $NGSqualityfile;
			}			
		} elsif ($informatfile eq 'fastq') {
			print "\nWhat is its name?\n";
			$NGSfastqfile = <STDIN>;
			chomp $NGSfastqfile;
			print "\nAre they paired or not?\n";
			my $paired = <STDIN>;
			chomp $paired;
			if ($paired eq 'yes' or $paired eq 'y' or $paired eq 'Yes' or $paired eq 'Y') {
				$pairedswitch = 1;
				print "\nWhat is its name?\n";
				$NGSfastq2file = <STDIN>;
				chomp $NGSfastq2file;
			} elsif ($paired eq 'no' or $paired eq 'No' or $paired eq 'n'  or $paired eq 'N') {
				$pairedswitch = 0;
			} else {
				die "\nI don't recognize your answer.\n";
			}
			print "\nHave your quality data got Phred+64 format (yes/no)?\n";
			my $phred64 = <STDIN>;
			chomp $phred64;
			if ($phred64 eq 'yes' or $phred64 eq 'Yes' or $phred64 eq 'y' or $phred64 eq 'Y') {
				$phred64switch = 1;
			} elsif ($phred64 eq 'no' or $phred64 eq 'No' or $phred64 eq 'n' or $phred64 eq 'N') {
				$phred64switch = 0;
			} else {
 				die "\nI don't recognize your answer\n";
 			}
		}
	} else {
		for (my $i = 3; $i < $nrparameters; $i++) {
			if ($ARGV[$i] =~ /^-informat\=(\S+)/) {
				$informatfile = $1;
			}
			if ($informatfile eq 'fasta') {	
				if ($ARGV[$i] =~ /^-fasta\=(\S+)/) {
					$NGSfastafile = $1;
				}
				if ($ARGV[$i] =~ /^-quality\=(\S+)/) {
					$NGSqualityfile = $1;
				}
			} elsif ($informatfile eq 'fastq') {
				if ($ARGV[$i] =~ /^-fastq\=(\S+)/) {
					$NGSfastqfile = $1;
				}
				if ($ARGV[$i] =~ /^-paired\=(\S+)/) {
					my $paired = $1;
					if ($paired eq 'yes' or $paired eq 'Yes' or $paired eq 'y' or $paired eq 'Y') {
						$pairedswitch = 1;
						if ($ARGV[$i] =~ /^-fastq2\=(\S+)/) {
							$NGSfastq2file = $1;
						}
					} elsif ($paired eq 'no' or $paired eq 'No' or $paired eq 'n' or $paired eq 'N') {
						$pairedswitch = 0;
					}
				}
				if ($ARGV[$i] =~ /^-phred64\=(\w+)/) {
					my $phred64 = $1;
					if ($phred64 eq 'yes' or $phred64 eq 'Yes' or $phred64 eq 'y' or $phred64 eq 'Y') {
						$phred64switch = 1;
					} elsif ($phred64 eq 'no' or $phred64 eq 'No' or $phred64 eq 'n' or $phred64 eq 'N') {
						$phred64switch = 0;
					}
				}
			}
		}
	}
	
	if ($pairedswitch eq "") {
		$pairedswitch = 0;
	}
	if ($phred64switch eq "") {
		$phred64switch = 0;
	}
	
	# Printing all basic statistics from NGS data
	if ($informatfile eq 'fasta') {
		if ($NGSfastafile =~ /\.fna$/) {
			$NGSrootnamefile = $NGSfastafile;
			$NGSrootnamefile =~ s/\.fna$//;
		} elsif ($NGSfastafile =~ /\.fas$/) {
			$NGSrootnamefile = $NGSfastafile;
			$NGSrootnamefile =~ s/\.fas$//;		
		} elsif ($NGSfastafile =~ /\.fasta$/) {
			$NGSrootnamefile = $NGSfastafile;
			$NGSrootnamefile =~ s/\.fasta$//;		
		}
		$NGSgraphdatafile = $NGSrootnamefile . ".gd";
	} elsif ($informatfile eq 'fastq') {
		$NGSgraphdatafile = $NGSfastqfile;
		$NGSgraphdatafile =~ s/\.fastq/\.gd/;
		$NGSrootnamefile = $NGSfastqfile;
		$NGSrootnamefile =~ s/\.fastq//;
	}
			
	print "\n########################################################################################";
	print "\n# Obtaining the basic statistics on the NGS dataset you have and the different graphs. #";
	print "\n########################################################################################\n";
	my $allparameters;
	if ($informatfile eq 'fasta') {
		if ($NGSfastafile ne '') {
			my $inputparameter1 = " -fasta '$NGSfastafile'";
			$allparameters .= $inputparameter1;
		} else {
			HelpingU();
		}
		if ($NGSqualityfile ne '') {
			my $inputparameter2 = " -qual '$NGSqualityfile'";
			$allparameters .= $inputparameter2;
		}
	} elsif ($informatfile eq 'fastq') {
		if ($NGSfastqfile ne '') {
			my $inputparameter3 = " -fastq '$NGSfastqfile'";
			$allparameters .= $inputparameter3;
		} else {
			HelpingU();
		}
		if ($NGSfastq2file ne '') {
			my $inputparameter4 = " -fastq2 '$NGSfastq2file'";
			$allparameters .= $inputparameter4;
		}
		if ($phred64switch != 0) {
			my $phred64parameter = " -phred64";
			$allparameters .= $phred64parameter;
		}
	}
	if ($NGSgraphdatafile ne '') {
		my $outparameter = " -graph_data '$NGSgraphdatafile' -out_good null -out_bad null -verbose";
		$allparameters .= $outparameter;
	}

	system("$prinseqprogram1 $allparameters");
	system("$prinseqprogram2 -i $NGSgraphdatafile -html_all -o $NGSrootnamefile"); # Creating the HTML information page
	unless(-e "PRINSEQ_GRAPHS" or mkdir "PRINSEQ_GRAPHS") {
		die "Unable to create PRINSEQ_GRAPHS in your current folder.\n";
	}
	system("$prinseqprogram2 -i $NGSgraphdatafile -png_all -o PRINSEQ_GRAPHS/$NGSrootnamefile"); # Creating all PNG files
	print "\nIt is done! You have the basic statistics in the HTML file $NGSrootnamefile.html. This information was analysed using PrinSeq (Schmieder & Edwards 2011).\nCitation:\nSchmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. Bioinformatics 27:863-4.\n"; 

} elsif ($option eq "-AM2") {

 	# Starting variables
 	my $informatfile = "";
 	my $NGSfastafile = "";
 	my $NGSqualityfile = "";
 	my $NGSfastqfile = "";
 	my $pairedswitch = "";
 	my $NGSfastq2file = "";
 	my $phred64switch = "";
 	my $minimumlength = "";
 	my $maximumlength = "";
 	my $minimumgc = "";
 	my $maximumgc = "";
 	my $minimumqual = "";
 	my $allowedns = "";
 	my $trimmingtails = "";
 	my $trimmingqual = "";	
 	my $filtermethod = "";
	my $filterthreshold = "";
	my $hmmreferencefile = "";
	my $nosplit = "";
	my $oligosfile = "";
	my $barcodedrevprim = "";
	my $designfile = "";
	my $allowedmismatches = "";
	my $hmmthr = "";
	my $samplefile = "";
	my $outname = "";
	my $cutoff = "";
	my $aligner = "";
	my $reformalswitch = "";
	my $nohaplswitch = "";
	my $minclustersize = "";
 
	# Batch mode or interactive mode?
	my $parameter2 = $ARGV[3];
	if (!defined($parameter2) or $parameter2 eq '') {
	print "\nFirst step: knowing all parameters to filtering and trimming sequences according to quality.\n";
	print "Be sure you see before the basic statistics.\n";
	print "What kind of file(/s) do you have (fastq or fasta file)?\n";
  		$informatfile = <STDIN>;
  		chomp $informatfile;
  		if ($informatfile eq 'fasta') {
  			print "\nWhat is its name?\n";
  			$NGSfastafile = <STDIN>;
  			chomp $NGSfastafile;
  			print "\nDo you have a quality file (yes/no)?\n";
  			my $qualexists = <STDIN>;
  			chomp $qualexists;
  			if ($qualexists eq 'yes' or $qualexists eq 'y' or $qualexists eq 'Yes' or $qualexists eq 'Y') {
  				print "\nWhat is its name?\n";
  				$NGSqualityfile = <STDIN>;
 				chomp $NGSqualityfile;
				print "\nDo you want to filter sequences according to a minimum quality?\n";
				my $answer1 = <STDIN>;
 				chomp $answer1;
				if ($answer1 eq 'yes' or $answer1 eq 'Yes' or $answer1 eq 'y' or $answer1 eq 'Y') { 
					print "\nWhat is the minimum average quality of the sequences to pass the filters?\n";
					$minimumqual = <STDIN>;
					chomp $minimumqual;
 				}
 				print "\nDo you want to trim sequences according to a quality threshold?\n";
				my $answer2 = <STDIN>;
 				chomp $answer2;
				if ($answer2 eq 'yes' or $answer2 eq 'Yes' or $answer2 eq 'y' or $answer2 eq 'Y') {
					print "\nHow many quality do you consider to remove bad nucleotides?\n";
					$trimmingqual = <STDIN>;
					chomp $trimmingqual;
				}
 			}			
		} elsif ($informatfile eq 'fastq') {
			print "\nWhat is its name?\n";
			$NGSfastqfile = <STDIN>;
			chomp $NGSfastqfile;
			print "\nAre they paired or not?\n";
			my $paired = <STDIN>;
			chomp $paired;
			if ($paired eq 'yes' or $paired eq 'y' or $paired eq 'Yes' or $paired eq 'Y') {
				$pairedswitch = 1;
				print "\nWhat is its name?\n";
				$NGSfastq2file = <STDIN>;
				chomp $NGSfastq2file;
			} elsif ($paired eq 'no' or $paired eq 'n' or $paired eq 'No' or $paired eq 'N') {
				$pairedswitch = 0;
			}
			print "\nHave your quality data got Phred+64 format (yes/no)?\n";
			my $phred64 = <STDIN>;
			chomp $phred64;
			if ($phred64 eq 'yes' or $phred64 eq 'y' or $phred64 eq 'Yes' or $phred64 eq 'Y') {
				$phred64switch = 1;
			} else {
				$phred64switch = 0;
			}
			print "\nDo you want to filter sequences according to a minimum quality (yes/no)?\n";
			my $answer1 = <STDIN>;
 			chomp $answer1;
			if ($answer1 eq 'yes' or $answer1 eq 'Yes' or $answer1 eq 'y' or $answer1 eq 'Y') { 
				print "\nWhat is the minimum average quality of the sequences to pass the filters?\n";
				$minimumqual = <STDIN>;
				chomp $minimumqual;
 			}
 			print "\nDo you want to trim sequences according to a quality threshold?\n";
			my $answer2 = <STDIN>;
 			chomp $answer2;
			if ($answer2 eq 'yes' or $answer2 eq 'Yes' or $answer2 eq 'y' or $answer2 eq 'Y') {
				print "How many quality do you consider to remove bad nucleotides?\n";
				$trimmingqual = <STDIN>;
				chomp $trimmingqual;
			}
		}
		print "\nDo you want to filter sequences according to minimum length (yes/no)?\n";
		my $answer3 = <STDIN>;
		chomp $answer3;
		if ($answer3 eq 'yes' or $answer3 eq 'y' or $answer3 eq 'Yes' or $answer3 eq 'Y') {
			print "\nWhat is the sequences minimum length (in bp) to pass the filters?\n";
			$minimumlength = <STDIN>;
			chomp $minimumlength;
		}
		print "\nDo you want to filter sequences according to maximum length (yes/no)?\n";
		my $answer3b = <STDIN>;
		chomp $answer3b;
		if ($answer3b eq 'yes' or $answer3b eq 'y' or $answer3b eq 'Yes' or $answer3b eq 'Y') {
			print "\nWhat is the sequences maximum length (in bp) to pass the filters?\n";
			$maximumlength = <STDIN>;
			chomp $maximumlength;
		}
		print "\nDo you want to filter sequences according to minimum GC content (yes/no)?\n";
		my $answer3c = <STDIN>;
		chomp $answer3c;
		if ($answer3c eq 'yes' or $answer3c eq 'y' or $answer3c eq 'Yes' or $answer3c eq 'Y') {
			print "\nWhat is the sequences minimum GC content (in percentage) to pass the filters?\nPlease, don't write the '%' symbol.\n";
			$minimumgc = <STDIN>;
			chomp $minimumgc;
		}
		print "\nDo you want to filter sequences according to maximum GC content (yes/no)?\n";
		my $answer3d = <STDIN>;
		chomp $answer3d;
		if ($answer3d eq 'yes' or $answer3d eq 'y' or $answer3d eq 'Yes' or $answer3d eq 'Y') {
			print "\nWhat is the sequences maximum GC content (in percentage) to pass the filters?\nPlease, don't write the '%' symbol.\n";
			$maximumgc = <STDIN>;
			chomp $maximumgc;
		}
		print "\nHow many uncertain nucleotides (Ns) do you allow in your sequences?\n";
		$allowedns = <STDIN>;
		chomp $allowedns;
		print "\nHow many homopolynucleotide tails (eg., poly-As) do you allow before trimming sequences in their extremes?\n";
		$trimmingtails = <STDIN>;
		chomp $trimmingtails;
		print "\nDo you want to remove homopolymers (yes/no)?\n";
		my $answer4 = <STDIN>;
		chomp $answer4;
		if ($answer4 eq 'yes' or $answer4 eq 'y' or $answer4 eq 'Yes' or $answer4 eq 'Y') {
			print "\nWhich 'homopolymers removal' method do you want to use (entropy or dust)?\n";
			$filtermethod = <STDIN>;
			chomp $filtermethod;
			print "\nWhich threshold do you want to use with $filtermethod?\n";
			$filterthreshold = <STDIN>;
			chomp $filterthreshold;
		}
		
 		print "\nSecond step: knowing the HMM filename to filtering sequences according to a previous reference.\n";
 		print "\nWhat is the reference HMM file or database name?\n";
 		$hmmreferencefile = <STDIN>;
 		chomp $hmmreferencefile;
 		print "\nWhat is the maximum allowed e-value (by default: 1e-10)?\n";
 		$hmmthr = <STDIN>;
 		chomp $hmmthr;
		print "\nThird step: knowing all parameters and the OLIGOS and DESIGN files to splitting sequences according to the samples and remove primers.\n";
 		print "\nDo you want to split sequences according to your samples (yes/no)?\n";
 		my $answer30 = <STDIN>;
 		chomp $answer30;
		if ($answer30 eq 'yes' or $answer30 eq 'Yes' or $answer30 eq 'y' or $answer30 eq 'Y') {
			$nosplit = 0;
			print "\nThis program needs the information of the primers. It is in a file called OLIGOS with has the following structure:\nseqadapfor\tSEQUENCE\nseqadaprev\tSEQUENCE\nforward\tSEQUENCE\nreverse\tSEQUENCE\nbarcode\tSEQUENCE\tBARID\n\nWhat is the OLIGOS file name?\n";
			$oligosfile = <STDIN>;
			chomp $oligosfile;
			print "\nDo you have barcoded reverse primers? (Yes or No)\n";
			my $brp = <STDIN>;
			chomp $brp;
			if ($brp eq 'yes' or $brp eq 'Yes' or $brp eq 'y' or $brp eq 'Y') {
				$barcodedrevprim = '-bry';
			} elsif ($brp eq 'no' or $brp eq 'No' or $brp eq 'n' or $brp eq 'N') {
				$barcodedrevprim = '-brn';
			} else {
				die "Unrecognized answer\n";
			}
			print "\nAnother requered file is DESIGN. It has a table with two or three columns (it depends if you has barcoded reverse primers or not) and as the following structure:\nBARID1\t(BARDID2)\tSAMPLEID\n\nWhat is the DESIGN file name?\n";
			$designfile = <STDIN>;
			chomp $designfile;
			print "\nHow many mismatches do you allow to detect primers in your sequences?\nIt can be an exact number or a percenteage. By default is 1%.\n";
			$allowedmismatches = <STDIN>;
			chomp $allowedmismatches;
		} else {
			$nosplit = 1;
		}
 		print "\nFourth step: knowing all parameters to denoising your dataset and filtering your dataset.\n";
 		print "\nWhich cutoff do you want to use to clustering sequences?\nThis number must be between 0.00 and 1.00! If you don't specify any number, the threshold could be calculated automatically.\n";
		$cutoff = <STDIN>;
		chomp $cutoff;
 		print "\nAnother QRS's filter is based on cluster's size. We consider that clusters with few sequences are 454 artifacts.\nWhich is the minimum cluster size do you want to allow? (By default, we consider three sequences as a good cluster size)\n";
		$minclustersize = <STDIN>;
		chomp $minclustersize;
 		print "\nWhat aligner do you want to use (clustalo, fsa, gramalign, kalign, muscle, mafft-fast, mafft-ginsi, mafft-linsi or prank)? By default, we are using prank.\n";
		$aligner = <STDIN>;
		chomp $aligner;
		print "\nDo you want to use ReformAl to fine-tune the alignment (i.e., as a meta-aligner)?\nBy default, it is yes.\n";
		my $metalign = <STDIN>;
		chomp $metalign;
		if ($metalign eq "" or $metalign eq 'yes' or $metalign eq 'y' or $metalign eq 'Yes' or $metalign eq 'Y') {
			$reformalswitch = 1;
		} elsif ($metalign eq 'no' or $metalign eq 'n' or $metalign eq 'No' or $metalign eq 'N') {
			$reformalswitch = 0;
		} else {
			die "I don't recognize the parameter.\n";
		}
		
		
	} else {
		for (my $i = 3; $i < $nrparameters; $i++) {
			if ($ARGV[$i] =~ /^-informat\=(\S+)/) {
				$informatfile = $1;
			}
			if ($informatfile eq 'fasta') {	
				if ($ARGV[$i] =~ /^-fasta\=(\S+)/) {
					$NGSfastafile = $1;
				}
				if ($ARGV[$i] =~ /^-quality\=(\S+)/) {
					$NGSqualityfile = $1;
				}
			} elsif ($informatfile eq 'fastq') {
				if ($ARGV[$i] =~ /^-fastq\=(\S+)/) {
					$NGSfastqfile = $1;
				}
				if ($ARGV[$i] =~ /^-paired\=(\S+)/) {
					my $paired = $1;
					if ($paired eq 'yes' or $paired eq 'Yes' or $paired eq 'y') {
						$pairedswitch = 1;
						if ($ARGV[$i] =~ /^-fastq2\=(\S+)/) {
							$NGSfastq2file = $1;
						}
					} elsif ($paired eq 'no' or $paired eq 'No' or $paired eq 'n') {
						$pairedswitch = 0;
					}
				}
				if ($ARGV[$i] =~ /^-phred64\=(\w+)/) {
					my $phred64 = $1;
					if ($phred64 eq 'yes' or $phred64 eq 'Yes' or $phred64 eq 'y') {
						$phred64switch = 1;
					} elsif ($phred64 eq 'no' or $phred64 eq 'No' or $phred64 eq 'n') {
						$phred64switch = 0;
					}
				}
			}
			if ($ARGV[$i] =~ /^-minlen\=(\d+)/) {
				$minimumlength = $1;
			}
			if ($ARGV[$i] =~ /^-maxlen\=(\d+)/) {
				$maximumlength = $1;
			}
			if ($ARGV[$i] =~ /^-mingc\=(\d+)/) {
				$minimumgc = $1;
			}
			if ($ARGV[$i] =~ /^-maxgc\=(\d+)/) {
				$maximumgc = $1;
			}
			if ($ARGV[$i] =~ /^-minqual\=(\d+)/) {
				$minimumqual = $1;
			}
			if ($ARGV[$i] =~ /^-trimqual\=(\d+)/) {
				$trimmingqual = $1;
			}
			if ($ARGV[$i] =~ /^-trimtails\=(\d+)/) {
				$trimmingtails = $1;
			}
			if ($ARGV[$i] =~ /^-allowns\=(\d+)/) {
				$allowedns = $1;
			}
			if ($ARGV[$i] =~ /^-filmet\=(\w+)/) {
				$filtermethod = $1;
			}
			if ($ARGV[$i] =~ /^-filthr\=(\d+)/) {
				$filterthreshold = $1;
			}
			if ($ARGV[$i] =~ /^-hmmfile=(\S+)/) {
				$hmmreferencefile = $1; 		
			}
			if ($ARGV[$i] =~ /^-hmmthr\=(\S+)/) {
				$hmmthr = $1;
 			}
 			if ($ARGV[$i] =~ /^-nosplit$/) {
				$nosplit = 1;
 			}
 			if ($ARGV[$i] =~ /^-oligos=(\S+)/) {
				$oligosfile = $1;
			}
			if ($ARGV[$i] =~ /^-design=(\S+)/) {
				$designfile = $1;
			}
			if ($ARGV[$i] =~ /^(-br[y|n])/) {
				$barcodedrevprim = $1;				
			} 
			if ($ARGV[$i] =~ /^-allowmis\=(\S+)/) {
				$allowedmismatches = $1;
			}
			if ($ARGV[$i] =~ /^-cutoff\=(\d\.\d+)/) {
				$cutoff = $1;
			}
			if ($ARGV[$i] =~ /^-minclustersize\=(\d+)/) {
				$minclustersize = $1;
			}
			if ($ARGV[$i] =~ /^-aligner\=(\S+)/) {
				$aligner = $1;
			}
			if ($ARGV[$i] =~ /^-reformal\=(\S+)/) {
				my $metalign = $1;
				if ($metalign eq 'yes' or $metalign eq 'y' or $metalign eq 'Yes' or $metalign eq 'Y') {
					$reformalswitch = 1;
				} elsif ($metalign eq 'no' or $metalign eq 'n' or $metalign eq 'No' or $metalign eq 'N') {
					$reformalswitch = 0;
				}
			}
		}
	}

	# Defining variables
	my $NGSrootnamefile = "";
	my $NGSgoodsequences = "";
	my $NGSlogfile = "";
	my $NGSgoodfastafile = "";
	my $NGSgoodfastqfile = "";
	my $NGSrevcompfastafile = "";
	my $NGStwicefastafile = "";
	my $NGSgrcfileslist = "";
	my $NGScontigs = "";
	my $HMMresults = "";
	my $NGShmmfilteredfastafile = "";
	my $tempalltrimmedderepfile = "";
	my $alltrimmedfile = "";
	my $alltrimmeduniquetable = "";
	my $alltrimmedderepfile = "";
	my $alltrimmeddereptable = "";
	my $alltrimmeddereplist = "";
	my $alltrimmedclustfile = "";
	my $alltrimmedclustsortedfile = "";
	my $alltrimmedclustlist = "";
	my $alltrimmednonchimerafile = "";
	my $alltrimmedalignment = "";
	my $alltrimmeddeclustered = "";
	my $badsequences = "fakeseqs";
	my $metNGSfile = "";
	my $avelength = "";
	
	# Starting names of files
	if ($informatfile eq 'fasta') {
		if ($NGSfastafile =~ /\.fna$/) {
			$NGSrootnamefile = $NGSfastafile;
			$NGSrootnamefile =~ s/\.fna$//;
		} elsif ($NGSfastafile =~ /\.fas$/) {
			$NGSrootnamefile = $NGSfastafile;
			$NGSrootnamefile =~ s/\.fas$//;		
		} elsif ($NGSfastafile =~ /\.fasta$/) {
			$NGSrootnamefile = $NGSfastafile;
			$NGSrootnamefile =~ s/\.fasta$//;		
		}
	} elsif ($informatfile eq 'fastq') {
		$NGSrootnamefile = $NGSfastqfile;
		$NGSrootnamefile =~ s/\.fastq//;
	}
	$NGSgoodsequences = $NGSrootnamefile . ".good";
	$NGSlogfile = $NGSrootnamefile . ".txt";

	if ($pairedswitch eq "") {
		$pairedswitch = 0;
	}
	if ($phred64switch eq "") {
		$phred64switch = 0;
	}
	if ($hmmthr eq "") {
		$hmmthr = 1e-10;
	}
	if ($nosplit eq "") {
		$nosplit = 0;
	}
	if ($allowedmismatches eq "") {
		$allowedmismatches = "1%";
	}
	if ($cutoff ne "") {
		if ($cutoff > 1.00 or $cutoff < 0.00) {
			HelpingU();
		}
	}
	if ($minclustersize eq "") {
		$minclustersize = 3;
	}

	# Filtering and trimming sequences
	print "\nFiltering and trimming sequences from the NGS data set according to your parameters.";
	print "\n------------------------------------------------------------------------------------\n";
	my $allparameters;
	if ($informatfile eq 'fasta') {
		if ($NGSfastafile ne '') {
			my $inputparameter1 = " -fasta '$NGSfastafile'";
			$allparameters .= $inputparameter1;
		} else {
			HelpingU();
		}
		if ($NGSqualityfile ne '') {
			my $inputparameter2 = " -qual '$NGSqualityfile'";
			$allparameters .= $inputparameter2;
		}
	} elsif ($informatfile eq 'fastq') {
		if ($NGSfastqfile ne '') {
			my $inputparameter3 = " -fastq '$NGSfastqfile'";
			$allparameters .= $inputparameter3;
		} else {
			HelpingU();
		}
		if ($NGSfastq2file ne '') {
			my $inputparameter4 = " -fastq2 '$NGSfastq2file'";
			$allparameters .= $inputparameter4;
		}
		if ($phred64switch == 1) {
			my $phred64parameter = " -phred64";
			$allparameters .= $phred64parameter;
		}
	}	
	if ($minimumqual ne '') {
		my $qualparameter1 = " -min_qual_mean $minimumqual";
		$allparameters .= $qualparameter1;
	}
	if ($trimmingqual ne '') {
		my $qualparameter2 = " -trim_qual_left $trimmingqual -trim_qual_right $trimmingqual";
		$allparameters .= $qualparameter2;
	}
	if ($trimmingtails ne '') {
		my $trimtailparameter = " -trim_tail_left $trimmingtails -trim_tail_right $trimmingtails";
		$allparameters .= $trimtailparameter;
	}
	if ($allowedns ne '') {
		my $allownsparameter = " -ns_max_n $allowedns";
		$allparameters .= $allownsparameter;
	}
	if ($minimumlength ne '') {
		my $lengthparameter = " -min_len $minimumlength";
		$allparameters .= $lengthparameter;
	}
	if ($maximumlength ne '') {
		my $lengthparameter2 = " -max_len $maximumlength";
		$allparameters .= $lengthparameter2;
	}
	if ($minimumgc ne '') {
		my $lengthparameter = " -min_gc $minimumgc";
		$allparameters .= $lengthparameter;
	}
	if ($maximumgc ne '') {
		my $lengthparameter2 = " -max_gc $maximumgc";
		$allparameters .= $lengthparameter2;
	}
	if ($filtermethod ne '' and $filterthreshold ne '') {
		my $homopolymerparameter = " -lc_method $filtermethod -lc_threshold $filterthreshold";
		$allparameters .= $homopolymerparameter;
	}
	if ($NGSgoodsequences ne '') {
		my $outparameter = " -out_format 1 -out_good '$NGSgoodsequences' -out_bad '$badsequences' -log '$NGSlogfile' -seq_case upper -verbose";
		$allparameters .= $outparameter;
	}
	system("$prinseqprogram1 $allparameters");
	
	# Starting variables
	$NGSgoodfastafile = $NGSgoodsequences . ".fasta";
	$NGSrevcompfastafile = $NGSgoodsequences . ".rc.fasta";
	$NGStwicefastafile = $NGSgoodsequences . ".twice.fasta";
	$NGSgrcfileslist = join(',',$NGSgoodfastafile,$NGSrevcompfastafile);
	$HMMresults = $NGStwicefastafile;
	$HMMresults =~ s/\.twice\.fasta$/\.twice\.table/;
	$NGShmmfilteredfastafile = $NGSgoodsequences . ".hmmfiltered.fasta";
	$tempalltrimmedderepfile = $NGShmmfilteredfastafile;
	$tempalltrimmedderepfile =~ s/\.fasta$/\.temptrimmed\.fasta/;
	$alltrimmedfile = $NGShmmfilteredfastafile;
	$alltrimmedfile =~ s/\.fasta$/\.trimmed\.fasta/;
	$alltrimmedderepfile = $alltrimmedfile;
	$alltrimmedderepfile =~ s/\.fasta$/\.derep\.fasta/;
	$alltrimmeddereptable = $alltrimmedderepfile . ".uctbl";
	$alltrimmeddereptable =~ s/\.fasta//;
	$alltrimmeddereplist = $alltrimmeddereptable;
	$alltrimmeddereplist =~ s/\.uctbl$/\.list/;
	$alltrimmedclustfile = $alltrimmedderepfile;
	$alltrimmedclustfile =~ s/\.fasta$/\.clust\.fasta/;
	$alltrimmeduniquetable = $alltrimmedclustfile . ".uctbl";
	$alltrimmedclustlist = $alltrimmeduniquetable;
	$alltrimmedclustlist =~ s/\.uctbl$/\.list/;
	$alltrimmedclustsortedfile = $alltrimmedclustfile;
	$alltrimmedclustsortedfile =~ s/\.fasta$/\.sorted\.fasta/;
	$alltrimmednonchimerafile = $alltrimmedclustsortedfile;
	$alltrimmednonchimerafile =~ s/\.fasta$/\.nonchimera\.fasta/;
	$alltrimmedalignment = $alltrimmednonchimerafile;
	$alltrimmedalignment =~ s/\.fasta$/\.aligned\.fasta/;
	if ($reformalswitch == 1) {
		$metNGSfile = $alltrimmedalignment;
		$metNGSfile =~ s/\.fasta$/\.reformal\.fasta/;
		$alltrimmeddeclustered = $metNGSfile;
	} else {
		$alltrimmeddeclustered = $alltrimmedalignment;
	}
	$alltrimmeddeclustered =~ s/\.fasta$/\.declustered\.fasta/;
	$alltrimmeddeclustered =~ s/\.derep//;
	$alltrimmeddeclustered =~ s/\.clust//;
	$alltrimmeddeclustered =~ s/\.sorted//;
		
	# Creating the Reverse Complementary sequences and making a whole dataset to analyse...
	ModifyingHeaders($NGSgoodfastafile);
	RevComSeq($NGSgoodfastafile, $NGSrevcompfastafile);
	Cat($NGStwicefastafile, $NGSgrcfileslist);

	# Filtering sequences according to the references
	print "\nFiltering your NGS sequences according to your HMM reference file.";
	print "\n------------------------------------------------------------------\n";
	system("$hmmerprogram3 --toponly --dfamtblout $HMMresults --cpu $ncpu $hmmreferencefile $NGStwicefastafile > /dev/null");
	HMMfiltering_SingleHMM($HMMresults,$NGStwicefastafile,$hmmthr);
	my $nrheaders = CountingHeaders($NGShmmfilteredfastafile);
	print "In this step, you recovered $nrheaders sequences.\n";
		
	# Splitting samples and removing barcodes.
	if ($nosplit == 0) {
		print "\nPutting sequences into samples and removing primers.";
		print "\n----------------------------------------------------\n";
		if ($oligosfile eq "" or $designfile eq "" or $barcodedrevprim eq "") {
			HelpingU();
		} else {
			Splitting($NGShmmfilteredfastafile,$oligosfile,$designfile,$barcodedrevprim,$allowedmismatches);
		}
		my $fileslist = "";
		my @files = glob('*.trimmed.fasta');
		foreach my $trimmedfile (@files) {
			if ($fileslist ne "") {
				$fileslist .= ",$trimmedfile";
			} else {
				$fileslist = $trimmedfile;
			}
		}
		Cat($alltrimmedfile, $fileslist);
		my $nrheaders2 = CountingHeaders($alltrimmedfile);
		print "In this step, you recovered $nrheaders2 sequences.\n";
	} else {
		Cat($alltrimmedfile, $NGShmmfilteredfastafile); 
	}
	
	# Clustering and denoising sequences
	print "\nDereplicating all your sequences.";
	print "\n---------------------------------\n";
	system("$usearchprogram -derep_fulllength $alltrimmedfile -output $tempalltrimmedderepfile -uc $alltrimmeddereptable -sizeout");
	system("$usearchprogram -sortbysize $tempalltrimmedderepfile -output $alltrimmedderepfile -quiet");
	unlink $tempalltrimmedderepfile;
	USEARCHanalysis($alltrimmeddereptable);

	if ($cutoff eq "") {
		$avelength = AveLength($alltrimmedfile);
		$cutoff = ($avelength - 1.000) / $avelength - 0.001;
	}
	print "\nClustering all your sequences at $cutoff (denoising step).";
	print "\n--------------------------------------------------------------\n";
	system("$usearchprogram -cluster_smallmem $alltrimmedderepfile -usersort -id $cutoff -uc $alltrimmeduniquetable -centroids $alltrimmedclustfile -sizein -sizeout -qmask none -idprefix $ncpu");
	USEARCHanalysis($alltrimmeduniquetable);
	my $perid = sprintf("%.8f", $cutoff * 100);

	# Filtering clusters with few sequences
	print "\nFiltering clusters by abundance.";
	print "\n--------------------------------\n";
	system("$usearchprogram -sortbysize $alltrimmedclustfile -output $alltrimmedclustsortedfile -minsize $minclustersize");

	# Removing chimaeras
	print "\nRemoving chimaeric sequences from all your sequences.";
	print "\n-----------------------------------------------------\n";
	system("$usearchprogram -uchime_denovo $alltrimmedclustsortedfile -nonchimeras $alltrimmednonchimerafile");

	#Aligning sequences...
	print "\nAligning all your sequences.";
	print "\n----------------------------\n";
	if ($aligner eq "kalign") {
		system("$kalignprogram -i $alltrimmednonchimerafile -o $alltrimmedalignment -f fasta");
	} elsif ($aligner eq "clustalo") {
		system("$clustaloprogram -i $alltrimmednonchimerafile -o $alltrimmedalignment -v --threads=$ncpu");
	} elsif ($aligner eq "fsa") {
		system("$fsaprogram --fast $alltrimmednonchimerafile > $alltrimmedalignment");
	} elsif ($aligner eq "gramalign") {
		system("$gramalignprogram -i $alltrimmednonchimerafile -o $alltrimmedalignment -f 2");
	} elsif ($aligner =~ /^mafft/) {
		if ($aligner eq "mafft-fast") {
 			system("$mafftprogram --ep 0.123 --thread $ncpu $alltrimmednonchimerafile > $alltrimmedalignment");
		} elsif ($aligner eq "mafft-linsi") {
			system("$mafftprogram --maxiterate 1000 --localpair --ep 0.123 --thread $ncpu $alltrimmednonchimerafile > $alltrimmedalignment");
		} elsif ($aligner eq "mafft-ginsi") {
			system("$mafftprogram --maxiterate 1000 --globalpair --ep 0.123 --thread $ncpu $alltrimmednonchimerafile > $alltrimmedalignment");
		}
	} elsif ($aligner eq "muscle") {
		system("$muscleprogram -in $alltrimmednonchimerafile -out $alltrimmedalignment");
	} elsif ($aligner eq "prank") {
		my $intermediatealignment = $alltrimmedalignment . ".best.fas";
		system("$prankprogram -d=$alltrimmednonchimerafile -o=$alltrimmedalignment");
		Cat($alltrimmedalignment,$intermediatealignment);
		unlink $intermediatealignment;
	}
	
	if ($reformalswitch == 1) {
		print "\nImproving the original alignment using ReformAl.\n";
		system("$reformalprogram -i $alltrimmednonchimerafile -o $metNGSfile -a $alltrimmedalignment");
		DeclusteringFASTAfile($metNGSfile,$alltrimmedclustlist,$alltrimmeddereplist);
		
	} else {
		DeclusteringFASTAfile($alltrimmedalignment,$alltrimmedclustlist,$alltrimmeddereplist);
	}
	
	# Removing intermediate files
	unlink glob("${badsequences}*");
	unlink glob('*.trimmed.fasta');
 	unlink glob('*.tbl');
 	unlink glob('*.allblastedtbl');
	unlink glob('*.list');
 	unlink $HMMresults;
 	unlink $NGShmmfilteredfastafile;
 	unlink $NGStwicefastafile;
 	unlink $NGSgoodfastafile;
 	unlink $NGSrevcompfastafile;
 	unlink $NGShmmfilteredfastafile;
	unlink $alltrimmeddereptable;
	unlink $alltrimmedderepfile;
 	unlink $alltrimmeduniquetable;
 	unlink $alltrimmedclustfile;
 	unlink $alltrimmednonchimerafile;
 	unlink $alltrimmedalignment;
	unlink $alltrimmedclustsortedfile;
 	if ($reformalswitch == 1) {
		unlink $metNGSfile
	}
	
	# Printing final results
	print "\nIt is done! Your data was analysed according to your starting parameters. ";
	if ($avelength ne "") {
		print "Moreover, all your sequences were analysed based on CDHIT-OTU (Li et al 2012) but with the following modifications.\n";
	}
	print "Your data was processed in PrinSeq (Schmieder & Edwards 2011) to filter and trimming according to length and quality. Then, HMMER (Finn et al 2011) was launched in order to detect your requested sequences. Later, this sequences were splitted according to the OLIGOS and DESIGN files and the primers were removed using CutAdapt (Martin 2011).\nAll your fasta files were clustered at $perid % ID (allowing to cluster sequences that are very similar and have 1-2 mismatches) using USEARCH (Edgar 2010). Then, chimaeric sequences were removed using UCHIME (Edgar et al. 2011). After that, all accepted sequences were aligned using ";
	if ($aligner eq "clustalo") {
		print "Clustal Omega (Sievers et al 2011)";
	} elsif ($aligner eq "fsa") {
		print "FSA (Bradley et al 2008)\n";
	} elsif ($aligner eq "gramalign") {
		print "GramAlign (Russell et al 2008)";
	} elsif ($aligner eq "kalign") {
		print "KAlign (Lassmann & Sonnhammer 2005)";
	} elsif ($aligner =~ /^mafft/) {
		print "MAFFT (Katoh & Standley 2013).";
	} elsif ($aligner eq "muscle") {
		print "MUSCLE (Edgar 2004)";
	} elsif ($aligner eq "prank" or $aligner eq "") {
		print "PRANK (Löytynoja & Goldman 2008)";
	}
	if ($reformalswitch == 1) {
		print " and automatically fine-tuned using ReformAl (Lyras and Metzler 2014)";
	}
	print ".\nCitation:\n";
	if ($aligner eq "fsa") {
		print "Bradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. PLoS Computational Biology. 5:e1000392.\n";
	}
	print "Finn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39:W29-W37.\n";
	if ($aligner eq "muscle") {
		print "Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5):1792-7.\n";
	}
	print "Edgar RC (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-1.\n";
	print "Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27(16):2194-200.\n";
	if ($aligner =~ /^mafft/) {
		print "Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n";
	} elsif ($aligner eq "kalign") {
		print "Lassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. BMC Bioinformatics 6:298\n";
	}
	if ($avelength ne "") {
	print "Li W, Fu L, Niu B, Wu S, Wooley J (2012) Ultrafast clustering algorithms for metagenomic sequence analysis. Briefings in Bioinformatics 13(6):656-68\n";
	}
	if ($aligner eq "prank") {
		print "Löytynoja A, Goldman N (2008) A model of evolution and structure for multiple sequence alignment. Philosophical Transactions of the Royal Society B 363(1512):3913-9.\n";
	}
	if ($reformalswitch == 1) {
		print "Lyras DP, Metzler D (2014) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach. Submitted\n";
	}
	print "Martin M (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17.\n";
	if ($aligner eq "gramalign") {
		print "Russell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. BMC Bioinformatics 9:306.\n";
	}
	print "Schmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. Bioinformatics 27:863-4.\n";
	if ($aligner eq "clustalo") {
		print "Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539\n";
	}

} elsif ($option eq "-AM3") {

 	# Starting variables
 	my $fastafile = "";
 	my $steps = "";
	my $outname = "";
	my $samplefile = "";

	# Batch mode or interactive mode?
	my $parameter2 = $ARGV[3];
	if (!defined($parameter2) or $parameter2 eq '') {
		print "\nWhat is the name of your fasta file?\n";
  		$fastafile = <STDIN>;
  		chomp $fastafile;
		print "\nWhich maximum number of steps do you want to perform Statistical Parsimony analysis?\nBy default, we used 3 differences as maximum to differenciate populations.\n";
		$steps = <STDIN>;
		chomp $steps;
		print "\nLast required file is SAMPLE. It has a table with two or three columns (it depends of your dataset) and as the following structure:\nSAMPLEID\tDATA1\t(DATA2)\n\nWhat is the SAMPLE file name?\n";
		$samplefile = <STDIN>;
		chomp $samplefile;
		print "\nHow do you name your output files?\n";
		$outname = <STDIN>;
		chomp $outname;
	} else {
		for (my $i = 3; $i < $nrparameters; $i++) {
			if ($ARGV[$i] =~ /^-fasta\=(\S+)/) {
					$fastafile = $1;
			}
			if ($ARGV[$i] =~ /^-nrsteps\=(\d+)/) {
				$steps = $1;
			}
			if ($ARGV[$i] =~ /^-sample\=(\S+)/) {
				$samplefile = $1;
			}
			if ($ARGV[$i] =~ /^-outfile\=(\S+)/) {
				$outname = $1;
			}
		}
	}

	if ($steps eq "") {
		$steps = 3;
	}
	if ($fastafile eq "") {
		HelpingU();
	}
	if ($outname eq "") {
		HelpingU();
	}

	# Analysing your sequences with TCS
	print "\nAnalysing your dataset using TCS (Statistical Parsimony)";
	print "\n--------------------------------------------------------\n";
	Fasta2Phylip($fastafile);
	if ($samplefile ne "") {
		if ($directory eq ".") {
				print "\nMaking TCS analysis in $direcname with $steps steps.\n\n";	
				system("$java -jar -Xmx${freememory}M $TCSfile $direcname $steps");
				TCSloglister($direcname,$samplefile,$fastafile,$outname);
 				unlink glob('*.phyli*');
		} else {
				my $joiningfolder = $direcname . "/" . $directory;
				print "\nMaking TCS analysis in $joiningfolder with $steps steps.\n\n";	
				system("$java -jar -Xmx${freememory}M $TCSfile $joiningfolder $steps");
				TCSloglister($joiningfolder,$samplefile,$fastafile,$outname);
 				unlink glob('*.phyli*');
		}
	} else {
		HelpingU();
	}
		
	# Printing final results
	print "\nIt is done! Your data was classified into different haplotypes according to your parameters. Your fasta file were analysed using Statistical Parsimony in TCS (Clement et al 2000) considering different haplotypes if they have $steps or more differences between them.\nCitation:\nTCS:\t'Clement M, Posada D, Crandall K (2000) TCS: a computer program to estimate gene genealogies. Molecular Ecology 9(10):1657-60.'\n";

} else {
	die "I don't recognize your parameters. Please, read the help ('QRS.pl -h' or 'QRS.pl --help') to know how this program works.";
}
close LOGFILE;

## Different subroutines
# Calculate average length routine
sub AveLength {

	# starting variables
	my $file = shift;
	my $count = 0;
	my $sum = 0;
	my ($min, $max);

	my $seq_io = Bio::SeqIO->new(-file => "$file", -format => 'fasta');
	while (my $seq = $seq_io->next_seq){
		my $len = $seq->length;	
		$count++;
		$sum += $len;
	}
	my $average = int($sum / $count + 0.5);
	return $average;
}

# Cat function
sub Cat {

	my $target = shift;
	my $files = shift;
	my @files = split(",",$files);
	open (SECONDFILE, "> $target") or die "Can't create $target: $!\n";
	foreach my $file (@files) {
		open (FIRSTFILE, $file) or die "Can't open $file: $!\n";
		while (<FIRSTFILE>) {
			print SECONDFILE $_;
		}
	}
	close FIRSTFILE;
	close SECONDFILE;
}

# Counting number of Headers routine
sub CountingHeaders {

	my $fastafile = shift;
	my $nrheaders = 0;
	open (FASTA, $fastafile) or die "There is no way to open $fastafile\n";
	while (my $line = <FASTA>) {
		if ($line =~ /^>/) {
			$nrheaders++;
		}
	}
	close FASTA;
	return $nrheaders;
}

# Declustering Fasta files routine
sub DeclusteringFASTAfile {

	my $ffile = shift;
	my $namefile = shift;
	my $derepfile = shift;

	my $newffile = $ffile;
	$newffile =~ s/\.fasta$/\.declustered\.fasta/;
	$newffile =~ s/\.clust//;
	$newffile =~ s/\.derep//;
	$newffile =~ s/\.sorted//;
	my (@seqNames, %nameSeq);
	my $seqName = "";
	my @splittedheader = ();
	my %oldsequencesheaders = ();
	my %newsequencesheaders = ();
	my %latestsequencesheaders = ();
	my %nrseqheaders = ();
	my %sample = ();
	my %nameofcluster = ();
	my %abundanceofcluster = ();

	open DEREPFILE, $derepfile or die "$derepfile can't be found\n";
	while (my $derepfile = <DEREPFILE>) {
		chomp $derepfile;
		my @derepline = split('\t',$derepfile);
		my $refName = $derepline[0];
		my $Synonims = $derepline[1];
		$oldsequencesheaders{$refName} = $Synonims;
	}
	close DEREPFILE;

	open NAMEFILE, $namefile or die "$namefile can't be found\n";
	while (my $nameline = <NAMEFILE>) {
		chomp $nameline;
		my @nameline = split('\t',$nameline);
		my $refName = $nameline[0];
		$refName =~ s/;size=\d+;//;
		my $anotherelements = $nameline[1];
		$anotherelements =~ s/;size=\d+;//g;
		$newsequencesheaders{$refName} = $anotherelements;
	}
	close NAMEFILE;

	foreach my $key (keys %newsequencesheaders) {
		my @array = split(',',$newsequencesheaders{$key});
		foreach my $element (@array) {
			$element =~ s/$element/$oldsequencesheaders{$element}/;
		}
		my $list = join(',',@array);
		my @newarray = split(',',$list);
		$latestsequencesheaders{$key} = [ @newarray ];
		$nrseqheaders{$key} = scalar(@newarray);
	}

	# Analysing FASTA file
	open FASTA, $ffile or die "$ffile can't be found\n";
	while (my $line = <FASTA>) {
		chomp $line;
		if ($line =~ /^\s*$/) {
			next;
		}
		if ($line =~ /^>(\S+);\S+;/) {
			$seqName = $1;
			$seqName =~ s/(\S+)\_\S+\_\(\S+\)/$1/;
			push(@seqNames, $seqName);
		} else {
			$nameSeq{$seqName} .= $line;
		}
	}
	close FASTA;

	# Analyzing the FASTA file
	open NEWFASTA, ">$newffile" or die "$newffile can't be created.\n";
	foreach my $seqName (@seqNames) {
		if (exists $latestsequencesheaders{$seqName}) {
			for (my $i = 0; $i < $nrseqheaders{$seqName}; $i++) {
				print NEWFASTA ">" . $latestsequencesheaders{$seqName}[$i] . "\n";
				print NEWFASTA $nameSeq{$seqName},"\n";
			}
		}
	}
	close NEWFASTA;

}

sub Fasta2Phylip {

	my $inputfilename = shift;
	my $outfilename = $inputfilename;
	$outfilename =~ s/\.fasta$/\.phylip/;
	
	my $in  = Bio::AlignIO->new(-file => "$inputfilename" , -format => "fasta");
	my $out = Bio::AlignIO->new(-file => "> $outfilename" , -format => "phylip", -interleaved => 0, -wrap_sequential => 0);
	
	while ( my $aln = $in->next_aln() ) {
		$out->write_aln($aln);
	}
}

# GNUv3 "window"
sub GNUv3 {
    die "\n
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.

                            Preamble

  The GNU General Public License is a free, copyleft license for software and other kinds of works.

  The licenses for most software and other practical works are designed to take away your freedom to share and change the works.  By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users.  We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors.  You can apply it to your programs, too.

  When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are designed to make sure that you have the freedom to distribute copies of free software (and charge for them if you wish), that you receive source code or can get it if you want it, that you can change the software or use pieces of it in new free programs, and that you know you can do these things.

  To protect your rights, we need to prevent others from denying you these rights or asking you to surrender the rights.  Therefore, you have certain responsibilities if you distribute copies of the software, or if you modify it: responsibilities to respect the freedom of others.

  For example, if you distribute copies of such a program, whether gratis or for a fee, you must pass on to the recipients the same freedoms that you received.  You must make sure that they, too, receive or can get the source code.  And you must show them these terms so they know their rights.

  Developers that use the GNU GPL protect your rights with two steps: (1) assert copyright on the software, and (2) offer you this License giving you legal permission to copy, distribute and/or modify it.

  For the developers' and authors' protection, the GPL clearly explains that there is no warranty for this free software.  For both users' and authors' sake, the GPL requires that modified versions be marked as changed, so that their problems will not be attributed erroneously to authors of previous versions.

  Some devices are designed to deny users access to install or run modified versions of the software inside them, although the manufacturer can do so.  This is fundamentally incompatible with the aim of protecting users' freedom to change the software.  The systematic pattern of such abuse occurs in the area of products for individuals to use, which is precisely where it is most unacceptable.  Therefore, we have designed this version of the GPL to prohibit the practice for those products.  If such problems arise substantially in other domains, we stand ready to extend this provision to those domains in future versions of the GPL, as needed to protect the freedom of users.

  Finally, every program is threatened constantly by software patents. States should not allow patents to restrict development and use of software on general-purpose computers, but in those that do, we wish to avoid the special danger that patents applied to a free program could make it effectively proprietary.  To prevent this, the GPL assures that patents cannot be used to render the program non-free.

  The precise terms and conditions for copying, distribution and modification follow.

                       TERMS AND CONDITIONS

  0. Definitions.

  'This License' refers to version 3 of the GNU General Public License.

  'Copyright' also means copyright-like laws that apply to other kinds of works, such as semiconductor masks.

  'The Program' refers to any copyrightable work licensed under this License.  Each licensee is addressed as 'you'. 'Licensees' and 'recipients' may be individuals or organizations.

  To 'modify' a work means to copy from or adapt all or part of the work in a fashion requiring copyright permission, other than the making of an exact copy.  The resulting work is called a 'modified version' of the earlier work or a work 'based on' the earlier work.

  A 'covered work' means either the unmodified Program or a work based on the Program.

  To 'propagate' a work means to do anything with it that, without permission, would make you directly or secondarily liable for infringement under applicable copyright law, except executing it on a computer or modifying a private copy.  Propagation includes copying, distribution (with or without modification), making available to the public, and in some countries other activities as well.

  To 'convey' a work means any kind of propagation that enables other parties to make or receive copies.  Mere interaction with a user through a computer network, with no transfer of a copy, is not conveying.

  An interactive user interface displays 'Appropriate Legal Notices' to the extent that it includes a convenient and prominently visible feature that (1) displays an appropriate copyright notice, and (2) tells the user that there is no warranty for the work (except to the extent that warranties are provided), that licensees may convey the work under this License, and how to view a copy of this License.  If the interface presents a list of user commands or options, such as a menu, a prominent item in the list meets this criterion.

  1. Source Code.

  The 'source code' for a work means the preferred form of the work for making modifications to it. 'Object code' means any non-source form of a work.

  A 'Standard Interface' means an interface that either is an official standard defined by a recognized standards body, or, in the case of interfaces specified for a particular programming language, one that is widely used among developers working in that language.

  The 'System Libraries' of an executable work include anything, other than the work as a whole, that (a) is included in the normal form of packaging a Major Component, but which is not part of that Major Component, and (b) serves only to enable use of the work with that Major Component, or to implement a Standard Interface for which an implementation is available to the public in source code form.  A 'Major Component', in this context, means a major essential component (kernel, window system, and so on) of the specific operating system (if any) on which the executable work runs, or a compiler used to produce the work, or an object code interpreter used to run it.

  The 'Corresponding Source' for a work in object code form means all the source code needed to generate, install, and (for an executable work) run the object code and to modify the work, including scripts to control those activities.  However, it does not include the work's System Libraries, or general-purpose tools or generally available free programs which are used unmodified in performing those activities but which are not part of the work.  For example, Corresponding Source includes interface definition files associated with source files for the work, and the source code for shared libraries and dynamically linked subprograms that the work is specifically designed to require, such as by intimate data communication or control flow between those subprograms and other parts of the work.

  The Corresponding Source need not include anything that users can regenerate automatically from other parts of the Corresponding Source.

  The Corresponding Source for a work in source code form is that same work.

  2. Basic Permissions.

  All rights granted under this License are granted for the term of copyright on the Program, and are irrevocable provided the stated conditions are met.  This License explicitly affirms your unlimited permission to run the unmodified Program.  The output from running a covered work is covered by this License only if the output, given its content, constitutes a covered work.  This License acknowledges your rights of fair use or other equivalent, as provided by copyright law.

  You may make, run and propagate covered works that you do not convey, without conditions so long as your license otherwise remains in force.  You may convey covered works to others for the sole purpose of having them make modifications exclusively for you, or provide you with facilities for running those works, provided that you comply with the terms of this License in conveying all material for which you do not control copyright.  Those thus making or running the covered works for you must do so exclusively on your behalf, under your direction and control, on terms that prohibit them from making any copies of your copyrighted material outside their relationship with you.

  Conveying under any other circumstances is permitted solely under the conditions stated below.  Sublicensing is not allowed; section 10 makes it unnecessary.

  3. Protecting Users' Legal Rights From Anti-Circumvention Law.

  No covered work shall be deemed part of an effective technological measure under any applicable law fulfilling obligations under article 11 of the WIPO copyright treaty adopted on 20 December 1996, or similar laws prohibiting or restricting circumvention of such measures.

  When you convey a covered work, you waive any legal power to forbid circumvention of technological measures to the extent such circumvention is effected by exercising rights under this License with respect to the covered work, and you disclaim any intention to limit operation or modification of the work as a means of enforcing, against the work's users, your or third parties' legal rights to forbid circumvention of technological measures.

  4. Conveying Verbatim Copies.

  You may convey verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice; keep intact all notices stating that this License and any non-permissive terms added in accord with section 7 apply to the code; keep intact all notices of the absence of any warranty; and give all recipients a copy of this License along with the Program.

  You may charge any price or no price for each copy that you convey, and you may offer support or warranty protection for a fee.

  5. Conveying Modified Source Versions.

  You may convey a work based on the Program, or the modifications to produce it from the Program, in the form of source code under the terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified it, and giving a relevant date.

    b) The work must carry prominent notices stating that it is released under this License and any conditions added under section 7.  This requirement modifies the requirement in section 4 to 'keep intact all notices'.

    c) You must license the entire work, as a whole, under this License to anyone who comes into possession of a copy.  This License will therefore apply, along with any applicable section 7 additional terms, to the whole of the work, and all its parts, regardless of how they are packaged.  This License gives no permission to license the work in any other way, but it does not invalidate such permission if you have separately received it.

    d) If the work has interactive user interfaces, each must display Appropriate Legal Notices; however, if the Program has interactive interfaces that do not display Appropriate Legal Notices, your work need not make them do so.

  A compilation of a covered work with other separate and independent works, which are not by their nature extensions of the covered work, and which are not combined with it such as to form a larger program, in or on a volume of a storage or distribution medium, is called an 'aggregate' if the compilation and its resulting copyright are not used to limit the access or legal rights of the compilation's users beyond what the individual works permit.  Inclusion of a covered work in an aggregate does not cause this License to apply to the other parts of the aggregate.

  6. Conveying Non-Source Forms.

  You may convey a covered work in object code form under the terms of sections 4 and 5, provided that you also convey the machine-readable Corresponding Source under the terms of this License, in one of these ways:

    a) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by the Corresponding Source fixed on a durable physical medium customarily used for software interchange.

    b) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by a written offer, valid for at least three years and valid for as long as you offer spare parts or customer support for that product model, to give anyone who possesses the object code either (1) a copy of the Corresponding Source for all the software in the product that is covered by this License, on a durable physical medium customarily used for software interchange, for a price no more than your reasonable cost of physically performing this conveying of source, or (2) access to copy the Corresponding Source from a network server at no charge.

    c) Convey individual copies of the object code with a copy of the written offer to provide the Corresponding Source.  This alternative is allowed only occasionally and noncommercially, and only if you received the object code with such an offer, in accord with subsection 6b.

    d) Convey the object code by offering access from a designated place (gratis or for a charge), and offer equivalent access to the Corresponding Source in the same way through the same place at no further charge.  You need not require recipients to copy the Corresponding Source along with the object code.  If the place to copy the object code is a network server, the Corresponding Source may be on a different server (operated by you or a third party) that supports equivalent copying facilities, provided you maintain clear directions next to the object code saying where to find the Corresponding Source.  Regardless of what server hosts the Corresponding Source, you remain obligated to ensure that it is available for as long as needed to satisfy these requirements.

    e) Convey the object code using peer-to-peer transmission, provided you inform other peers where the object code and Corresponding Source of the work are being offered to the general public at no charge under subsection 6d.

  A separable portion of the object code, whose source code is excluded from the Corresponding Source as a System Library, need not be included in conveying the object code work.

  A 'User Product' is either (1) a 'consumer product', which means any tangible personal property which is normally used for personal, family, or household purposes, or (2) anything designed or sold for incorporation into a dwelling.  In determining whether a product is a consumer product, doubtful cases shall be resolved in favor of coverage.  For a particular product received by a particular user, 'normally used' refers to a typical or common use of that class of product, regardless of the status of the particular user or of the way in which the particular user actually uses, or expects or is expected to use, the product.  A product is a consumer product regardless of whether the product has substantial commercial, industrial or non-consumer uses, unless such uses represent the only significant mode of use of the product.

  'Installation Information' for a User Product means any methods, procedures, authorization keys, or other information required to install and execute modified versions of a covered work in that User Product from a modified version of its Corresponding Source.  The information must suffice to ensure that the continued functioning of the modified object code is in no case prevented or interfered with solely because modification has been made.

  If you convey an object code work under this section in, or with, or specifically for use in, a User Product, and the conveying occurs as part of a transaction in which the right of possession and use of the User Product is transferred to the recipient in perpetuity or for a fixed term (regardless of how the transaction is characterized), the Corresponding Source conveyed under this section must be accompanied by the Installation Information.  But this requirement does not apply if neither you nor any third party retains the ability to install modified object code on the User Product (for example, the work has been installed in ROM).

  The requirement to provide Installation Information does not include a requirement to continue to provide support service, warranty, or updates for a work that has been modified or installed by the recipient, or for the User Product in which it has been modified or installed.  Access to a network may be denied when the modification itself materially and adversely affects the operation of the network or violates the rules and protocols for communication across the network.

  Corresponding Source conveyed, and Installation Information provided, in accord with this section must be in a format that is publicly documented (and with an implementation available to the public in source code form), and must require no special password or key for unpacking, reading or copying.

  7. Additional Terms.

  'Additional permissions' are terms that supplement the terms of this License by making exceptions from one or more of its conditions. Additional permissions that are applicable to the entire Program shall be treated as though they were included in this License, to the extent that they are valid under applicable law.  If additional permissions apply only to part of the Program, that part may be used separately under those permissions, but the entire Program remains governed by this License without regard to the additional permissions.

  When you convey a copy of a covered work, you may at your option remove any additional permissions from that copy, or from any part of it.  (Additional permissions may be written to require their own removal in certain cases when you modify the work.)  You may place additional permissions on material, added by you to a covered work, for which you have or can give appropriate copyright permission.

  Notwithstanding any other provision of this License, for material you add to a covered work, you may (if authorized by the copyright holders of that material) supplement the terms of this License with terms:

    a) Disclaiming warranty or limiting liability differently from the terms of sections 15 and 16 of this License; or

    b) Requiring preservation of specified reasonable legal notices or author attributions in that material or in the Appropriate Legal Notices displayed by works containing it; or

    c) Prohibiting misrepresentation of the origin of that material, or requiring that modified versions of such material be marked in reasonable ways as different from the original version; or

    d) Limiting the use for publicity purposes of names of licensors or authors of the material; or

    e) Declining to grant rights under trademark law for use of some trade names, trademarks, or service marks; or

    f) Requiring indemnification of licensors and authors of that material by anyone who conveys the material (or modified versions of it) with contractual assumptions of liability to the recipient, for any liability that these contractual assumptions directly impose on those licensors and authors.

  All other non-permissive additional terms are considered 'further restrictions' within the meaning of section 10.  If the Program as you received it, or any part of it, contains a notice stating that it is governed by this License along with a term that is a further restriction, you may remove that term.  If a license document contains a further restriction but permits relicensing or conveying under this License, you may add to a covered work material governed by the terms of that license document, provided that the further restriction does not survive such relicensing or conveying.

  If you add terms to a covered work in accord with this section, you must place, in the relevant source files, a statement of the additional terms that apply to those files, or a notice indicating where to find the applicable terms.

  Additional terms, permissive or non-permissive, may be stated in the form of a separately written license, or stated as exceptions; the above requirements apply either way.

  8. Termination.

  You may not propagate or modify a covered work except as expressly provided under this License.  Any attempt otherwise to propagate or modify it is void, and will automatically terminate your rights under this License (including any patent licenses granted under the third paragraph of section 11).

  However, if you cease all violation of this License, then your license from a particular copyright holder is reinstated (a) provisionally, unless and until the copyright holder explicitly and finally terminates your license, and (b) permanently, if the copyright holder fails to notify you of the violation by some reasonable means prior to 60 days after the cessation.

  Moreover, your license from a particular copyright holder is reinstated permanently if the copyright holder notifies you of the violation by some reasonable means, this is the first time you have received notice of violation of this License (for any work) from that copyright holder, and you cure the violation prior to 30 days after your receipt of the notice.

  Termination of your rights under this section does not terminate the licenses of parties who have received copies or rights from you under this License.  If your rights have been terminated and not permanently reinstated, you do not qualify to receive new licenses for the same material under section 10.

  9. Acceptance Not Required for Having Copies.

  You are not required to accept this License in order to receive or run a copy of the Program.  Ancillary propagation of a covered work occurring solely as a consequence of using peer-to-peer transmission to receive a copy likewise does not require acceptance.  However, nothing other than this License grants you permission to propagate or modify any covered work.  These actions infringe copyright if you do not accept this License.  Therefore, by modifying or propagating a covered work, you indicate your acceptance of this License to do so.

  10. Automatic Licensing of Downstream Recipients.

  Each time you convey a covered work, the recipient automatically receives a license from the original licensors, to run, modify and propagate that work, subject to this License.  You are not responsible for enforcing compliance by third parties with this License.

  An 'entity transaction' is a transaction transferring control of an organization, or substantially all assets of one, or subdividing an organization, or merging organizations.  If propagation of a covered work results from an entity transaction, each party to that transaction who receives a copy of the work also receives whatever licenses to the work the party's predecessor in interest had or could give under the previous paragraph, plus a right to possession of the Corresponding Source of the work from the predecessor in interest, if the predecessor has it or can get it with reasonable efforts.

  You may not impose any further restrictions on the exercise of the rights granted or affirmed under this License.  For example, you may not impose a license fee, royalty, or other charge for exercise of rights granted under this License, and you may not initiate litigation (including a cross-claim or counterclaim in a lawsuit) alleging that any patent claim is infringed by making, using, selling, offering for sale, or importing the Program or any portion of it.

  11. Patents.

  A 'contributor' is a copyright holder who authorizes use under this License of the Program or a work on which the Program is based.  The work thus licensed is called the contributor's 'contributor version'.

  A contributor's 'essential patent claims' are all patent claims owned or controlled by the contributor, whether already acquired or hereafter acquired, that would be infringed by some manner, permitted by this License, of making, using, or selling its contributor version, but do not include claims that would be infringed only as a consequence of further modification of the contributor version.  For purposes of this definition, 'control' includes the right to grant patent sublicenses in a manner consistent with the requirements of this License.

  Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims, to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contents of its contributor version.

  In the following three paragraphs, a 'patent license' is any express agreement or commitment, however denominated, not to enforce a patent (such as an express permission to practice a patent or covenant not to sue for patent infringement).  To 'grant' such a patent license to a party means to make such an agreement or commitment not to enforce a patent against the party.

  If you convey a covered work, knowingly relying on a patent license, and the Corresponding Source of the work is not available for anyone to copy, free of charge and under the terms of this License, through a publicly available network server or other readily accessible means, then you must either (1) cause the Corresponding Source to be so available, or (2) arrange to deprive yourself of the benefit of the patent license for this particular work, or (3) arrange, in a manner consistent with the requirements of this License, to extend the patent license to downstream recipients. 'Knowingly relying' means you have actual knowledge that, but for the patent license, your conveying the covered work in a country, or your recipient's use of the covered work in a country, would infringe one or more identifiable patents in that country that you have reason to believe are valid.

  If, pursuant to or in connection with a single transaction or arrangement, you convey, or propagate by procuring conveyance of, a covered work, and grant a patent license to some of the parties receiving the covered work authorizing them to use, propagate, modify or convey a specific copy of the covered work, then the patent license you grant is automatically extended to all recipients of the covered work and works based on it.

  A patent license is 'discriminatory' if it does not include within the scope of its coverage, prohibits the exercise of, or is conditioned on the non-exercise of one or more of the rights that are specifically granted under this License.  You may not convey a covered work if you are a party to an arrangement with a third party that is in the business of distributing software, under which you make payment to the third party based on the extent of your activity of conveying the work, and under which the third party grants, to any of the parties who would receive the covered work from you, a discriminatory patent license (a) in connection with copies of the covered work conveyed by you (or copies made from those copies), or (b) primarily for and in connection with specific products or compilations that contain the covered work, unless you entered into that arrangement, or that patent license was granted, prior to 28 March 2007.

  Nothing in this License shall be construed as excluding or limiting any implied license or other defenses to infringement that may otherwise be available to you under applicable patent law.

  12. No Surrender of Others' Freedom.

  If conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, they do not excuse you from the conditions of this License.  If you cannot convey a covered work so as to satisfy simultaneously your obligations under this License and any other pertinent obligations, then as a consequence you may not convey it at all.  For example, if you agree to terms that obligate you to collect a royalty for further conveying from those to whom you convey the Program, the only way you could satisfy both those terms and this License would be to refrain entirely from conveying the Program.

  13. Use with the GNU Affero General Public License.

  Notwithstanding any other provision of this License, you have permission to link or combine any covered work with a work licensed under version 3 of the GNU Affero General Public License into a single combined work, and to convey the resulting work.  The terms of this License will continue to apply to the part which is the covered work, but the special requirements of the GNU Affero General Public License, section 13, concerning interaction through a network will apply to the combination as such.

  14. Revised Versions of this License.

  The Free Software Foundation may publish revised and/or new versions of the GNU General Public License from time to time.  Such new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.

  Each version is given a distinguishing version number.  If the Program specifies that a certain numbered version of the GNU General Public License 'or any later version' applies to it, you have the option of following the terms and conditions either of that numbered version or of any later version published by the Free Software Foundation.  If the Program does not specify a version number of the GNU General Public License, you may choose any version ever published by the Free Software Foundation.

  If the Program specifies that a proxy can decide which future versions of the GNU General Public License can be used, that proxy's public statement of acceptance of a version permanently authorizes you to choose that version for the Program.

  Later license versions may give you additional or different permissions.  However, no additional obligations are imposed on any author or copyright holder as a result of your choosing to follow a later version.

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM 'AS IS' WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

  17. Interpretation of Sections 15 and 16.

  If the disclaimer of warranty and limitation of liability provided above cannot be given local legal effect according to their terms, reviewing courts shall apply local law that most closely approximates an absolute waiver of all civil liability in connection with the Program, unless a warranty or assumption of liability accompanies a copy of the Program in return for a fee.

                     END OF TERMS AND CONDITIONS

            How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest possible use to the public, the best way to achieve this is to make it free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest to attach them to the start of each source file to most effectively state the exclusion of warranty; and each file should have at least the 'copyright' line and a pointer to where the full notice is found.

    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

Also add information on how to contact you by electronic and paper mail.

  If the program does terminal interaction, make it output a short notice like this when it starts in an interactive mode:

    <program>  Copyright (C) <year>  <name of author>
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate parts of the General Public License.  Of course, your program's commands might be different; for a GUI interface, you would use an 'about box'.

  You should also get your employer (if you work as a programmer) or school, if any, to sign a 'copyright disclaimer' for the program, if necessary. For more information on this, and how to apply and follow the GNU GPL, see <http://www.gnu.org/licenses/>.

  The GNU General Public License does not permit incorporating your program into proprietary programs.  If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.  If this is what you want to do, use the GNU Lesser General Public License instead of this License.  But first, please read <http://www.gnu.org/philosophy/why-not-lgpl.html>.
\n";
}

# Help "window"
sub HelpingU {
    system("$perldoc $0");
    exit;
}

# Filtering HMMER results function
sub HMMfiltering_SingleHMM {

	my $table = shift;
	my $fastafile = shift;
	my $cutoff = shift;
	
	# Starting new variables
	my $newfastafile = $fastafile;
	$newfastafile =~ s/\.twice\.fasta/\.hmmfiltered\.fasta/;
	my %data = ();
	my (@seqNames, %nameSeq);
	my $seqName = "";
	my %evalue = ();
	
	# Saving the list
	open TABLE, "$table" or die "There isn't any ". $table. " in your files\n";
	while (<TABLE>) {
		chomp;
		my $target = $_;
		if ($target =~ /^#/) {
			next;
		} else {
			my @line = split('\s+',$target);
			if ($line[8] eq "+") {
				my $reference = $line[0];
				$evalue{$reference} = sprintf("%e", $line[4]);
				$data{$reference} = 1;
			}
		}
	}       
	close TABLE;
	
	# Saving the original fasta file in memory
	open FASTA, $fastafile or die "couldn't open $fastafile: $!\n";
	while (my $line = <FASTA>) {
		chomp $line;
		if ($line =~ /^\s*$/) {
			next;
		}
		if ($line =~ /^>(\S+)/) {
			$seqName = $1;
			if (exists $data{$seqName}) {
				if ($evalue{$seqName} < $cutoff) {
					push(@seqNames, $seqName);
				}
			}
		} else {
			$nameSeq{$seqName} .= $line;
		}
	}
	close FASTA;

	open NEWFASTA, ">$newfastafile" or die "Can't create $newfastafile\n";
	foreach $seqName (@seqNames) {
		print NEWFASTA ">",$seqName, "\n";
		print NEWFASTA $nameSeq{$seqName},"\n";
	}
	close NEWFASTA;
}

# Modifying headers in original fasta file
sub ModifyingHeaders {
	my $inFastaFile = shift;
	my (@seqNames, %nameSeq);
	my $seqName = "";

	# Putting the modified fasta file into memory
	open (ORIGFASTA, $inFastaFile) or die "Couldn't open $inFastaFile\n";
	my $line = "";
	while ($line = <ORIGFASTA>) {
		chomp $line;
		if ($line =~ /^\s*$/) {
			next;
		}
		if ($line =~ /^>(.*)/) {
			$seqName = $1;
			$seqName =~ s/\./\_/g;
# 			$seqName =~ s/\s+/\_/g;
			push(@seqNames, $seqName);
		} else {
			$nameSeq{$seqName} .= $line;
		}
	}
	close ORIGFASTA;

	# Making the modified fasta file
	open MODFASTA, ">$inFastaFile" or die "couldn't open $inFastaFile: $!\n";
	foreach $seqName (@seqNames) {
		print MODFASTA ">",$seqName,"\n";
		print MODFASTA $nameSeq{$seqName},"\n";
	}
	close MODFASTA;
}

# Reverse Complementary sequences function
sub RevComSeq {
	my $inFastaFile = shift;
	my $outFastaFile = shift;
	my (@seqNames, %nameSeq);
	my $seqName = "";

	# Putting the modified fasta file into memory
	open (ORIGFASTA, $inFastaFile) or die "Couldn't open $inFastaFile\n";
	my $line = "";
	while ($line = <ORIGFASTA>) {
		chomp $line;
		if ($line =~ /^\s*$/) {
			next;
		}
		if ($line =~ /^>(.*)/) {
			$seqName = "RC_" . $1;
			$seqName =~ s/\./\_/g;
			push(@seqNames, $seqName);
		} else {
			$nameSeq{$seqName} .= $line;
		}
	}
	close ORIGFASTA;

	# Making the Reverse Complementary fasta file
	open RCFASTA, ">$outFastaFile" or die "couldn't open $outFastaFile: $!\n";
	foreach $seqName (@seqNames) {
		my $seq = reverse $nameSeq{$seqName};
		$seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy[]/TVGHCDKNYSAABWXRtvghcdknysaabwxr][/;
		print RCFASTA ">",$seqName, "\n";
		print RCFASTA $seq,"\n";
	}
	close RCFASTA;
}

# Splitting samples and removing barcodes function
sub Splitting {

	# Reading all input parameters
	my $fastafile = shift;
	my $oligosfile = shift;
	my $designfile = shift;
	my $barcodereverse = shift;
	my $mismatches = shift;
	my $bcdrvs = "";
	
	if ($barcodereverse eq '-brn') {
		$bcdrvs = 0;
	} else {
		$bcdrvs = 1;
	}

	# Starting internal parameters
	my %barcode = ();
	my %sample = ();
	my @listmids = ();
	my %constructedfwdprimers = ();
	my %constructedrvsprimers = ();
	my %rcprimer = ();
	my $fwd = "";
	my $rvs = "";
	my $barname = "";
	my ($seqadapfor, $seqadaprev);
	
	# Saving the OLIGOS file in memory
	open OLIGOS, "$oligosfile" or die "There isn't ".$oligosfile; # Reading the oligos file
	while (my $oligosline=<OLIGOS>) {
		if ($oligosline =~ /^#/) {
			next;
		} elsif ($oligosline =~ /^seqadapfor\t(\S+)/) {
			$seqadapfor = $1;
		} elsif ($oligosline =~ /^seqadaprev\t(\S+)/) {
			$seqadaprev = $1;
		} elsif ($oligosline =~ /^forward\t(\S+)/) {
			$fwd = $1;
		} elsif ($oligosline =~ /^reverse\t(\S+)/) {
			$rvs = $1;
		} elsif ($oligosline =~ /^barcode\t(\S+)\t(\S+)/) {
			$barname = $2;
			$barcode{$barname} = $1;
		} else {
			next;
		}
	}
	close OLIGOS;
	 
	# Saving the DESIGN file in memory
	my @designline = ();
	open DESIGN, "$designfile" or die "There isn't ".$designfile; # Reading the design file
	while (my $designline=<DESIGN>) {
		if ($designline =~ /^#/) {
			next;
		} else {
			chomp($designline);		
			@designline = split(/\t/, $designline);
			if ($bcdrvs == 0) {
				my $midfwd = $designline[0];
				$sample{$midfwd} = $designline[1];
			} else {
				my $midfwd = $designline[0];
				my $midrvs = $designline[1];
				$sample{$midfwd}{$midrvs} = $designline[2];
			}
		}
	}
	close DESIGN;
	 
	# Reconstruction of primers
	foreach my $mid (keys %barcode) {
		push(@listmids, $mid); # Generating a list of mids
		if (defined $seqadapfor) {
			$constructedfwdprimers{$mid} = $seqadapfor.$barcode{$mid}.$fwd;
		} else {
			$constructedfwdprimers{$mid} = $barcode{$mid}.$fwd;
		}
		if (defined $seqadaprev) {
			if ($bcdrvs == 1) {
				$constructedrvsprimers{$mid} = $seqadaprev.$barcode{$mid}.$rvs;
			} else {
				$constructedrvsprimers{$mid} = $seqadaprev.$rvs;
			}
		} else {
			if ($bcdrvs == 1) {
				$constructedrvsprimers{$mid} = $barcode{$mid}.$rvs;
			} else {
				$constructedrvsprimers{$mid} = $rvs;
			}
		}
		$rcprimer{$mid} = reverse_complement($constructedrvsprimers{$mid});
	}
	 
	# Processing FASTA file
	my $record_seq = "";
	my $number;
	foreach my $mid (@listmids) {
		if ($bcdrvs == 1) {
			foreach my $dim (@listmids) {
				if (exists $sample{$mid}{$dim}) {
					my $newfastafile = "$sample{$mid}{$dim}.fasta";
					open FASTA, "$fastafile" or die "There isn't ".$fastafile; # Reading the fasta file
					open TEMPFILE, ">$newfastafile" or die "Can't create ".$newfastafile;
					$number = 0;
					while (my $line=<FASTA>) {
						if ($line =~ /^>/) {
							next;
						} else {
							$record_seq = $record_seq . $line;
							if ($record_seq =~ /^$constructedfwdprimers{$mid}/ and $record_seq =~ /$rcprimer{$dim}$/ and exists $sample{$mid}{$dim}) {
								$number++;
								print TEMPFILE ">" . $sample{$mid}{$dim} . "_" . $number . "\n";
								print TEMPFILE $record_seq;
								$record_seq = "";
							} elsif (defined $mismatches) {
								if (amatch($constructedfwdprimers{$mid},[$mismatches],$record_seq) and amatch($rcprimer{$dim},[$mismatches],$record_seq) and exists $sample{$mid}{$dim}) {
									$number++;
									print TEMPFILE ">" . $sample{$mid}{$dim} . "_" . $number . "\n";
									print TEMPFILE $record_seq;
									$record_seq = "";
								} else {
									$record_seq = "";
								}
							} else {
								$record_seq = "";
							}
						}
					}
					close TEMPFILE;
				} else {
				}
			}
		} else {
			my $newfastafile = "$sample{$mid}.fasta";
			open FASTA, "$fastafile" or die "There isn't ".$fastafile; # Reading the fasta file
			open TEMPFILE, ">$newfastafile" or die "Can't create ".$newfastafile;
			$number = 0;
			while (my $line=<FASTA>) {
				if ($line =~ /^>/) {
					next;
				} else {
					$record_seq = $record_seq . $line;
					if ($record_seq =~ /^$constructedfwdprimers{$mid}/ and $record_seq =~ /$rcprimer{$mid}$/ and exists $sample{$mid}) {
						$number++;
						print TEMPFILE ">" . $sample{$mid} . "_" . $number . "\n";
						print TEMPFILE $record_seq;
						$record_seq = "";
					} elsif (defined $mismatches) {
						if (amatch($constructedfwdprimers{$mid},[$mismatches],$record_seq) and amatch($rcprimer{$mid},[$mismatches],$record_seq) and exists $sample{$mid}) {
							$number++;
							print TEMPFILE ">" . $sample{$mid} . "_" . $number . "\n";
							print TEMPFILE $record_seq;
							$record_seq = "";
						} else {
							$record_seq = "";
						}
					} else {
						$record_seq = "";
					}
				}
			}
			close TEMPFILE;
		}
	}
#	$number = 0;
	close FASTA;
	
	# Removing complete primer of each sequence
	foreach my $mid (keys %constructedfwdprimers) {
		if ($bcdrvs == 1) {
			foreach my $dim (keys %constructedrvsprimers) {
				if (exists $sample{$mid}{$dim}) {
					my $tempfastafile = "$sample{$mid}{$dim}.fasta";
					my $partfastafile = "$sample{$mid}{$dim}.part.fasta";
					my $newfastafile = "$sample{$mid}{$dim}.trimmed.fasta";
					system("$cutadaptprogram -g \^$constructedfwdprimers{$mid} $tempfastafile > $partfastafile 2> /dev/null");
					system("$cutadaptprogram -b $rcprimer{$dim} $partfastafile > $newfastafile 2> /dev/null");
					unlink $tempfastafile;
					unlink $partfastafile;
				} else {
				}
			}
		} else {
			my $tempfastafile = "$sample{$mid}.fasta";
			my $partfastafile = "$sample{$mid}.part.fasta";
			my $newfastafile = "$sample{$mid}.trimmed.fasta";
			system("$cutadaptprogram -g \^$constructedfwdprimers{$mid} $tempfastafile > $partfastafile 2> /dev/null");
			system("$cutadaptprogram -b $rcprimer{$mid} $partfastafile > $newfastafile 2> /dev/null");
			unlink $tempfastafile;
			unlink $partfastafile;
		}
	}
	
	sub reverse_complement {
		my $dna = shift;
		my $revcomp = reverse($dna);
		$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
		return $revcomp;
	}
}

# TCS log interpreter routine
sub TCSloglister {

	#Starting variables
	my $dirname = shift;
	my $samplefile = shift;
	my $fastafile = shift;
	my $nameoutfile = shift;
	my $abstempfile = $dirname . "/" . $nameoutfile . ".abs.list";
	my $absresultsfile = $dirname . "/" . $nameoutfile . ".abs.csv";
	my $haplotypesfastafile = $dirname . "/" . $nameoutfile . ".fasta";
	my @logfiles = ();
	my %samples = ();
	my ($seqName, @seqNames, %nameSeq);

	open SAMPLEFILE, $samplefile or die "$samplefile can't be found\n";
	while (my $line=<SAMPLEFILE>) {
		if ($line =~ /^#/) {
			next;
		} else {
			my @line = split('\t',$line);
			my $sample = $line[0];
			$samples{$sample} = 1;
		}
	}
	close SAMPLEFILE;

	@logfiles = glob($dirname . '/*phylip.log');
	foreach my $tcslogfile (@logfiles) {

		# Starting the variables
		my %acceptedids = ();
		my %haplotypes = ();
		my %haplolength = ();
		my %ids = ();
		my $network;
		my $id;
		my %nrseqspernet = ();
		my %totalnrseqs = ();
		my %abshash = ();
		my %net_list = ();
		my %representative = ();
		
		# Analyzing the TCS log file
		open TCSFILE, $tcslogfile or die "$tcslogfile can't be found\n";
		while (my $line=<TCSFILE>) {
			if ($line =~ /^\s-\s(\S+)\s+\:/) {
				my $seqid = $1;
				$acceptedids{$seqid} = 1;
				$line =~ s/-//;
				$line =~ s/\://;
				$line =~ s/\s+/\t/g;
				$line =~ s/^\t//;
				$line =~ s/\t$//;
				my @line = split('\t',$line);
				$haplotypes{$seqid} = [ @line ];
				$haplolength{$seqid} = scalar (@line);
			} elsif ($line =~ /^\*\*\*\sNetwork\s(\d+)/) {
				$network = $1;	
			} elsif ($line =~ /^(\S+)\s+weigth\s\=\s\d\.\d*/) {
				$id = $1;
				if (exists $acceptedids{$id}) {
					for my $sample (sort keys %samples) {
						for (my $j = 0; $j < $haplolength{$id}; $j++) {
							if ($haplotypes{$id}[$j] =~ /^$sample/) {
								if (exists $ids{$sample}{$network}) {
									$ids{$sample}{$network} .= ",$haplotypes{$id}[$j]";
								} else {
									$ids{$sample}{$network} = $haplotypes{$id}[$j];
								}
							}
						}
					}
				}
			} elsif ($line =~ /^Biggest\soutgroup\sprobability\sis\s(\S+)/) {
				$id = $1;
				$representative{$id} = $network;
			}
		}
		close TCSFILE;
		
		open ABSTEMPFILE, ">> $abstempfile" or die "$abstempfile can't be created\n";
		for my $sample (sort { lc($a) cmp lc($b) } keys %ids) {
			$totalnrseqs{$sample} = 0;
			for my $net (sort { $a <=> $b } keys $ids{$sample}) {
				my @ids = split(',', $ids{$sample}{$net});
				$nrseqspernet{$sample}{$net} = scalar @ids;
				print ABSTEMPFILE "$sample\t$net\t$nrseqspernet{$sample}{$net}\n";
			}
		}
		close ABSTEMPFILE;

		# Converting the absolute tempfile into a absolute matrixfile
		open ABSTEMPFILE, "$abstempfile" or die "Can't exist $abstempfile\n";
		while (my $listline = <ABSTEMPFILE>) {
			chomp $listline;
			my @line = split ("\t", $listline);
			$abshash{$line[0]}{$line[1]} = $line[2];
			$net_list{$line[1]} = 1;
		}
		close ABSTEMPFILE;
		unlink $abstempfile;

		# Printing the absolute matrix
		open ABSRESULTSFILE, "> $absresultsfile" or die "Can't create $absresultsfile\n"; # Reading the list
		print ABSRESULTSFILE "sample";
		my @net_list = sort { $a <=> $b } keys (%net_list);
		foreach my $net (@net_list) {
			print ABSRESULTSFILE "\t$net";
		}
		print ABSRESULTSFILE "\n";

		foreach my $sample (sort { lc($a) cmp lc($b) } keys %abshash) {
			print ABSRESULTSFILE "$sample";
			foreach my $net (@net_list){
				if (exists $abshash{$sample}{$net}) {
					print ABSRESULTSFILE "\t$abshash{$sample}{$net}";
				} else {
					print ABSRESULTSFILE "\t0";
				}
			}
			print ABSRESULTSFILE "\n";
		}
		close ABSRESULTSFILE;

 		# Printing a new FASTA file with only TCS-types sequences
 		open FASTAFILE, "$fastafile" or die "$fastafile can't be opened\n";
 		while (my $line = <FASTAFILE>) {
 			chomp $line;
 			if ($line =~ /^\s*$/) {
 				next;
 			}
 			if ($line =~ /^>(\S+)/) {
 				$seqName = $1;
 				$seqName =~ s/{^.{10}}.*$/$1/;
 				if (exists $representative{$seqName}) {
 					push(@seqNames, $seqName);
 				} else {
 				}
 			} else {
 				$nameSeq{$seqName} .= $line;
 			}
 		}
 		close FASTAFILE;
 
 		open NEWFASTA, ">$haplotypesfastafile" or die "$haplotypesfastafile can't be created\n";
		foreach (sort {$representative{$a} <=> $representative{$b}} keys %representative) {
			my $element = $_;
			foreach my $seqName (@seqNames) {
				if ($seqName eq $element) {
					print NEWFASTA ">Network_" . $representative{$element} . "\n";
					print NEWFASTA $nameSeq{$element},"\n";
				}
			}
		}
 		close NEWFASTA;
	}
}

# Analyzing USEARCH tables
sub USEARCHanalysis {
	my $usearchtable = shift;
	my $newlistfile = $usearchtable;
	$newlistfile =~ s/\.uctbl$/\.list/;

	my %abcltr = ();
	my %nmcltr = ();
	my %listcltr = ();

	# Analyzing the UCLUST clusters file
	open TABLE, "$usearchtable" or die "Can't exist ".$usearchtable;
	while (my $line=<TABLE>) {
		if ($line =~ /^#/) {
			next;
		} elsif ($line =~ /^C/) {
			my @line = split('\s+',$line);
			my $repseq = $line[8];
			$nmcltr{$repseq} = $line[1];
			$abcltr{$repseq} = $line[2];
		} else {		
		}
	}
	close TABLE;
	
	# Re-analyzing the UCLUST clusters file
	open TABLETWICE, "$usearchtable" or die "Can't exist ".$usearchtable;
	while (my $line=<TABLETWICE>) {
		if ($line =~ /^#/) {
			next;
		} elsif ($line =~ /^S\s+/) {
			my @line = split('\s+',$line);
			if (exists $nmcltr{$line[8]}) {
				$listcltr{$line[8]} = $line[8];
			}
		} elsif ($line =~ /^H/) {
			my @line = split('\s+',$line);
			my $newseq = $line[8];
			if (exists $nmcltr{$line[9]}) {
				if ($listcltr{$line[9]} ne "") {
					$listcltr{$line[9]} .= ",$newseq";
				} else {
					$listcltr{$line[9]} = $line[9];
				}
			}
		} else {		
		}
	}
	close TABLETWICE;

	open NEWLIST, ">$newlistfile" or die "Can't exist ".$newlistfile;
	for my $key (keys %listcltr) {
		print NEWLIST "$key\t$listcltr{$key}\n";
	}
	close NEWLIST;

}

## Help in POD format

=head1 NAME

Quantification of Representative Sequence (QRS.pl)

=head1 DESCRIPTION

QRS is a script written in Perl 5.14.2 that allows analysing 454 pyrosequencing platform sequences automatically to study population structure. This program may operate either in batch processing mode where all sequences are automatically analysed in an unsupervised way or it may interact with the user at various checkpoints if no parameters have been specified prior to execution.

=head1 USAGE

The program QRS can be executed as batch command or as an interactive program if no parameters are defined. To understand how QRS works in batch mode, we provide several use-case scenarios.

=head2 Creating a reference HMM file (-RM)

=head3 Basic parameters

When this option is used, the batch command can be executed using only basic parameters:
QRS.pl -RM -folder=test -reffiles=test.fasta

In this case, we call QRS to create a reference HMM file of test.fasta in the folder test using default parameters (the aligner is PRANK and QRS will call ReformAl to fine-tune the alignment).
However, QRS can create more than a single HMM file if you specify different reference files joined with a single dash (-) in the -reffiles parameter, like in the following example, where the program will create two HMM files using default parameters:
QRS.pl -RM -folder=test -reffiles=test1.fasta-test2.fasta

=head3 Alignment options

Another option that you can modify is the multiple sequence alignment program. You can choose the program with the -aligner parameter. PRANK is the default aligner program, like in this example where QRS will create a single reference HMM file:
QRS.pl -RM -folder=test -aligner=prank -reffiles=test.fasta

To employ a different aligner for the sequence alignment task simply modify the -aligner parameter. In the following example, QRS will call MUSCLE to align a reference file to create a single reference HMM file:
QRS.pl -RM -folder=test -aligner=muscle -reffiles=test.fasta

Here, we use MAFFT with an iterative global alignment algorithm called G-INS-i to perform an accurate alignment:
QRS.pl -RM -folder=test -aligner=mafft-ginsi -reffiles=test.fasta

The last option you may modify is whether to use ReformAl to fine-tune the alignment or not. By default, this option is enabled (-reformal=yes), like in the following example, where QRS will create two reference HMM files after calling GramAlign to align your reference sequences:
QRS.pl -RM -folder=test -aligner=gramalign -reformal=yes -reffiles=test1.fasta-test2.fasta

If you want to disable this option, you have to type -reformal=no. In this example, QRS will align your three reference files using the default aligner (i.e., PRANK) and will not curate the alignments before creating three reference HMM files:
QRS.pl -RM -folder=test -reformal=no -reffiles=test1.fasta-test2.fasta-test3.fasta

In the following example, we call QRS to align four reference files using MAFFT with their iterative local alignment algorithm (L-INS-i) and we do not want to fine-tune the alignments:
QRS.pl -RM -folder=test -aligner=mafft-linsi -reformal=no -reffiles=test1.fasta-test2.fasta-test3.fasta-test4.fasta

=head2 Describing the 454 data set (-AM1)

=head3 Basic parameters

When this option is used, the batch command cannot be executed using only basic parameters. Instead, you have to specify if you have a 454 FASTA file or a 454 FASTQ file. In this example the input file is a FASTA file: 
QRS.pl -AM1 -folder=test -informat=fasta -fasta=file.fna

Here, the input file is a FASTQ file:
QRS.pl -AM1 -folder=test -informat=fastq -fastq=file.fastq

=head3 Another input files
If you have a FASTA file as input, you can add a quality file to do all basic statistics on the base quality using the parameter -quality:
QRS.pl -AM1 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual

If you have a paired FASTQ file, you have to add the parameter -paired=yes and the name of the paired FASTQ:
QRS.pl -AM1 -folder=test -informat=fastq -fastq=file.fastq -paired=yes -fastq2=file2.fastq

=head2 Processing a 454 data set (-AM2)

=head3 Basic parameters

Before executing this part, you have to evaluate the characteristics of your data set according to the previous step (see AM1 section). As all 454 pyrosequencing data sets are different depending on the case study, there are no default parameters to filter and trim sequences according to length, quality and complexity. 
However, in a hypothetical case that you do not need to filter your data set, the input parameters are the input files added as describe before (see AM1 section). Moreover, you have to add the reference HMM file (created in -RM step, see RM section for more details), oligos and design files (read below to know more about these files) and specify if you have reverse barcoded primers or not. In the following example, we use a 454 FASTA file with quality file and our data set has reverse barcoded primers (-bry): 
QRS.pl -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -bry

The oligos file is a plain text file that has two columns splitted by tabs (except if the label is barcode, where there are a third column). The first column contains always a label that indicates forward or reverse adapter, forward or reverse primer and barcode. The second contains the DNA sequence for each element and, in the case of barcode label, the third one indicates the barcode ID (MID1, MID2, MID3...). An example of this file is as follows:
		   seqadapfor	SEQUENCE
		   seqadaprev	SEQUENCE
		   forward	SEQUENCE
		   reverse	SEQUENCE
		   barcode	SEQUENCE	BARID

The design file is also a plain text file that has two or three columns (depending of the use of reverse barcoded primers) splitted by tabs. An example of these file is as follows:
		   BARID1	(BARID2)	SAMPLEID

Here, BARID1 is the Barcode ID for the forward primer, BARID2 is the Barcode ID for the reverse primer (if exists) and SAMPLEID is the name of the sample.

=head3 Filtering by length

As it is said in AM2.1 section, there are no default parameters to filter and trim sequences according to length, quality and complexity because it depends of your data set. In the following example calls QRS to accept sequences that are greater than 300 bp in a 454 FASTQ file and have no reverse barcoded primers:
QRS.pl -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -minlen=300 -oligos=oligos.csv -design=design.csv -brn

In this example, QRS is executed to accept sequences that have between 355 and 500 bp (the data is given as 454 FASTA file without quality file):
QRS.pl -AM2 -folder=test -informat=fasta -fasta=file.fna -hmmfile=test.hmm -minlen=355 -maxlen=500 -oligos=oligos.csv -design=design.csv -bry

=head3 Filtering by GC content

Another possible filter is based on the GC content. In the following example QRS is called to accept sequences that have more than 60% GC content in a 454 FASTQ file and have no reverse barcoded primers:
QRS.pl -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -mingc=60 -oligos=oligos.csv -design=design.csv -brn

In this example, QRS is executed to accept sequences that have between 40 and 50% GC content (the data is given as 454 FASTA file without quality file):
QRS.pl -AM2 -folder=test -informat=fasta -fasta=file.fna -hmmfile=test.hmm -mingc=40 -maxgc=50 -oligos=oligos.csv -design=design.csv -bry

=head3 Filtering and trimming according to quality

We offer another filter based on base quality. In this example, QRS accepts only sequences that have at least a mean sequence quality value of 25:
QRS.pl -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -minqual=25 -oligos=oligos.csv -design=design.csv -bry

If you have a FASTA file, it is convenient to have a quality file to filter by base quality. In the following example, QRS will filter an input 454 FASTA file by length (accepting only sequences that are between 375 and 480 bp long) and quality (rejecting sequences that have a mean sequence quality score smaller than 28). This data set has no reverse barcoded primers:
QRS.pl -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -minlen=375 -maxlen=480 -minqual=28 -brn

The following option is used to trim bases that have an insufficient quality score. In the following example, QRS will trim all nucleotides with a quality value less than 25 in the beginning and at the end of the sequences in a non-reverse barcoded 454 data set:
QRS.pl -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -trimqual=25 -brn

=head3 Trimming poly-A/Ts tails

It is feasible to trim poly-A/Ts tails. In this example, QRS trims the regions that have more than three followed A/T at the beginning or the end of the sequences:
QRS.pl -AM2 -folder=test -informat=fastq -FASTQ=file.fastq -hmmfile=test.hmm -trimtails=3 -oligos=oligos.csv -design=design.csv -bry

=head3 Filtering by ambiguous characters

The number of allowed ambiguous characters in your data set can be modified with the parameter -allowns. In the following example, QRS accepts all sequences that have at most one ambiguous character in the sequence:
QRS.pl -AM2 -folder=folder -informat=fastq -fastq=file.fastq -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowns=1 -bry

However, it is not recommended to allow ambiguous characters as they might indicate bad quality nucleotides (Huse et al. 2007). In the following example, QRS removes all sequences that have more than one ambiguous character:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowns=0 -bry

=head3 Filtering by complexity

Finally, if you want to remove all low complexity sequences like homopolymers, you can do this with-filmet and -filthr parameters. The first argument defines the method to remove this kind of sequences (DUST (Morgulis et al. 2006) or Entropy-based filter) and the second argument defines the threshold for the filtering method. For more details on these filtering methods, see PrinSeq manual (http://prinseq.sourceforge.net/manual.html#QCCOMPLEXITY). In the following example, we use DUST to remove all sequences with a complexity greater than 7:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -filmet=dust -filthr=7 -bry

In this example, we use entropy to remove all sequences with a complexity smaller than 70:
QRS.pl -AM2 -folder=folder -informat=FASTA -FASTA=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -filmet=entropy -filthr=70 -bry

In this example, a 454 data set (that consists in a FASTQ file) is filtered by length (accepting only sequences that are between 390 and 500 bp long), base quality (removing all sequences with mean sequence quality score less than 28), and low-complexity based on the DUST filter (considering all sequences with values greater than 5 as low complexity sequences):
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minlen=390 -maxlen=500 -minqual=28 -filmet=dust -filthr=5 -bry

=head3 Classifying sequences as a specific marker
You can also modify the maximum allowed e-value to classify a sequence as a specific marker using a HMM profile (-hmmthr). By default, this argument is set to 10-10 as then similar results as with the usage of BLAST are achieved (Altschul et al. 1990). We do not recommend changing this value. However, you can change the value like in the following example where QRS classifies a sequence as a specific marker if the e-value is smaller than 10-25:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -hmmthr=1E-25 -bry

=head3 Assigning sequences to samples

If you want to modify the maximum number of allowed mismatches to detect the barcodes in your data set, you can use the parameter -allowmis. In this example, QRS classifies a non reverse barcoded data set in different samples with an exact match of the barcodes:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowmis=0 -brn

You can define a percentage to define the maximum value for mismatches. By default, we consider a maximum percentage of mismatches of 1% to detect the barcodes. Although we do not recommend modifying this value, it is possible to change it like in this example, where QRS considers a maximum percentage of 2% to detect the barcodes in a reverse barcoded data set:
QRS.pl -AM2 -folder=folder -informat=FASTQ -FASTQ=file.FASTQ -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowmis=2% -bry

=head3 Clustering step

Another parameter you can modify is -cutoff. This parameter is used to cluster sequences in order to de-noise your data set, i. e., to remove all spurious nucleotides across all sequences. By default, QRS calculates this value according to the CD-HIT-OTU algorithm (Li et al. 2012) but you can define the value by a number between 0.00 and 1.00. For example, if you want to cluster all sequences that have 99.5% of similarity, you can type:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -cutoff=0.995 -brn

The parameter -minclustersize specifies the minimum size of cluster to consider that analysed sequences are not sequencing artifacts. By default, this filter is enabled and all clusters that have at least three sequences are considered as good clusters. If you want to disable this argument, you have to type -minclustersize=0, like in this example:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minclustersize=0 -bry

=head3 Alignment step

The last options you may modify concern the use of the alignment program and ReformAl to fine-tune the alignment. For more information, see the RM examples.

=head3 The big example

In the following example, QRS retrieves all sequences that have more than 350 bp, a mean sequence quality score of 30, no ambiguous characters, and a sequence complexity greater than 79 according to the entropy-based filter. This data set has reverse paired barcodes and the script discards all clusters that have only one sequence. Finally, QRS uses KAlign to align all accepted sequences and this alignment is not fine tuned:
QRS.pl -AM2 -folder=folder -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minlen=350 -minqual=30 -allowns=0 -filmet=entropy -filthr=79 -bry -minclusterzise=2 -aligner=kalign -reformal=no

=head2 Classifying all 454 sequences into TCS-types (-AM3)

=head3 Basic parameters

When this option is used, the batch command can be executed with only basic parameters, like in the following example:
QRS.pl -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -outfile=allsamples

In the previous example, QRS will use myaligneddata alignment to run TCS for the data and then retrieve a TCS-types FASTA file called allsamples.fasta and a frequencies matrix file called allsamples.abs.csv. The program makes use of some information from the samples.csv file. This file which is a plain text file that contains the following information always splitted by tabs:

SAMPLEID	DATA1	(...)	DATAN

Here, SAMPLEID is the name of the sample and DATA are data you want to put here (Place, year, species...). 

=head3 Statistical parsimony parameters

You can modify the maximum number of steps you want to use to determine TCS-types with the -nrsteps parameter. By default, QRS considers that two sequences belong to the same TCS-type if it has three differences or less (-nrsteps=3), like in the following example:
QRS.pl -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -nrsteps=3 -outfile=allsamples

In this example, the maximum number of steps to determine if two sequences belong from the same TCS-type is five differences:
QRS.pl -AM3 -folder=test -fasta=myaligneddata.fasta -sample=sample.csv -nrsteps=5 -outfile=allsamples

=head1 REQUIREMENTS

Before using this pipeline, the following Perl modules and programs should be installed:

* Perl modules:
- BioPerl (Bio::AlignIO and Bio::SeqIO) (Stajich et al. 2002)
- File::Find::Rule
- File::Which
- JSON
- String::Approx
- Sys::CPU
- Sys::MemInfo

* Programs:
- PrinSeq (Schmieder & Edwards 2011): it is used to generate summary statistics of sequence and quality data and to filter, reformat and trim 454 data. This program is freely available at http://prinseq.sourceforge.net under the GPLv3 licence.
- CutAdapt (Martin 2011): a Python script that removes primers and adapters. This program is publicly available at http://code.google.com/p/cutadapt under the MIT licence.
HMMER 3.1b (Finn et al 2011): it is used to search sequences using probabilistic models called profile hidden Markov models (profile HMM). This program is freely available at http://hmmer.janelia.org/ under GPLv3 licence.
- USEARCH (Edgar 2010): it is used to cluster sequences and detect chimaeras according to the UCHIME algorithm (Edgar et al. 2011). This program is available after requiring it according to the authors' page (http://www.drive5.com/usearch/).
- At least one of the following aligners:

* Clustal Omega (Sievers et al. 2011): a hybrid aligner that mixes progressive multiple sequence alignments with the use of Markov models. This program is freely available at http://www.clustal.org/omega/ under GPLv2 licence.
* FSA (Bradley et al. 2009): a probabilistic aligner that uses the sequence annealing technique (Schwartz & Pachter 2007) for constructing a multiple alignment from pairwise homology estimation. FSA is publicly available at http://fsa.sourceforge.net/ under GPL licence.
* GramAlign (Russell et al. 2008): a time-efficient progressive aligner which estimates distances according to the natural grammar present in nucleotide sequences. This program is publicly available at http://bioinfo.unl.edu/gramalign.php
* KAlign (Lassmann & Sonnhammer 2005): a progressive aligner based on Wu-Manber string-matching algorithm. KAlign is publicly available at http://msa.sbc.su.se under GPLv2 licence.
* MAFFT (Katoh & Standley 2013): an iterative/progressive alignment program that is publicly available at http://mafft.cbrc.jp/alignment/software/ under BSD licence.
* MUSCLE (Edgar 2004): an iterative aligner that it is available at http://www.drive5.com/muscle/.
* PRANK (Loytynoja & Goldman 2005): a probabilistic multiple alignment program based on maximum likelihood methods used in phylogenetics that considers evolutionary distances between sequences. This program is publicly available at http://code.google.com/p/prank-msa under GPLv3 licence.

- TCS (Clement et al. 2000): a Java computer program to estimate phylogenetical networks and intraspecific genealogies according to statistical parsimony (Templeton et al. 1992). This program is freely available at http://darwin.uvigo.es/software/tcs.html under GPLv2 licence.

Although the aforementioned aligners have already been tested with the QRS, other aligners may be added to use in this pipeline. If you want to work with another alignment program, please feel free to contact with the author (tortuero@bio.lmu.de) to include it in the source code of QRS. 
Finally, an optional recommended requirement might be the installation of ReformAlign (Lyras and Metzler 2014). ReformAlign is a recently proposed profile-based meta-alignment approach that aims to fine-tune existing alignments via the employment of standard profiles.

=head1 HISTORY

v 1.0.0 - First version of the program. It should work fine.

=head1 REFERENCES

* Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) Basic local alignment search tool. I<Journal of Molecular Biology> B<215:>403-10.

* Bradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. I<PLoS Computational Biology> B<5:>e1000392.

* Clement M, Posada D, Crandall K (2000) TCS: a computer program to estimate gene genealogies. I<Molecular Ecology> B<9(10):>1657-60.

* Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. I<Nucleic Acids Research> B<32(5):>1792-7.

* Edgar RC (2010) Search and clustering orders of magnitude faster than BLAST. I<Bioinformatics> B<26(19):>2460-1.

* Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. I<Bioinformatics> B<27(16):>2194-200.

* Finn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. I<Nucleic Acids Research> B<39:>W29-W37.

* Huse SM, Huber JA, Morrison HG, Sogin ML, Welch DM (2007) Accuracy and quality of massively parallel DNA pyrosequencing. I<Genome Biology> B<8:>R143.

* Katoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. I<Nucleic Acids Research> B<33:>511-18.

* Katoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. I<Nucleic Acids Research> B<30:>3059-66.

* Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. I<Molecular Biology and Evolution> B<30(4):>772-80.

* Lassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. I<BMC Bioinformatics> B<6:>298

* Li W, Fu L, Niu B, Wu S, Wooley J (2012) Ultrafast clustering algorithms for metagenomic sequence analysis. I<Briefings in Bioinformatics> B<13(6):>656-68.

* Loytynoja A, Goldman N (2005) An algorithm for progressive multiple alignment of sequences with insertions. I<Proceedings of the National Academy of Sciences of USA> B<102:>10557-62.

* Lyras DP, Metzler (2014) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach". I<Submitted>. 

* Martin M (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. I<EMBnet.journal> B<17>.

* Morgulis A, Gertz EM, Schaffer AA, Agarwala R (2006) A fast and symmetric DUST implementation to mask low-complexity DNA sequences. I<Journal of Computational Biology> B<13:>1028-40.

* Russell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. I<BMC Bioinformatics> B<9:>306.

* Schmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. I<Bioinformatics> B<27:>863-4.

* Schwartz AS, Pachter L (2007) Multiple alignment by sequence annealing. I<Bioinformatics> B<23:>e24-9.

* Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Soeding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. I<Molecular Systems Biology> B<7:>539

* Stajich JE, Block D, Boulez K, et al. (2002) The Bioperl toolkit: Perl modules for the life sciences. I<Genome Research> B<12:>1611-8.

* Templeton AR, Crandall KA, Sing CF (1992) A cladistic analysis of phenotypic associations with haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram estimation. I<Genetics> B<132:>619-33.

=head1 AUTHOR

Enrique Gonzalez-Tortuero L<tortuero@bio.lmu.de>

=cut
