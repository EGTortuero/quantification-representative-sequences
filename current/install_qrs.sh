#! /bin/bash

# How to install QRS
# Tested in Ubuntu 14.10.

echo "Automatic installation of QRS 1.4.0"
echo "-----------------------------------"
mkdir qrs
cd qrs
mkdir programs
cd programs

echo "First step: installing pre-requisites (libraries and languages)"
sudo apt-get install r-base-dev # Installation of R headers according to CRAN (with this command, R should be installed automatically).
sudo apt-get install python3-pip # Installation of PIP (Python libraries installer)
sudo apt-get install libjson-perl libcairo-perl # JSON and Cairo libraries
sudo apt-get install libargtable2-dev # Libraries needed for Clustal Omega
sudo apt-get install csh
wget http://cran.r-project.org/src/contrib/ade4_1.6-2.tar.gz
sudo R CMD INSTALL ade4_1.6-2.tar.gz # NOTE: This library is needed in order to avoid matters with SeqinR

echo "Second step: installing python modules"

# Installation of scipy and numpy (OBSERVATION: numpy is automatically installed because scipy requires numpy)
sudo apt-get install python3-scipy

# Installation of rpy2, psutil and cutadapt via python3
sudo pip3 install rpy2
sudo pip3 install psutil
sudo pip3 install cutadapt

# Installation of biopython via setuptools
sudo easy_install3 -f http://biopython.org/DIST/ biopython

echo "Third step: installing R libraries"

# Installation of ape
wget http://cran.r-project.org/src/contrib/ape_3.2.tar.gz
sudo R CMD INSTALL ape_3.2.tar.gz

# Installation of seqinr
wget http://cran.r-project.org/src/contrib/seqinr_3.1-3.tar.gz
sudo R CMD INSTALL seqinr_3.1-3.tar.gz

echo "Fourth step: installing external programs"
echo "NOTE: You must to request a copy of USEARCG from the Author webpage"
# mv usearch8.0_linux32 usearch
# sudo cp usearch /usr/local/bin/

echo "PrinSeq:"
wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
tar xvfz prinseq-lite-0.20.4.tar.gz
cd prinseq-lite-0.20.4
chmod +x prinseq-graphs-noPCA.pl prinseq-lite.pl
sudo cp prinseq-graphs-noPCA.pl /usr/local/bin/
sudo cp prinseq-lite.pl /usr/local/bin/
cd ..

echo "PEAR:"
wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-src.tar.gz
tar xvfz pear-0.9.6-src.tar.gz
cd pear-0.9.6-src
./configure
make
sudo make install
cd ..

echo "HECTOR:"
wget http://sourceforge.net/projects/hector454/files/HECTOR/hector.tar.gz
tar zxvf hector.tar.gz
make clean
make
sudo cp hector /usr/local/bin/
cd ..

echo "HMMER:"
wget http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-x86_64.tar.gz
tar zvfx hmmer-3.1b1-linux-intel-x86_64.tar.gz
cd hmmer-3.1b1-linux-intel-x86_64
./configure
make
sudo make install
cd ..

echo "Clustal Omega:"
wget http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz
tar zxvf clustal-omega-1.2.1.tar.gz
cd clustal-omega.1.2.1.tar.gz
./configure
make
sudo make install
cd ..

echo "FSA and MUMMER:"
# Installing MUMMER (NOTE: Requirement for long sequences)
wget http://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz
tar zvfx MUMmer3.23.tar.gz
cd MUMmer3.23
make check
sudo make install
sudo cp annotate combineMUMs delta-filter dnadiff exact-tandems gaps mapview mgaps mummer mummerplot nucmer nucmer2xfig promer repeat-match run-mummer1 run-mummer3 show-aligns show-coords show-diff show-snps show-tiling /usr/local/bin/
cd ..

# Installing FSA
wget http://sourceforge.net/projects/fsa/files/fsa-1.15.9.tar.gz
tar xvzf fsa-1.15.9.tar.gz 
cd fsa-1.15.9 
./configure 
make 
sudo make install 
cd ..

echo "GramAlign:"
wget http://bioinfo.unl.edu/downloads/GramAlign3_00.zip
unzip GramAlign3_00.zip
cd GramAlign3_00/src
make
sudo cp GramAlign /usr/local/bin
cd ../..

echo "KAlign:"
wget http://msa.sbc.su.se/downloads/kalign/current.tar.gz
mkdir Kalign; cp current.tar.gz Kalign/; cd Kalign
tar zxvf current.tar.gz
./configure
make
sudo make install
cd ..

echo "MAFFT:"
wget http://mafft.cbrc.jp/alignment/software/mafft-7.215-with-extensions-src.tgz
tar zxvf mafft-7.215-src.tgz
cd mafft-x.x/core/
make clean
make
sudo make install
cd ../extensions/
make clean
make
sudo make install
cd ../..

echo "MUSCLE:"
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz
tar zxvf muscle3.8.31_src.tar.gz
cd muscle3.8.31/src
make
sudo cp muscle /usr/local/bin/
cd ../..

echo "Opal:"
wget http://opal.cs.arizona.edu/opal.tgz
tar zxvf opal.tgz
cd opal_2.1.3
sudo cp opal Opal.jar predict_structure.pl /usr/local/bin/
cd ..

echo "PicXAA:"
wget http://www.ece.tamu.edu/~bjyoon/picxaa/PicXAA-v1.03.tar.gz
tar zxvf PicXAA-v1.03.tar.gz
cd PicXAA-v1.03
make clean
make
sudo cp picxaa /usr/local/bin
cd ..

echo "PRANK:"
sudo apt-get install exonerate # NOTE: PRANK needs Exonerate for anchoring sequences. Highly recommended!
wget http://wasabiapp.org/download/prank/prank.source.140603.tgz
tar zvxf prank.source.140603.tgz
cd prank-msa/src
make
sudo cp prank /usr/local/bin
cd ../..

echo "ReformAlign:"
wget http://evol.bio.lmu.de/_statgen/software/reformalign/ReformAlign.zip
unzip ReformAlign.zip
cd WEBSITE/Source\ Code\ and\ Makefiles/
cp Makefile.unix Makefile
make
sudo cp bin/Release/ReformAlign /usr/local/bin
cd ../..

echo "Fifth: Final step"

# Removing installation folder
cd ..
rm -rf programs

# Executing qrs_1.4
chmod +x qrs_1.4.py

echo "It is DONE"

