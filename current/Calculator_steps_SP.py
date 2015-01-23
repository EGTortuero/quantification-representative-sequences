#! /usr/bin/python3

### Importing Python modules

# Core modules
from functools import reduce
import glob, itertools, multiprocessing, operator, os, platform, re, string, subprocess, sys

# External modules
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from scipy.integrate import quad
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

# Defining internal functions

# Preliminar functions
def numeratorintegrand(teoq, j, m, b):
	return teoq * 2 * teoq ** (int(j) - 1) * (1 - teoq) ** (2 * int(m) + 1) * (1 - teoq / int(b)) * (2 - teoq * ((int(b) + 1) / int(b) )) ** (int(j) - 1) * (1 - 2 * teoq * (1 - teoq / int(b)))

def denominatorintegrand(teoq, j, m, b):
	return 2 * teoq ** (int(j) - 1) * (1 - teoq) ** (2 * int(m) + 1) * (1 - teoq / int(b)) * (2 - teoq * ((int(b) + 1) / int(b))) ** (int(j) - 1) * (1 - 2 * teoq * (1 - teoq / int(b)))

def prod(iterable):
	return reduce(operator.mul, iterable, 1)

## Print header of the program
print('''
#########################################
#      CALCULATOR STEPS v. 0.0.0a     	#
#---------------------------------------#
#   Author: Enrique Gonzalez-Tortuero   #
#########################################
''')

### Interpreting the parameters
nrparameters = len(sys.argv)
if nrparameters > 1:
	fastafile = sys.argv[1]
	limitid = sys.argv[2].replace(r'-limitsteps=','')
	b = sys.argv[3].replace(r'-model=','')
else:
	sys.exit("Calculator_steps_SP.py [fasta file] -limitsteps=[0.00-1.00] -model=[1/2/3]\nWhere 1 is K2P model and 3 is JC model according to Templeton et al (1992).")

# Knowing the alignment length
if os.path.isfile(fastafile):
	origfasta = open(fastafile, "rU")
else:
	sys.exit("%s can not be found\n" % fastafile)
for record in AlignIO.parse(origfasta, "fasta"):
	lengthDNA = record.get_alignment_length()
	print("\n%s have %i bp length." % (fastafile, lengthDNA))
origfasta.close()

## Calculate probability of parsimony until a given limit (by default, limitid=0.95).
rpy2.robjects.numpy2ri.activate()
Pind = []
Ppars = 0
CorPpars = 0
                
if limitid != None or limitid != "":

	# Calculating parsimony probability until a given limit
	print("\nCalculating the number of steps needed until your cut-off\n---------------------------------------------------------\n")
	print("Nr_steps\tPars_prob")
	for j in range(1, lengthDNA):
		m = lengthDNA - j # m is the number of invariant sites, and j is the number of changes in DNA
		A = quad(numeratorintegrand, 0, 1, args = (j, m, int(b)))
		B = quad(denominatorintegrand, 0, 1, args = (j, m, int(b)))
		if Pind != []:
			Pind.extend([1 - A[0] / B[0]])
			Pitmanq.extend([A[0] / B[0]])
		else:
			Pind = [1 - A[0] / B[0]]
			Pitmanq = [A[0] / B[0]]
		if prod(Pind[0:j]) < float(limitid):
			break
		else:
			Ppars = prod(Pind[0:j])
			limitstep = j
			print("%i\t%.6f" % (limitstep, Ppars))

else:
	sys.exit("Calculator_steps_SP.py [fasta file] -limitsteps=[0.00-1.00] -model=[1/2/3]\n Where 1 is K2P model and 3 is JC model according to Templeton et al (1992).")
