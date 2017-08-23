__author__ = 'Dante'

import sys
import os
from optparse import OptionParser
import numpy
import openbabel as ob
from pybel import Outputfile
from pybel import readfile

def get_Molfiles(molFileLocation, startCompFile):
	'''Grab all molfiles from a folder and separate them from the cofactors.

	Arguments: molFileLocation: folder in question
			   startCompFile: file listing start compounds.

	Returns: startComp: start compounds, sans cofactors
		     genMols: generated molfiles
	'''

	compList = os.listdir(molFileLocation)
	intList = [int(x.strip('.mol')) for x in compList].sort()
	allMolList = [os.path.join(molFileLocation, x + '.mol') for x in intList]

	#Read start compounds file for cofactors and non-cofactors.
	f = open(startCompFile, 'r')
	startComp = []
	cof = []
	for x in f:
		if 'Cofactor' in x:
			cof.append(x.strip('\n'))
		elif '.mol' in x:
			startComp.append(x.strip('\n'))
		else:
			pass

	f.close()

	m = len(allMolList)
	n = len(cof)

	#The `new` molfiles are the ones not in the start compounds file, so we segregate them here.
	if m > 0:
		genMols = allMolList[m:n]
	elif m == 0:
		genMols = allMolList[n:]
	else:
		sys.exit('Invalid input compounds file.')

	return genMols, startComp

def make_LibrariesFp(genMols, startComp, target, outName, fpForm):
	'''Messy way to make an sdf and subsequent dictionary of TC values.
	   Don't worry about most of the arguments, they end up defaulted.

	Arguments:  genMols: list of generated molfiles
	            startComp: list of starting molfiles
	            target: molfile of target compound
	            outname: name of library
	            fpForm: fignerprint type

	Returns:  baseTc: scalar TC value that is the basis for elimination
		      molTcs: Tc value fo generated compound with the target.
	'''

	genMolLibrary = Outputfile("sdf", "%s_MolLibrary.sdf" % outName, True)
	startLibrary = Outputfile("sdf", "%s_StartLibrary.sdf" % outName, True)

	#Compress into sdf file.
	for x in genMols:
		for y in readfile('mol', x):
			genMolLibrary.write(y)
	genMolLibrary.close()

	startCompFiles = [os.path.join('Molfiles', x.rstrip('\n')) for x in startComp]
	for x in startCompFiles:
		for y in readfile('mol', x):
			startLibrary.write(y)
	startLibrary.close()

	targetFile = readfile('mol', target).next()

	#Suppress error messages for neatness.
	error_log_check = hasattr(ob, "obErrorLog")
	if error_log_check:
		ob.obErrorLog.ClearLog()
		ob.obErrorLog.SetOutputLevel(-1)

	#Generate fingerprints and calculate TCs...then make them into a dictionary.
	molFps = [x.calcfp(fpForm) for x in readfile("sdf", "%s_MolLibrary.sdf" % outName)]
	startFps = [x.calcfp(fpForm) for x in readfile("sdf", "%s_StartLibrary.sdf" % outName)]
	targetFp = targetFile.calcfp(fpForm)

	startTcs = [x | targetFp for x in startFps]
	baseTc = numpy.average(startTcs)
	molTcs = {genMols[i]: targetFp | x for i, x in enumerate(molFps)}

	return baseTc, molTcs

def remove_LowTc(tcBasis, tcDict, t):
	'''Elimination function.

	Arguments: tcBasis: base TC
	           tcDict: Dictionary of TCs
	           t: tolerance
	'''

	i = 0

	for k, v in tcDict.iteritems():
		if v <= t * tcBasis:
			os.remove(k)
			i += 1

def sim_index_main(location, startCompFile, target, outName, fpForm, t):
	'''Main script. For easy wrapper calling.

	Arguments: See above functions.
	'''


	genMols, startComp = get_Molfiles(location, startCompFile)
	basis, molTcs = make_LibrariesFp(genMols, startComp, target, outName, fpForm)
	remove_LowTc(basis, molTcs,t)

if __name__ == '__main__':

	usage = "Usage: %prog [options] location/ startcompfile.dat targetmoleculepath.mol"
	parser = OptionParser(usage=usage)
	parser.add_option('-f', '--fingerprint', dest='finger', default='FP4', help='Fingerprint format, choose FP2, FP3, or FP4')
	parser.add_option('-n', '--nameoutput', dest='name', default='Net_Gen', help='Name of output sdf libraries')
	parser.add_option('-t', '--tolerance', dest='tol', default=1, help='Tolerance for compound TC removal')
	(options, args) = parser.parse_args()

	if options.finger not in ['FP4', 'FP2', 'fp2', 'fp4', 'FP3', 'fp3']:
		sys.exit('Fingerprint type is invalid\nValid types are FP2, FP3, or FP4')

	if options.tol < 0:
		sys.exit('Tolerance is negative\nChoose a value greater than 0')

	if options.tol > 1:
		print 'Warning! Tolerance given could result in a loss of a generation.'

	sim_index_main(options.loc, args[0], args[1], options.name, options.finger, float(options.tol))















