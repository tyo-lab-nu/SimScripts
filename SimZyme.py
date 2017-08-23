__author__ = 'Dante'

import DBParse
from optparse import OptionParser

if __name__ == '__main__':

	usage = "Usage: %prog [options] substrateSmiles outputname [ecNumbers]"
	parser = OptionParser(usage=usage)
	parser.add_option('-b', '--DBDumpFile', dest='DB', default='DB_download.txt', help='DB download file')
	parser.add_option('-e', '--ecFlag', dest='ec', action='store_false', default=True, help='Add to use mapping results')
	parser.add_option('-i', '--inchiDumpFile', dest='inchi', default='inchiFile.json', help=' default inchi dump json from BMK-React')
	parser.add_option('-s', '--smilesDumpFile', dest='smiles', default='DBCompFile.json', help='Smiles dump json')
	parser.add_option('-x', '--severFlag', dest='sever', action='store_true', default=False, help='Add to separate the DB file')
	(options, args) = parser.parse_args()

	results_dump = DBParse.SimZyme(options.DB, args[0], options.smiles, options.inchi, args[1], args[2:], options.ec, options.sever)