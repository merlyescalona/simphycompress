import argparse,datetime,logging,os,sys, refselector, msatools
import loggingformatter as lf
import numpy as np
import random as rnd
################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=0
FIX_VERSION=1
PROGRAM_NAME="simphy-compress-dataset"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
INSTITUTION="University of Vigo, Spain."
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
LINE="--------------------------------------------------------------------------------"
################################################################################
# python LociReferenceSelection.py -p <prefix> -SF <simphy_path> -o outout -m method
################################################################################
# Logger init
################################################################################
ch = logging.StreamHandler()
loggerFormatter=lf.MELoggingFormatter(\
	fmt="%(asctime)s - %(levelname)s:\t%(message)s",\
	datefmt="%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
ch.setLevel(logging.NOTSET)
APPLOGGER=logging.getLogger("simphy-compress-dataset")
APPLOGGER.addHandler(ch)
################################################################################
def createLogFile():
	formatString=""
	if platform.system()=="Darwin":
		formatString="%(asctime)s - %(levelname)s (%(module)s:%(lineno)d):\t%(message)s"
	else:
		formatString="%(asctime)s - %(levelname)s (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s"
	fh=logging.FileHandler(\
		"{0}/{2}.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
			os.getcwd(),\
			datetime.datetime.now(),\
			PROGRAM_NAME[0:-3].upper()\
			)\
		)
	fh.setLevel(logging.DEBUG)
	formatter=logging.Formatter(formatString)
	fh.setFormatter(formatter)
	APPLOGGER.addHandler(fh)
################################################################################
# Handling parameters
################################################################################
def handlingParameters():
	parser = argparse.ArgumentParser(\
        prog="{0} (v.{1}.{2}.{3})".format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        description=\
'''
\033[1m
================================================================================
SimPhy compress dataset
================================================================================
\033[0m
Description:
============

This program allows to compress the number of files and file sizes of a SimPhy run.

Assumptions:
============

- We are working under a SimPhy simulation scenario. Meaning, it follows hierarchical
SimPhy's folder structure and sequence labeling.
- SimPhy folder contains all filles from a SimPhy run (dbs, params...s)

General:
========

Output:
=======

		''',\
			epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
			add_help=False
		)
	requiredGroup= parser.add_argument_group('Required arguments')
	requiredGroup.add_argument('-s','--simphy-path',metavar='<path>', type=str,\
		help='Path of the SimPhy folder.', required=True)
	requiredGroup.add_argument('-ip','--input-prefix', metavar='<input_prefix>', type=str,\
		help='Prefix of the FASTA filenames.', required=True)
	requiredGroup.add_argument('-o','--output', metavar="<output_path>",type=str,\
		help="Path where output will be written.", required=True)
	requiredGroup.add_argument('-ofn','--output-folder-name', metavar="<folder_name>",type=str,\
		help="Name of the output folder. If a folder with such name exist in output_path/, folder name "+\
		"will have a suffix (a number) according to the number of existent folders + 1 (output_path/output_folder_name_n)", required=True)
	optionalGroup= parser.add_argument_group('Optional arguments')
	optionalGroup.add_argument('-n','--nsize',metavar='<N_seq_size>', type=int,\
		default=-1,\
		help="Number of N's that will be introduced to separate the sequences selected. "+\
            "If the parameter is not set, the output file per replicate will be a multiple alignment sequence file, "+\
            "otherwise, the output will be a single sequence file per replicate consisting of a concatenation  "+\
            "of the reference sequences selected separated with as many N's as set for this parameter.")
	optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
		choices=LOG_LEVEL_CHOICES, default="INFO",\
		help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(LOG_LEVEL_CHOICES,LOG_LEVEL_CHOICES[1]))
	informationGroup= parser.add_argument_group('Information arguments')
	informationGroup.add_argument('-v', '--version',\
		action='version',\
		version='Version {0}.{1}.{2}'.format(VERSION,MIN_VERSION,FIX_VERSION),\
		help="Show program's version number and exit")
	informationGroup.add_argument('-h', '--help',\
		action='store_true',\
		help="Show this help message and exit")
	try:
		tmpArgs = parser.parse_args()
	except:
		sys.stdout.write("\n\033[1m{}\033[0m\n".format(LINE))
		APPLOGGER.error("Something happened while parsing the arguments.")
		APPLOGGER.error("Please verify. Exiting.\n{}".format(LINE))

		parser.print_help()
		sys.exit(-1)
	return tmpArgs

################################################################################
# MAIN
################################################################################
def main():
	try:
		cmdArgs = handlingParameters()
		APPLOGGER.setLevel(cmdArgs.log.upper())
		APPLOGGER.debug("Args. introduced: {}".format(cmdArgs))
		prog = refselector.ReferenceSelection(cmdArgs)
		prog.run()
	except refselector.NRSException as ex:
	    if ex.expression:
	        APPLOGGER.info("REFSELECTOR finished properly.")
	        APPLOGGER.info("Elapsed time (ETA):\t{0}".format(ex.time))
	        APPLOGGER.info("Ending at:\t{0}".format(datetime.datetime.now().strftime("%a, %b %d %Y. %I:%M:%S %p")))
	        sys.exit()
	    else:
	        APPLOGGER.error(ex.message)
	        APPLOGGER.error("Elapsed time (ETA):\t{0}".format(ex.time))
	        APPLOGGER.error("Ending at:\t{0}".format(datetime.datetime.now().strftime("%a, %b %d %Y. %I:%M:%S %p")))
	        sys.exit(-1)
	except KeyboardInterrupt:
	    APPLOGGER.error("{0}{1}\nProgram has been interrupted.{2}\nPlease run again for the expected outcome.\n{3}\n".format("\033[91m","\033[1m","\033[0m",LINE))
	    sys.exit(-1)

if __name__=="__main__":
	main()
