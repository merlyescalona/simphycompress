import datetime,logging,os,sys, msatools,re, glob, gzip, shutil
import numpy as np
import random as rnd

APPLOGGER=logging.getLogger("simphy-compress-dataset")
class NRSException(Exception):
	def __init__(self, expression, message, time):
		self.expression = expression
		self.message = message
		self.time= time

class SimPhyCompressDataset:
	startTime=None
	path=""
	projectName=""
	inputprefix=""
	nsize=-1
	numLociPerReplicate=[]
	numLociPerReplicateDigits=[]
	numReplicates=0
	numReplicatesDigits=0

	def __init__(self, args):
		self.startTime=datetime.datetime.now()
		self.endTime=None
		APPLOGGER.info("SCD Started")
		# Variable initialization
		self.inputprefix=args.input_prefix
		self.nsize=args.nsize
		########################################################################
		# Checking correctness of the given paths
		if (args.simphy_path[-1]=="/"):
			self.projectName=os.path.basename(args.simphy_path[0:-1])
		else:
			self.projectName=os.path.basename(args.simphy_path)
		self.path=os.path.abspath(args.simphy_path)


	def checkArgs(self):
		APPLOGGER.info("Checking arguments...")
		simphydir=os.path.exists(self.path)
		APPLOGGER.info("\tSimPhy: {}".format(simphydir))
		########################################################################
		if simphydir:
			APPLOGGER.info("SimPhy folder exists:\t{0}".format(simphydir))
		else:
			ex="SimPhy folder does not exist.\nPlease verify. Exiting."
			raise NRSException(False, ex, datetime.datetime.now()-self.startTime)
		fileList=os.listdir(self.path)
		for index in range(0,len(fileList)):
			fileList[index]=os.path.abspath(os.path.join(self.path,fileList[index]))
		APPLOGGER.info("\tIdentifying replicates...")
		# check how many of them are dirs
		for item in fileList:
			baseitem=os.path.basename(item)
			if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
				self.numReplicates=self.numReplicates+1
		self.numReplicatesDigits=len(str(self.numReplicates))
		########################################################################
		# check if at least one
		if not (self.numReplicates>0):
			ex="Number of replicates/folders:\t{0} [Required at least 1]".format(self.numReplicates>0)
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting."\
			)
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
		APPLOGGER.info("\tDone!")
		APPLOGGER.info("Number of replicates:\t{0}".format(self.numReplicates))
		########################################################################
		if self.nsize > -1:
			APPLOGGER.info("Reference sequence separation will be of {} Ns.".format(self.nsize))
		else:
			APPLOGGER.info("Sequences are concatenated without Ns sequences.")
		########################################################################
		self.numLociPerReplicate=[0 for item in range(0,self.numReplicates)]
		self.numLociPerReplicateDigits=[0 for item in range(0,self.numReplicates)]
		self.getNumLociPerReplicate()
		self.getFilteredReplicates()
		########################################################################

	def getFilteredReplicates(self):
		self.filtered=[(idx+1) for idx in range(0,self.numReplicates) \
		 	if not self.numLociPerReplicate[idx]==0]

	def getNumLociPerReplicate(self):
		for index in range(0, self.numReplicates):
			repID=index+1
			APPLOGGER.debug("Replicate {0}/{1} ".format(repID, self.numReplicates))
			curReplicatePath=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			)
			fileList=glob.glob("{0}/{1}_*_TRUE.fasta".format(\
				curReplicatePath,self.inputprefix))
			prefixLoci=len(fileList)
			APPLOGGER.info(index)
			self.numLociPerReplicate[index]=prefixLoci
			self.numLociPerReplicateDigits[index]=len(str(prefixLoci))

	def iterateOverReplicates(self):
		APPLOGGER.debug("IterateOverReplicate")
		for index in range(0, self.numReplicates):
			repID=index+1
			APPLOGGER.debug("Replicate {}/{}".format(repID, self.numReplicates))
			if repID in self.filtered:
				self.compressGeneTrees(repID)
				self.concatLoci(repID)
				self.removeUnzippedFasta(repID)
			else:
				self.compressGeneTrees(repID)


	def concatLoci(self,repID):
		data_true=None
		data=None
		msa=None
		msa_true=None
		nsequence=""
		if self.nsize > 0:
			nsequence="".join(["N"]*self.nsize)
		# missing n sequence - stop at concatenating.
		for locID in range(0, self.numLociPerReplicateDigits[repID-1]):
			inputfile_true=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				"{0}_{1:0{2}d}_TRUE.fasta".format(self.inputprefix,(locID+1), self.numLociPerReplicateDigits[repID-1])\
			)
			inputfile=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				"{0}_{1:0{2}d}.fasta".format(self.inputprefix,(locID+1), self.numLociPerReplicateDigits[repID-1])\
			)
			if locID==0:
				data=msatools.parseMSAFileWithDescriptions(inputfile)
				data_true=msatools.parseMSAFileWithDescriptions(inputfile_true)
			else:
				msa=msatools.parseMSAFileWithDescriptions(inputfile)
				msa_true=msatools.parseMSAFileWithDescriptions(inputfile_true)
				for key in msa.keys():
					data[key]="{}{}{}".format(data[key],nsequence,msa[key])
					data_true[key]="{}{}{}".format(data_true[key],nsequence,msa_true[key])
		self.writemsainplace(repID,"fasta",data)
		self.writemsainplace(repID,"true",data_true)

	def generateException(self,ex,exc_type,exc_obj, exc_tb):
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
			ex,exc_type,fname, exc_tb.tb_lineno,\
			"Please verify. Exiting."\
		)
		raise NRSException(False, message, datetime.datetime.now()-self.startTime)

	def removeUnzippedFasta(self,repID):
		outfolder=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID,self.numReplicatesDigits))
		APPLOGGER.info("(INPLACE) Revoving FASTA files from ST: {}  ({})".format(repID, outfolder))
		fileList=glob.glob(\
			"{0}/{1}_*_TRUE.fasta".format(\
				outfolder,\
				self.inputprefix\
			)\
		)
		for item in fileList:
			try:
				os.remove(item)
			except Exception as ex:
				APPLOGGER.warning("File cannot be removed ({})\n\t{}".format(item,ex.strerror))

		fileList=glob.glob(\
			"{0}/{1}_*.fasta".format(\
				outfolder,\
				self.inputprefix\
			)\
		)
		for item in fileList:
			try:
				os.remove(item)
			except Exception as ex:
				APPLOGGER.warning("File cannot be removed ({})\n\t{}".format(item,ex.strerror))

	def writemsainplace(self,repID,folder,msa):
		APPLOGGER.info("(INPLACE) Writing selected ST: {0}".format(repID))
		filename=""
		if folder == "true":
			filename="_{}".format(folder.upper())
		outname=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID,self.numReplicatesDigits),\
			"{0}_{1:0{2}d}{3}.fasta.gz".format(\
				self.inputprefix,\
				repID,\
				self.numReplicatesDigits,\
				filename\
			)\
		)
		with gzip.GzipFile(outname, 'wb') as zp:
			for key in msa.keys():
				zp.write(">{}\n{}\n".format(key,msa[key]))

	def compressGeneTrees(self,repID):
		outfolder=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID,self.numReplicatesDigits))
		APPLOGGER.info("(INPLACE) Compressing gene tree files from ST: {}  ({})".format(repID, outfolder))
		fileList=glob.glob(\
			"{0}/g_trees*.trees".format(\
				outfolder\
			)\
		)
		fileList=sorted(fileList)
		APPLOGGER.info("Len filelist: {}".format(len(fileList)))
		outname=os.path.join(outfolder,"gtrees.trees.gz")
		with gzip.GzipFile(outname, 'wb') as zp:
			for gtreefile in fileList:
				with open(gtreefile, "rb") as f:
					tree=f.readline()
					name,_=os.path.splitext(os.path.basename(gtreefile))
					# APPLOGGER.info("{}\t{}".format(name,tree))
					zp.write("{}\t{}".format(name,tree))
		for item in fileList:
			try:
				os.remove(item)
			except Exception as ex:
				APPLOGGER.warning("File cannot be removed ({})\n\t{}".format(item,ex.strerror))

	def run(self):
		"""
		Run process of the program.
		"""
		self.checkArgs()
		self.iterateOverReplicates()
		raise NRSException(True,"",datetime.datetime.now()-self.startTime)
