
import random as rnd

APPLOGGER=logging.getLogger("simphy-compress-dataset")
class NRSException(Exception):
	def __init__(self, expression, message, time):
		self.expression = expression
		self.message = message
		self.time= time

class SimPhyCompressDataset:
	startTime=None
import datetimeoutnameitertools,logging,os,sys, msatools,re, glob, gzip, shutil
import numpy as np
	path=""
	projectName=""
	inputprefix=""
	output=""
	outputFolderName=""
	nsize=-1
	numLociPerReplicate=[]
	numLociPerReplicateDigits=[]
	numReplicates=0
	numReplicatesDigits=0

	def __init__(self, args):
		self.startTime=datetime.datetime.now()
		self.endTime=None
		APPLOGGER.info(\
			"{0}".format(\
			"SCD started",\
		))
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

		output=os.path.abspath(args.output)
		outputFolderName=args.output_folder_name

		if (os.path.exists(os.path.join(output,outputFolderName))):
			listdir=os.listdir("{}".format(os.path.dirname(output)))
			counter=0
			for item in listdir:
				if outputFolderName in item:
					counter+=1
			if not counter == 0: outputFolderName+="_{0}".format(counter+1)
		self.output=os.path.join(os.path.dirname(output),outputFolderName)
		########################################################################
		# Generation of the output folder
		try:
			os.mkdir(self.output)
			APPLOGGER.info("Generating output folder:\t{}".format(self.output))
		except:
			APPLOGGER.info("Output folder ({0}) exists. ".format(self.output))

	def checkArgs(self):
		APPLOGGER.info("Checking arguments...")
		APPLOGGER.info("\tSimPhy...")
		simphydir=os.path.exists(self.path)
		########################################################################
		if simphydir:
			APPLOGGER.info("SimPhy folder exists:\t{0}".format(simphydir))
		else:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			ex="SimPhy folder does not exist."
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting.")
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
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
		self.getNumLociPerReplicate()
		self.getFilteredReplicates()
		########################################################################

	def getFilteredReplicates(self):
		self.filtered=[(idx+1) for idx in range(0,self.numLociPerReplicate) if not self.numLociPerReplicate[idx]==0]

	def getNumLociPerReplicate(self):
		for index in range(0, self.numReplicates):
			repID=index+1
			APPLOGGER.debug("Replicate {0}/{1} ".format(repID, self.numReplicates))
			curReplicatePath=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			)
			fileList=glob.glob("{0}/{1}_*_TRUE.fasta".format(curReplicatePath,self.inputprefix))
			prefixLoci=len(fileList)
			self.numLociPerReplicate[index]=prefixLoci
			self.numLociPerReplicateDigits[index]=len(str(prefixLoci))

	def iterateOverReplicate(self):
		APPLOGGER.debug("IterateOverReplicate")
		filtered=self.filterReplicatesMatchingIndPerSpeciesAndPloidy(self.ploidy)
		for index in range(0, self.numReplicates):
			if repID in self.filtered:
				APPLOGGER.debug("Replicate {0}/{1} ".format(repID, self.numReplicates))
				curReplicatePath=os.path.join(\
					self.path,\
					"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				)
				fileList=glob.glob("{0}/{1}_*.fasta".format(curReplicatePath,self.inputprefix))
				prefixLoci=len(fileList)
				APPLOGGER.info("Method chosen: {0}".format(self.method))
				print (self.numLociPerReplicateDigits)
				####################################################################
				self.concat(repID)

	def concat(repID):
		data_true=None
		data=None
		msa=None
		msa_true=None
		nsequence=""
		if (self.nsize > 0):
			nsequence=["N"]*self.nsize
		# missing n sequence - stop at concatenating.
		for locID in range(0, self.numLociPerReplicateDigits[repID-1]):
			inputfile_true=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				"{0}_{1:0{2}d}_TRUE.fasta".format(self.inputprefix,locID, self.numLociPerReplicateDigits),\
			)
			inputfile=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				"{0}_{1:0{2}d}.fasta".format(self.inputprefix,locID, self.numLociPerReplicateDigits),\
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

		self.writemsa(repID,"fasta",data)
		self.writemsa(repID,"true",data_true)

	def generateFolderStructure(self):
		out=os.path.join(self.output,self.outputFolderName,"fasta")
		out_true=os.path.join(self.output,self.outputFolderName,"true")
		out_trees=os.path.join(self.output,self.outputFolderName,"gtrees")
		out_rest=os.path.join(self.output,self.outputFolderName,"strees")
		try:
			os.makedirs(out)
		except:
			# error
			ex="Folder {} exist...".format(out)
			exc_type, exc_obj, exc_tb = sys.exc_info()
			self.generateException(ex,exc_type)
		try:
			os.makedirs(out_true)
		except:
			# error
			ex="Folder {} exist...".format(out_true)
			exc_type, exc_obj, exc_tb = sys.exc_info()
			self.generateException(ex,exc_type)

		try:
			os.makedirs(out_trees)
		except:
			# error
			ex="Folder {} exist...".format(out_trees)
			exc_type, exc_obj, exc_tb = sys.exc_info()
			self.generateException(ex,exc_type)

		try:
			os.makedirs(out_rest)
		except:
			# error
			ex="Folder {} exist...".format(out_rest)
			exc_type, exc_obj, exc_tb = sys.exc_info()
			self.generateException(ex,exc_type)

	def generateException(self,ex,exc_type):
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
			ex,exc_type,fname, exc_tb.tb_lineno,\
			"Please verify. Exiting."\
		)
		raise NRSException(False, message, datetime.datetime.now()-self.startTime)


	def writemsa(self,repID,folder,msa):
		APPLOGGER.info("Writing selected ST: {0}", repID)
		outname=os.path.join(\
			self.output,\
			folder,\
			"{0}_{1:0{2}d}_{3:0{4}d}.fasta.gz".format(\
				self.outputprefix,\
				repID,self.numReplicatesDigits\
			)\
		)
		with gzip.GzipFile(outname, 'wb') as zp:
			for item in msa.keys():
				zp.write(">{}\n{}\n".format(key,msa[key]))

	def concatSelectedLoci(self,index,locID,description,sequence):
		"""
		BEDFILE: replicateID startPOS endPOS locID
		"""
		repID=index+1
		APPLOGGER.info("Writing selected loci {1} from ST: {0}".format(repID,locID))
		outname=os.path.join(\
			self.output,\
			"{0}_{1:0{2}d}.fasta".format(\
				self.outputprefix,\
				repID,\
				self.numReplicatesDigits\
			)
		)
		newDes=">{0}:{1:0{2}d}".format(\
			self.projectName,\
			repID,\
			self.numReplicatesDigits\
		)
		nsequence="".join("N" for item in range(0,self.nsize))
		# I'm assuming that if the file does not exist it will be created
		fullseq="{0}{1}".format(sequence,str(nsequence))
		if os.path.exists(outname):
			with open(outname, 'a+') as f:
				f.seek(-1,2)
				if locID == self.numLociPerReplicate[repID-1]:
					f.write('\n'.encode())
				else:
					f.write(fullseq.encode())
		else:
			f=open(outname,"a+")
			f.write(">{0}\n{1}".format(description, fullseq))
			f.close()

	def run(self):
		"""
		Run process of the program.
		"""
		self.checkArgs()
		self.iterateOverReplicate()
		raise NRSException(True,"",datetime.datetime.now()-self.startTime)
