#!/usr/bin/env python3


# Library installation notes:
# plotly:
#   pip3  install chart_studio
#   pip3 install plotly
#   or, for updates
#   pip3 install plotly --upgrade
# scikit-learn:
# ref: https://www.kaggle.com/c/titanic/discussion/6801
#   python3 -m pip install scikit-learn
#   this might be necessary: pip3 install scipy
#   or try running: python3 -m pip install scikit-learn
# or use use if trouble installing on linux: python3 -m  pip install scikit-learn --user

#
# scipy regression demo: https://scikit-learn.org/stable/auto_examples/linear_model/plot_ols.html#sphx-glr-auto-examples-linear-model-plot-ols-py


import os, sys, math
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
#import plotly.plotly as py
import chart_studio.plotly as py

import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np



DATE = "24 July 2020"
VERSION = "2_ii"
AUTHOR = " Oliver Bonham-Carter"
AUTHORMAIL = "obonhamcarter@allegheny.edu"
#OUTPUT_DIR = "/tmp/0out/" # all results are saved in this local directory
OUTPUT_DIR = "0out/" # all results are saved in this local directory
INPUT_DIR = "data/"

# the file containing the names of genes to be used for normalizing
NORM_FILE = "normNames_i.csv"

# the file containing the names of the genes that we are studying.
PICK_THESE = "pickThese.csv" # contains the genes to study in the datasets. note: this file must include all data files and those used for normalizing.
THRESH = "thresh.csv" # contains a list of thresholds of r-squared values to study in heatmaps.

# the below line is to exclude particular files
IGNORE_FILES_list = [".DS_Store", "MANIFEST.txt",NORM_FILE, "~lock","annotations.txt",".gz",".html",PICK_THESE]


def help():
	h_str = "   "+DATE+" | version: "+VERSION+" |"+AUTHOR+" | "+AUTHORMAIL
	print("  "+len(h_str) * "-")
	print(h_str)
	print("  "+len(h_str) * "-")

	print("\n\tThe GeneExPy program to perform linear regression over GDP datasets.")

	print("""\n\tLibrary installation notes:
	plotly:
	pip3 install plotly, or try running python3 -m pip install scikit-learn
	scikit-learn:
	python3 -m pip install scikit-learn, maybe necessary: pip3 install scipy
	""")
	print("\t+ \U0001f600  USAGE: programName <any key to launch>")
	print("\t+ INPUT directory: (your data files are here)     : ",INPUT_DIR)
	print("\t+ OUTPUT directory: (your output is placed here)  : ",OUTPUT_DIR)
	print("\n\t+ Note:\n\t   Use parameter <<heatmap>> or <rsqu>>> to ensure that\n\t   Plotly and the statistical libraries have been correctly installed.")
	print("\n\t Note: the data directory cannot handle subdirectories holding data. Please \n\tplace the text files into this data directory without using subdirectories.")
#end of help()


class Wrangler:
	#""" Class to wrangle the data: to convert files to usable data for analysis"""
	def __init__(self):
		""" initiation method for Wrangler class"""
		#print("   Wrangler Class __init__()")
		self.file_list = [] # holds each file and diretory
		self.lower_list =  [g.lower() for g in IGNORE_FILES_list]
		#self.ensNums_list = [] # holds the ensNums for building a mtrix
		#self.exp_list = [] # holds the expression for each dataset
		self.raw_dic = {} # the matrix in a dictionary. Key is dataset, value is list of expressions in order of file
		# Note:
		#self.raw_dic[ds_str][0] # contains the ensNumbers for this dataset
		#self.raw_dic[ds_str][1] # contains the raw expressions for this dataset
		#self.raw_dic[ds_str][2] # contains the rawlogs of expressions for this dataset
		#self.raw_dic[ds_str][3] # contains the AVGg1 normalizing
		#self.raw_dic[ds_str][4] # contains the AVGg2 normalizing
		#self.raw_dic[ds_str][5] # contains the AVGg3 normalizing
		#self.raw_dic[ds_str][6] # contains the tubb normalizing
		#self.raw_dic[ds_str][7] # contains the tuba1a normalizing
		self.picker_dic = {} # a dic of genes to study from the PICK_THESE file
		self.normNames_dic = {} # a listing of the names of genes to be used for averaging and then normalizing
		self.groupGene_dic = {} # contains the names of the genes of each group. exp of these used for avgs.
		self.thresh_list = [] # listing of thresholds to capture r-squared values for study.
		self.rSquMat_dic = {} # dic to contain the r-squared matrices. key is name.html, value is [[x_list], [y_list], [z_list]]
		self.threshold_dic = {} # dic to contain the threshold specific heatmaps. prepared by filterMatrix()
	#end of __init__()


	def openTextFile(self, inFile):
		#print("openTextFile(): ",inFile)
		#read a list and then restun a dic
		try: # is there a file?
			#data = open(inFile).read().lower() # return a string
			data = open(inFile).readlines()
			return data
		except IOError:
			print("No file found... \"",inFile,"\" Exiting")
			sys.exit(1)
	#end of openTextFile()


	def isFileInIgnoreList(self, fName_str):
		"""function to determine whether a file is to be ignored"""
		#lower_list =  [g.lower() for g in self.IGNORE_FILES_list]
		fName_str = fName_str.lower()
		#print(" fName_str :" ,fName_str)
		#        print("self.lower_list",self.lower_list)
		for i in self.lower_list:
			if i in fName_str:
			#print(i,"found in ",fName_str)
				return True # ignore the file
		return False # load the file, not in the list
#end of isFileInIgnore()


	def getFileListing(self,corpusDir):
		""" method to grab all files not on the ignore list """
		#self.file_list = [] # holds each file and diretory
		for root, dirs, files in os.walk(corpusDir):
			for file1 in files:
				if self.isFileInIgnoreList(file1) == False:
					#print("loading ",file1)
					dataFile = os.path.join(root, file1)
					self.file_list.append(dataFile)
				else:
					print("\t- Ignoring file at this step: ",file1)
		return self.file_list
	#end of getFileListing


	def getParams(self):
		"""Method to open the picker file (PICK_THESE) to find the genes to study and the thresholds file (THRESH) for an r-squared focus. """
		self.picker_dic = {} # a dic of genes to study from the PICK_THESE file
		self.thresh_list = [] # listing of thresholds to capture r-squared values for study.

# pickThese file format:
# UBE2V2,ENSG00000169139.10
# FAAP20 (C1orf86),ENSG00000162585.15
# DNAJA4,ENSG00000140403.11
# PSMD4,ENSG00000159352.14
# POLE3,ENSG00000148229.11


# include the PICK_THESE genes to use.
# file format
# ABR,ENSG00000159842.13
# ACTR8,ENSG00000113812.12
# APLF,ENSG00000169621.8
# APTX,ENSG00000137074.17

		d = open(PICK_THESE)
		for i in d:
			isplit_list = i.strip().split(",") # positions 0 and 1 for ensNums and expressions, resp
			self.picker_dic[isplit_list[1]] = isplit_list[0].strip()#replace(" ","")# ensnum (key) gene name (value)
#        printer(self.picker_dic)


# include the NORM_FILE genes to use.
# file format
# WDR77,ENSG00000116455.12,1
# USP39,ENSG00000168883.18,2
# CDC5L,ENSG00000096401.7,2

		d = open(NORM_FILE)
		for i in d:
			isplit_list = i.strip().split(",") # positions 0 and 1 for ensNums and expressions, resp
			self.picker_dic[isplit_list[1]] = isplit_list[0].strip()#replace(" ","")# ensnum (key) gene name (value)
#        printer(self.picker_dic)

# threshold values; genes having less than or equal to these values
# file format:
# 0.1
# 0.2
# 0.3
		d = open(THRESH)
		for i in d: # each line of the file
			val_float = float(i)
			self.thresh_list.append(val_float)
#        printer(self.thresh_list)

# end of getParams()


	def getRawMatrix(self):
		"""Method to load each file in the file_list and then create a huge dictionary to make matrix for working"""
		#printer(self.picker_dic)
		for f in self.file_list:# for each file
			ensNums_list = [] # contains the ensemble numbers from the current file
			exp_list = [] # contains the expression values from current file
			rawlogs_list = [] # contains the log of the exp. 0's are placed instead of log(0)
			m_list = [] # contains all the above lists as a list.
			d = self.openTextFile(f)
			#            counter = 0
			for i in d: # in the file itself
				m_list = []
				isplit_list = i.split() # positions 0 and 1 for ensNums and expressions, resp
				#print("\ti in d: isplit :",isplit_list, "file: ",f)
				# prints ['ENSG00000153561.11', '16.1267035586'] file:  data/d578e27f-537c-4aaa-8903-6ffe68346276.FPKM.txt

# raw data
				if isplit_list[0] in self.picker_dic: #if ensNum in self.picker_dic, then keep
					#print(isplit_list[0],"in the dictionary: RAW")
					#                    counter += 1
					ensNums_list.append(isplit_list[0])
					exp_list.append(isplit_list[1])
#                print("total :",counter)
# rawlogs
					try: #replace the counter list with a log-normed of raw data
						rawlogs_list.append(math.log(float(isplit_list[1]),math.exp(1)))
					except ValueError:
						rawlogs_list.append(0)
				#counter  += 1
			m_list = [ensNums_list, exp_list, rawlogs_list]
			#print(m_list[0:10]) # show what the data looks like...
			#print("FFFFFFFF:", f)
			#ff = f.replace(" ","") # remove the spaces in the filename, good for dictionary keys
			self.raw_dic[f] = m_list
		#print("keys (files) ")
		#print(self.raw_dic.keys())

#end of getRawMatrix(file_list)

	def getNormNamesMatrix(self, normFilename):
		"""Method to load the file containing the names of files to use for normalizing factor creation. Returns a matrix of these gene names in file. """
#file format:
# humanGene	EnsNum	Group
# WDR77	ENSG00000116455.12	1
# USP39	ENSG00000168883.18	2
# CDC5L	ENSG00000096401.7	2
# CASC3	ENSG00000108349.13	1

		counter_list = [] # contains the position in the file
		counter = 0
		humanGene_list = [] # contains the human gene names
		ensNums_list = [] # contains the ensemble numbers from the current file
		group_list = [] # contains the group number for each gene. The group is to deterine which set the gene is to be placed for averaging
		m_list = [] # contains all the above lists as a list.
		d = self.openTextFile(normFilename)
		for i in d: # in the file itself
			isplit_list = i.split() # positions 0 and 1 for ensNums and expressions, resp
			#print("\tisplit :",isplit_list,type(isplit_list))
			if len(isplit_list) > 1:
				humanGene_list.append(isplit_list[0])
				ensNums_list.append(isplit_list[1])
				group_list.append(isplit_list[2])
				counter_list.append(counter)
				counter  += 1
			if len(isplit_list) == 1:
				headers_list = isplit_list[0].split(",")
				#print("\theaders_list :", headers_list)
				humanGene_list.append(headers_list[0])
				ensNums_list.append(headers_list[1])
				group_list.append(headers_list[2])
				counter_list.append(counter)
				counter += 1
		m_list = [humanGene_list, ensNums_list, group_list, counter_list]
		#print(m_list[0:10])
		self.normNames_dic[normFilename] = m_list
		#print(self.normNames_dic.keys())

#end of getNormNamesMatrix(NormFileName)


# TODO: compare the dataset names across the whole set of files to make sure that all expressions are for the same ensNum.
# Not sure that this is necessary?
	def compareGeneOrder(self):
		print("CompareGeneOrder()")
	#end of compareGeneOrder()


	def getGroups(self):
		""" Method to determine the number of groups. Creates self.group_set to contain groups and self.gene_dic to hold the genes from the group."""
		#print("\traw data :",self.raw_dic.keys())
		#print(self.raw_dic.keys())
		#print("normNames")
		#print(self.normNames_dic)

		# make a dictionary of lists (values) with key (groups)
		group_dic = {}
		#        print(self.normNames_dic[NORM_FILE])

		main_list = self.normNames_dic[NORM_FILE]
		#        print("&&&&", main_list)
		#print("\n0",main_list[0]) # humanGene
		#print("\n1",main_list[1]) # ensNum
		#print("\n2",main_list[2]) # group

		# How many groups are there?
		self.group_set = set() # contains element to represent each group
		for i in main_list[2]:
			self.group_set.add(i)
		self.group_set.discard("Group") # note: discard removes the member element but does nothing if element is not in set
		#print("self.group_set :",self.group_set)

		# go through each group, pull ensNums of each to make a list
		self.groupGene_dic = {}
		my_list = []
		myGeneGroup_list = []
		for i in self.group_set:
			myGeneGroup_list = []
			#print( "\tGroup: ",i)
			for j in range(len(main_list[2])): # look at each position in the group list
				#print(main_list[2][j], type(main_list[2][j]))
				if main_list[2][j] == i: # group is correct
			# go through the raw data set, pull group
			#print(main_list[2][j], "= main_list[2][j] == i =",i)
					myGeneGroup_list.append(main_list[1][j])
			#print("  myGeneGroup_list",myGeneGroup_list)
			self.groupGene_dic[int(i)] = myGeneGroup_list
		# print("\t+self.gene_dic : ",self.gene_dic)

# end of getGroups()

	def getNormFactor(self, ds_str, group_str):
		"""Method to get the raw expressions from a specified dataset (ds_str) and pull the expression values of the genes defined in the group (see NORM_FILE). Averages are then calculated of these expressions."""
		#print("  getNormFactor()")
		#print("\t: Dataset: ",ds_str, type(ds_str))
		#print("\t: group  : ",group_str, type(group_str))

		# what genes are we using for this group?
		#print("self.groupGene_dic : ",self.groupGene_dic, self.groupGene_dic.keys(), type(self.groupGene_dic))

		myGene_list = self.groupGene_dic[int(group_str)]
		#print("\t getNormFactor() myGene_list",myGene_list)

		# Collect the expression values for these genes in the specified dataset
		#printer(myGene_list)
		#print(self.raw_dic[ds_str])

		setGenes_list = self.raw_dic[ds_str][0] # contains the ensNumbers for this dataset
		genesExp_list = self.raw_dic[ds_str][1] # contains the expressions for this dataset
		loc_list = [] # contains the locations of the ensNums in setGenes_list. locations should be same in the expressions list
		exp_dic = {} # containt the exp (value) for the ensNum (key)

		for i in myGene_list:
			#print(" ++++ Now searching for : ", i, "in ds =", ds_str)
			loc_list.append(setGenes_list.index(i))
		# print(loc_list)

		for i in loc_list: # the index are locations
			#print("This i :",i)
			exp_dic[i] = genesExp_list[i]
			#printer(exp_dic)

    # get an average
		sum = 0
		for i in exp_dic:
			sum += float(exp_dic[i])
		avg_float = sum/len(exp_dic)
		#print("\t+ getNormFactor() avg_float: ",avg_float, " for group: ",group_str)
		return avg_float
# end of getNormFactor()

	def getNormMatrix(self,group_int):
		""" returns a matrix for which the variables have been normalized according to the genes listed in the NORM_FILE"""
		#self.raw_dic[ds_str][0] # contains the ensNumbers for this dataset
		#self.raw_dic[ds_str][1] # contains the expressions for this dataset
		#self.raw_dic[ds_str][2] # contains the counters for this dataset

		print("\n\tgetNormMatrix(): Creating normalizations for group:", group_int)
		# go through all datasets to make a list of normalized values. add this list to the raw_dic for the associated key
		headers_list = self.raw_dic.keys()
		norm_list = [] # contains the normalized values
		for i in headers_list: # the datasets, ex: "data/d578e27f-537c-4aaa-8903-6ffe68346276.FPKM.txt"
			norm_list = [] # reset the list
			#print("\t Key = ",i)# show the current dataset
			normFactor = self.getNormFactor(i, group_int)
			#print("\t+ NormFactor:",normFactor)
			if normFactor == 0: # normFactor not present?
				normFactor = 1

			l_dic = self.raw_dic[i] # the dictionary for the dataset. [0]: ensNums, [1]: expVal, [2]: rawlogs
			#ensNum_list = self.raw_dic[i][0]
			expVal_list = self.raw_dic[i][1] # [0] is ensNums, [1] is the raw expressions in self.raw_dic
			for exp in expVal_list: # need to add this new list to the raw dictionary for current key
				try: # if valueErrors exist
					norm_list.append(math.log((float(exp)/normFactor),math.exp(1)))
				except ValueError:
					norm_list.append(0)

			#print(" norm_list: ",norm_list,i)
			# add the list to self.raw_dic
			#print("BEFORE ADDED length of raw_dic:", len(self.raw_dic[i]),"group_int : ",group_int)
			self.raw_dic[i] = self.raw_dic[i] + [norm_list]
			#print("AFTER ADDED length of raw_dic:", len(self.raw_dic[i]),"group_int : ",group_int)

# note: adding another list to a key in a dictionary
#       h_dic = {}
#       h_dic[1] = [["one"],[1]]
#       h_dic[1] = h_dic[1] + [["un"]]

#end of getNormMatrix()


	def makeHeatmap1(self, group_int):
		"""A method to produce lists for Plotly heatmaps from the self.raw_dic dataset. """
# note:
#self.raw_dic[ds_str][0] # contains the ensNumbers for this dataset
#self.raw_dic[ds_str][1] # contains the raw expressions for this dataset
#self.raw_dic[ds_str][2] # contains the rawlogs of expressions for this dataset
#self.raw_dic[ds_str][3] # contains the AVGg1 normalizing
#self.raw_dic[ds_str][4] # contains the AVGg2 normalizing
#self.raw_dic[ds_str][5] # contains the AVGg3 normalizing
#self.raw_dic[ds_str][6] # contains the tubb normalizing
#self.raw_dic[ds_str][7] # contains the tuba1a normalizing


# note:
#g_dic = {}
#g_dic[1] = [[0],[1],[2],[3],[4]]
#g_dic[1][0]
#[0]

# note:
#first_list = []
#my_list = [1,2,3]
#first_list.append(my_list)
#my_list = [4,5,6]
#first_list.append(my_list)
#first_list
#[[1, 2, 3], [4, 5, 6]]


		#print("\t+ makeHeatmap1()")

		# go through all datasets to make a list for each type of data.
		headers_list = self.raw_dic.keys()
		plotThis_list = [] # the matrix of data for heatmaps. This takes the 'z' coordinate


		# build the lists of data. Each dataset is in own list. and has same format shown below
		# z = [[ds1_list], [ds2_list],...,[ds10_list]]

		my_list = []
		x_list = [] # contains the x-axis labels
		y_list = [] # containt the y-axis labels

		for h in headers_list:
			#print("\t Key = ",h, type(h), "for group :",group_int)# show the current dataset
			y_list.append(str(h[5:13])) # contains the shortened name of dataset

			# For each dataset key, pull the hth list of the dictionary's value.
			#            print( "BEFORE Adding: len(my_list) : ",len(my_list))
			my_list.append(self.raw_dic[h][group_int])
#            print( "AFTER Adding: len(my_list) : ",len(my_list))

		x_list = self.raw_dic[h][0]
		z_list = my_list
		#        myFname_str = "/tmp/" + str(group_int)+"_raw" + ".html"
		myFname_str = OUTPUT_DIR + str(group_int)+"_raw" + ".html"
		self.drawHeatmap(x_list, y_list, z_list, myFname_str)

    #end of makeHeatmap1()

# TODO: build a working list of lists to do some regression. Maybe, the giant matrix is not nexessary afterall because
# we have all the data in the self.matrix_dic.


	def drawHeatmap(self, x_list, y_list, z_list, myFname_str):
		""" Method to create the heatmaps from list inputs"""
# ref: https://plot.ly/python/heatmaps/

# debugging
#        printer(x_list)
#        print("  x_list size :",len(x_list)),
#        input("  The above is x_list: enter to continue")

#        printer(y_list)
#        print("  y_list size :",len(y_list)),
#        input("  The above is y_list: enter to continue")

#        printer(z_list)
#        print("  z_list size :",len(z_list))
#        input("  The above is z_list: enter to continue")

		trace = go.Heatmap(z = z_list, x = x_list, y = y_list)
		data=[trace]
		print("              +Saving heatmap file: ",myFname_str)
		plot(data, filename = myFname_str)

#end of drawHeatmap()


	def getRsquaredHeatmap(self, whichList_str):
		""" Across all datasets, make a list of all expression values for each gene. These lists will be used to get r-squared values from linear regression. Use values in self.raw_dic. Note: whichList_str is the column value for the key in self.raw_dic. """

		print("\n\tgetRsquaredHeatmap(): creating the R-squared values for the heatmaps: ",whichList_str )# which list in value of self.raw_dic.

		headers_list = self.raw_dic.keys()
		geneName_list = [] # contains the ensNums in the current dataset
		geneNamePos_dic = {} # ensNum (keys) pos(value); used to grab all gene expressions of same ensNum.
		# get all gene names:

		for i in headers_list: # the datasets, ex: "data/d578e27f-537c-4aaa-8903-6ffe68346276.FPKM.txt"
			geneName_list = self.raw_dic[i][0]
			for name in geneName_list:
				#print("\t Data set : ",i, " name : ",name) # lists the datasets as filenames.
				genePos_int = self.raw_dic[i][0].index(name)
				#print("\tpos: ",genePos_int)
				geneNamePos_dic[name] = genePos_int
			break # since all genename lists will be same, make a list from one iteration.
		#printer(geneNamePos_dic)

		exp_list = [] # contains the expressions for a particular gene across all datasets.
		geneExp_dic = {} # key: ensNum, value, a list of all expressions of this gene across datasets
		for ensNum in geneNamePos_dic: # p is ensNum, dic value is the position of expression in current dataset
			#print("p = ",p,", geneNamePos_dic[p] = ",geneNamePos_dic[p])
			exp_list = []
# need to gather a specific gene expression from each dataset
			for i in headers_list: # for each dataset...
				l_list = self.raw_dic[i][whichList_str] # make a list of these values.

#                print("l_list[geneNamePos_dic[p]:", l_list[geneNamePos_dic[p]])
#                exp_list.append(l_list[geneNamePos_dic[p]])
				exp_list.append(float(l_list[geneNamePos_dic[ensNum]]))
			geneExp_dic[ensNum] = exp_list # key: ensNum, value, a list of all expressions of this gene across datasets
#        print("geneExp_dic: ",geneExp_dic, len(geneExp_dic))

		#perform the r_squared calculation.
		counter = 0
		upperbound_int = len(geneExp_dic) * len(geneExp_dic)
		x_list = [] # x axis values for heatmap
		y_list = [] # y axis values for heatmap
		z_list = [] # r-sqrt values for heatmap
		miniz_list = [] # contains the row by row
		for x in geneExp_dic: # "i's are EnsNum values
			x_list.append(x)
			y_list.append(x)
			miniz_list = []
			for y in geneExp_dic:

				#print("geneExp_dic[x] x geneExp_dic[y] <-", geneExp_dic[x],geneExp_dic[y])
				res_float = self.getRSquaredScore(geneExp_dic[x],geneExp_dic[y])

# note: for x in range(100000):print("Progress {:2.1%}".format(x / 10), end="\r")

				print("{:5} of {:5}".format(counter, upperbound_int), end = '\r')
#                print("res_float : ",res_float)

				counter += 1
				miniz_list.append(res_float) # for this j value
			z_list.append(miniz_list)

#        myFname_str = "/tmp/" + str(whichList_str) + "_rsqu.html"
		myFname_str = OUTPUT_DIR + str(whichList_str) + "_rsqu.html"
		self.drawHeatmap(x_list, y_list, z_list, myFname_str)
		self.filterMatrix(x_list, y_list, z_list, whichList_str) #check the matrix in light of applied threshold constraints
		#end of getRsquaredHeatmap(self, whichList_str):

	def getRSquaredScore(self, d1_list, d2_list):
		""" get the R-squared value for a linear model between two lists"""
# reference:
# >>> from scipy import stats
# >>> help(stats)
		#print("getRSquaredScore()\n",d1_list,"and",d2_list)
		slope, intercept, r_value, p_value, std_err = stats.linregress(d1_list, d2_list)
		#print("slope: %f    intercept: %f" % (slope, intercept))
		res_float = r_value**2
		#print("R-squared: %f" % res_float)
		return res_float

		#end of geRSquareScore(d1_list, d2_list)

	def filterMatrix(self, x_list, y_list, z_list, whichList_str):
		""" Method to build a matrix according to min / max criteria for each value. The thresholds represent the upperbounds of values to keep."""
		print("\n\tfilterMatrix(): Creating heatmaps by threshold values.")
#        printer(self.thresh_list)
# debugging
#        printer(x_list)
#        print("  x_list size :",len(x_list)),
#        input("  The above is x_list: enter to continue")

#        printer(y_list)
#        print("  y_list size :",len(y_list)),
#        input("  The above is y_list: enter to continue")

#        printer(z_list)
        # print("  z_list size :",len(z_list))
        # print("  z_list[0] size :",len(z_list[0]))
        # input("  The above is z_list: enter to continue")

# overall idea: all values that are not less than the threshold will become zeros.

#note: y and z lists have to be same sizes.
		for tVal in self.thresh_list:
			#print("  Current value of tVal",tVal,type(tVal))
#            # clone the list.
			miniz_list = []
			zz_list = []
			for x in range(len(z_list)):
				my_list = z_list[x]
				for y in range(len(my_list)):
					if my_list[y] <= tVal:
						miniz_list.append(my_list[y])
					else:
						miniz_list.append(0)
				zz_list.append(miniz_list)
				miniz_list = []

#            myFname_str = "/tmp/" + str(whichList_str)+"_thresh_"+ str(tVal) +".html"
			myFname_str = OUTPUT_DIR + str(whichList_str)+"_thresh_"+ str(tVal) +".html"
			self.threshold_dic[myFname_str] = [x_list, y_list, z_list]
			self.drawHeatmap(x_list, y_list, zz_list, myFname_str)
#example:
#z_list = [['a',1, 2], ['b',3, 4], ['c',5, 6]]
#for i in range(len(z_list)):
#      print("i = ",i,z_list[i])
#      my_list = z_list[i]
#      for j in my_list:
#          print("j :",j)

    #end of filterMatrix(x_list, y_list, z_list, whichList_str)

	def compHeatmaps(self, name1_str, name2_str, m1_list, m2_list):
		""" make a comparison of all heatmaps stored in two lists of lists (m1_list and m2_list). We are compairng the values of the z_list (3rd list) to determine if both are non-zero or zero values."""

		#note:
		# s_list = [['xa','xb','xc'],['ya','yb','yc'],['za','zb','zc']]
		# sd_dic = {}
		# sd_dic["one"] = s_list
		# sd_dic
		# -> {'one': [['xa', 'xb', 'xc'], ['ya', 'yb', 'yc'], ['za', 'zb', 'zc']]}

		# sd_dic["one"][2]
		# -> ['za', 'zb', 'zc']

		# s_list = [['xxa','xxb','xxc'],['yya','yyb','yyc'],['zza','zzb','zzc']]
		# sd_dic["two"] = s_list
		# len(sd_dic)
		# -> 2
# input lists
# m1_list =  [['ENSG00000172137.17', 'ENSG00000167700.7', 'ENSG00000240423.1', 'ENSG00000060642.9'], ['ENSG00000172137.17', 'ENSG00000167700.7', 'ENSG00000240423.1', 'ENSG00000060642.9'],
# [[1.0, 0.2, 0.0, 4.0], [0.5, 6.0, 0.7, 0.8], [0.123, 0.345, 456.0, 567.0], [3.0, 0, 3.0, 0]]]
#
# m2_list =  [['ENSG00000172137.17', 'ENSG00000167700.7', 'ENSG00000240423.1', 'ENSG00000060642.9'], ['ENSG00000172137.17', 'ENSG00000167700.7', 'ENSG00000240423.1', 'ENSG00000060642.9'],
# [[4.0, 0.0, 0, 2.0], [2.0, 5.0, 0.0, 1.0], [0.0, 0.350, .750, .870], [0.340, 0.670, 0.120, 0.13]]]

		print("\n\t + compHeatmaps()")
		#print("m1_list = ",m1_list)
		#print("m2_list = ",m2_list)

		x_list = m1_list[0]
		y_list = m1_list[1]

		m1_list = m1_list[2]
		m2_list = m2_list[2]
		zrow_list = []
		zz_list = [] # contains the binary. 1 is both are non-zero, 0 one is non-zero
		#counter = 1
		for i in range(len(m1_list)):
			rows1_list = m1_list[i] # a row of the whole
			rows2_list = m2_list[i] # a row of the whole
			#print("rows1_list :",rows1_list)
			#print("rows2_list :",rows2_list)

			for k in range(len(rows1_list)):
				#print("k from rows1_list :",rows1_list[k])
				#print("k from rows2_list :",rows2_list[k])
				#print("Counter : ",counter,"\n")
				#counter += 1
				COMPARISON_THRESH = 0.7
				#                if rows1_list[k] != 0 and rows2_list[k] != 0:
				#                if rows1_list[k] < COMPARISON_THRESH  and rows2_list[k] < COMPARISON_THRESH:
				if rows1_list[k] >= COMPARISON_THRESH  and rows2_list[k] >= COMPARISON_THRESH:
					zrow_list.append(1)
				else:
					zrow_list.append(0)

			zz_list.append(zrow_list)
			zrow_list = []
		#print(zz_list)

#format the name of the html file output
		myFname_str = str(name1_str) + "_and_"+str(name2_str)
		myFname_str = myFname_str.replace(OUTPUT_DIR,"").replace(".html","")+"_comparison.html"
		myFname_str = OUTPUT_DIR + myFname_str

		self.drawHeatmap(x_list, y_list, zz_list, myFname_str)
		# end of compHeatmaps()


##############################################################
# end of Wrangler class
##############################################################




def heatmapDemo():
	#import plotly.plotly as py
	import chart_studio as py
	import plotly.graph_objs as go
	from plotly.offline import download_plotlyjs, init_notebook_mode, plot,iplot


# Add colour scales. see ref: https://plot.ly/python/colorscales/

	trace = go.Heatmap(z=[[1, 20, 30, 50, 1], [20, 1, 60, 80, 30], [30, 60, 1, -10, 20]],
	               x=['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
	               y=['Morning', 'Afternoon', 'Evening'])
	data=[trace]
	# check the output dir
	tmp_dir = checkDataDir(OUTPUT_DIR)

#    plot(data, filename = '/tmp/labelled-heatmap.html')
	myFname_str = OUTPUT_DIR+'labelled-heatmap.html'
	print("              +Saving heatmap file: ",myFname_str)
	plot(data, filename = myFname_str)

	#end of heatmapDemo()





def RSquaredScoreDemo():
	""" get the R-squared value for a linear model between"""

	print("\n\tRSquaredScoreDemo(): \n\tIf this calculation is performed, then your libraries are correctly installed.")
	xx = [153.0, 549953.824845, 11.3589638951, 15.6436361402, 16.4413743419, 16.1267035586, 22.5096888809, 2123.0, 151243.109382, 16.0546348885]
	yy = [0.0, 8335741.97334, 194.355883931, 458.169453256, 266.523165476, 762.964596253, 339.411793374, 0.0, 5098152.58481, 213.000491758]
	print("\n\tData: ")
	print("\t x= ",xx)
	print("\n\t y=",yy)
	print("\n\t r2_score:", r2_score(yy, xx))
	slope, intercept, r_value, p_value, std_err = stats.linregress(yy, xx)
	print("\t Slope: %f    intercept: %f" % (slope, intercept))
	print("\t R-squared: %f" % r_value**2)

#end of RSquareScore()

def saveFile(in_str):
# not currently used in this version...
	"""Save the markdown string as a text file"""
	if len(in_str) > 0:
		try:
			tmp_dir = checkDataDir(OUTPUT_DIR)
			fname = "out.md"
			filename = OUTPUT_DIR + fname
			f = open(filename, "w")
			f.write(in_str)
			f.close()
			print(" + Saved md file of results in <",filename,"> ")
		except IOError:
			print("  Problem saving file... incorrect permissions?!")
	# end of saveFile()

def printer(inThing):
	""" prints things cleanly"""
	if "list" in str(type(inThing)):
		for i in range( len(inThing) ):
			print("\t",i,":", inThing[i])
	if "dict" in str(type(inThing)):
		counter = 0
		for i in inThing:
			print("\t",counter," |  ",i,":",inThing[i])
			counter += 1
#end of printer()


def checkDataDir(dir_str):
#function to determine whether a data output directory exists.
#if the directory doesnt exist, then it is created

# not currently used in this version...

	try:
		os.makedirs(dir_str)
		#print("  PROBLEM: output_dir doesn't exist")
		print("  * Creating :",dir_str)
		return 1

	except OSError:
		return 0
#end of checkDataDir()




def begin(task_str):
	"""Driver function of program"""
	print("\n\t Welcome to geneExPy!\n\t A heatmap generator of expression data.")
	#print(" Task :",task_str)
	# get current directory
	#cwd = os.getcwd()
	if task_str == "heatmap":
		print("\n\t+ Running heatmap demo...")
		heatmapDemo()
		exit()

	if task_str == "rsqu":
		print("\n\t+ Running r-square geneeration demo...")
		RSquaredScoreDemo()
		exit()


# check the output dir
	tmp_dir = checkDataDir(OUTPUT_DIR)

# check the input dir
	tmp_dir = checkDataDir(INPUT_DIR)

# define the class and determine which files are for data and which are to be used to normalize
	print("\n\t+ The input data  : ",INPUT_DIR)
	k = Wrangler()
	fileListing_list = k.getFileListing(INPUT_DIR)  # get a listing of the files out there in the dataInput dir

# show which files we are loading?
	print("\n\t+ Building matrix from these files :")
	printer(fileListing_list)

#load the data; create a data structure to contain the data called raw_dic
	#print("\n\t+ begin() Opening the data files ...")
	k.getParams() # open the picker file and the thresholds file.

	k.getRawMatrix() # the list is the task and the file. here no parameter is necessary

#load the genes to be used to make normalizing factors
	#print("\n\t+ begin() Opening the normNames file ...")
	k.getNormNamesMatrix(NORM_FILE) # find names of genes to average, second part of list is file to open.
	#print("\t+",k.normNames_dic.keys())


# determine which averaging groups there are. This information is in the NORM_FILE.
	#print("\n\t+ begin() getGroup, determining groups")
	k.getGroups()
	#print("\n\t+ group set, k.group_set  : ",  k.group_set)
	print("\n\t+ Normalizing by genes of groups: ")
	printer(k.groupGene_dic) # keys: group, values: list of genes

######################################################
# make normalizations
######################################################

# These below methods do not handle the rsquared values.
# note the value as parameter is the group number defined in the NORM_FILE
# and determines the normalizing factor for the dataset.
# the self.raw_dic has the three initial lists for each key: ensNums, raw, rawLogs.

# The group argument has been defined in the NORM_FILE. The argument is the value for the group from the NORM_FILE (group).
	print("\tNumber of iterations for k.groupGene_dic :",len(k.groupGene_dic))
	for i in range(1,len(k.groupGene_dic)+1,1):
		#print("k.groupGene_dic value :",i)
		k.getNormMatrix(i)



# manually
#    k.getNormMatrix(1) # group 1 from NORM_FILE,
#    k.getNormMatrix(2) # group 2,
#    k.getNormMatrix(3) # group 3,

#    k.getNormMatrix(4) # group 4,
#    k.getNormMatrix(5) # group 5,
#    k.getNormMatrix(6) # group 6,




######################################################
# make some heatmaps of raw and normalized data
######################################################


# note: for this function, the first and second "groups" are the raw and the rawlogs lists. The groups defined in the NORM_FILE begin at value 3. Use norm_list group_number + 2 to make a heatmap of the group from NORM_FILE


	print("\tNumber of iterations for k.makeHeatmap1 :",len(k.groupGene_dic)+2)
	# Note for below: add two for raw and rawlogs, then 1 to offset the range that starts at 1.

	for i in range(1,len(k.groupGene_dic)+2+1,1):
		#print("k.makeHeatmap1(i) value :",i)
		k.makeHeatmap1(i)


# manually
#    k.makeHeatmap1(1) # make a heatmap of raw data (not normalized)
#    k.makeHeatmap1(2) # rawLogs

#    k.makeHeatmap1(3) # beginning of the groups defined in the NORM_FILE
#    k.makeHeatmap1(4) # group defined in the NORM_FILE

#    k.makeHeatmap1(5) # group defined in the NORM_FILE
#    k.makeHeatmap1(6) # group defined in the NORM_FILE
#    k.makeHeatmap1(7) # group defined in the NORM_FILE






######################################################
# R-squared analysis
# Note: whichList_str is the column value for the key in self.raw_dic.
######################################################

	#how many columns do we have for making heatmaps?
	for i in k.raw_dic:
#        print(i,len(k.raw_dic[i]))
		columns_int = len(k.raw_dic[i])
		break



	print("\tNumber of iterations for k.getRsquaredHeatmap :",columns_int)
	for i in range(1,columns_int,1):
		k.getRsquaredHeatmap(i)



# manually
#    k.getRsquaredHeatmap(1) # first column or raw
#    k.getRsquaredHeatmap(2) # second column or rawlogs

#    k.getRsquaredHeatmap(3) # group
#    k.getRsquaredHeatmap(4) # group
#    k.getRsquaredHeatmap(5) # group

	print("\n\tComparison of threshold heatmaps: binary outputs.")
	#print("\tKeys:",k.rSquMat_dic.keys())

	#key_list = [i for i in k.rSquMat_dic.keys()]
	key_list = [i for i in k.threshold_dic.keys()]

#make permutation of the list
#import itertools
#list(itertools.permutations([1, 2, 3]))


	print("key_list: ",key_list)
	donePairs_dic = {} # store which pairs have already been processed
	p_list = []
	for i in range(len(key_list)):
		for j in range(len(key_list)):
			p_list = []
			#print("begin() key i :",i, key_list[i])
			#print("begin() key j : ",j, key_list[j])
			p_list.append(i)
			p_list.append(j)
			p_list = sorted(p_list)
			p_str = str(sorted(p_list))
			if p_str not in donePairs_dic:
				donePairs_dic[p_str] = 1
				#print("Comparison of this pair: ",p_list)
				try:
					k.compHeatmaps(key_list[i], key_list[j], k.threshold_dic[key_list[p_list[0]]], k.threshold_dic[key_list[p_list[1]]])
				except TypeError:
					pass
			else:
				pass


######################################################


######################################################
#end of begin()



if __name__ == '__main__':

	if len(sys.argv) == 2: # one parameter
		begin(sys.argv[1])#,sys.argv[2])#,sys.argv[3], sys.argv[4]),sys.argv[5])
	else:
		help()
		sys.exit(0)
