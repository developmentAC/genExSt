#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# # Oliver ######## . # #       # .  #
# # Bonham ######## . # #    # .  #
# # Carter ######## . # # # .  #
# ################# . # # ############
# ################# . # # # .  #
# ################# . # #    # .  #
# ################# . # #       # .  #
#
#
#
#
#

import streamlit as st
import genExSt_web_helperCode as hc


# DATE = "1 August 2020"
# VERSION = "1_i"
# AUTHOR = "Oliver Bonham-Carter"
# AUTHORMAIL = "obonhamcarter@allegheny.edu"


def begin():

	"""The driver function of the thing. """

	st.text(hc.banner0_str)
	st.sidebar.title("GenExSt; Gene Expression Analysis")

	st.sidebar.subheader(f"Date: {hc.DATE}, Ver: {hc.VERSION}")
	st.sidebar.text("\U0001F415 \U0001F631 \U0001f5ff \U0001F608 \U0001f600 ")
 	# Create a text element and let the reader know the data is loading.
	#

	# # Get directory to load the data files.
	# dataDir_str = hc.grabFile()
	# k = hc.Wrangler()
	#
	# fileListing_list = k.getFileListing(dataDir_str)  # get a listing of the files out there in the dataInput dir
	# # st.write(fileListing_list)


	# try:
	# 	data = hc.load_big_data(myFile_str)
	# 	# create a dictionary having headers as keys and values as lists of column data.
	# 	data_dic = hc.createMasterDataDic(data)
	# except:
	# 	st.sidebar.error("No data entered...")



# menu system
	doThis_sb = st.sidebar.selectbox(
		"What are we doing with this data?",
		[
			"ReadMe",
			# "Show Data Files",
			"Demo Heatmap",
			"Demo RSquared",
			"Gene Expression Analysis"
			],
	)
	if doThis_sb == "ReadMe":
		with open("README.md") as readme_file:
			st.markdown(readme_file.read())


	# if doThis_sb == "Show Data Files":
	# 	try:
	# 		hc.showData(fileListing_list)
	# 	except Exception:
	# 		st.write("Please wait")
	#

	if doThis_sb == "Gene Expression Analysis":
		st.write("Gene Expression")
		hc.geneExprSetup()


	if doThis_sb == "Demo Heatmap":
		hc.heatmapDemo()


	if doThis_sb == "Demo RSquared":
		hc.RSquaredScoreDemo()



	hc.writer("\U0001F415 WOO WOO!! \U0001F415")
	hc.writer(" ok :","computer")
		# end of begin()
#
#
# 	# progress_bar = st.progress(0)
# 	#
# 	# for i in range(100):
# 	# 	# Update progress bar.
# 	# 	progress_bar.progress(i)
#



begin()
