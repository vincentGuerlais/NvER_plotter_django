#DBfill allows to fill the NvERTx Database ;
#Use it on files with columns separated with tabulations (export table from R for example)
#how to use it : 
#	1. run the shell via python (in the /dER_plotter/nemVec_ER dir) : python manage.py shell
#	2. import this script in python : import DBfill
#	3. execute with the name of the file and the name of the table : DBfill.DBFill('fileName','table')
#	the tables are : Regen_cpm ; Fasta ; Embryo_cpm ; Annotation ; Regen_SE ; Embryo_SE ; you can check for their name in the admin panel
#	4. It may take a while depending on the size of your file 'done filling the DB' should print when done
#The tables Annotation and Embryo_cpm can handle missing data. they should be replaced by 'NA'

##################
### Var
##################

import os
from ER_plotter.models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Regen_SE, Embryo_SE, Regen_log_SE

##################
### Functions
##################

#Split the data on the tabs
def dataSplit(dataLine) :
	return dataLine[:-1].split('\t')

#Create the query to fill the database
def queryCreate(splitLine, DBTableName):
	if DBTableName == 'Regen_cpm' :
		query = "Regen_cpm(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", regen_anc_UC = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", regen_anc_0HPA = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", regen_anc_2HPA = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", regen_anc_4HPA = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_anc_8HPA = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_anc_12HPA = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_anc_16HPA = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_anc_20HPA = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", regen_anc_24HPA = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_anc_36HPA = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_anc_48HPA = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_anc_60HPA = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_anc_72HPA = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", regen_anc_96HPA = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_anc_120HPA = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_anc_144HPA = " + splitLine[16]
		query += ").save()"

	elif DBTableName == 'Fasta' :
		query = "%s(nvertx_id = '%s', fasta_sequence= '%s').save()" % (DBTableName,splitLine[0],splitLine[1])
		#if splitLine[0] == '>' :

	elif DBTableName == 'Embryo_cpm' :
		query = "Embryo_cpm(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", warner_anc_24HPF = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", warner_anc_48HPF = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", warner_anc_72HPF = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", warner_anc_96HPF = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", warner_anc_120HPF = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", warner_anc_144HPF = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", warner_anc_168HPF = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", warner_anc_192HPF = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", fischer_anc_0HPF = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", fischer_anc_1HPF = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", fischer_anc_2HPF = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", fischer_anc_3HPF = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", fischer_anc_4HPF = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", fischer_anc_5HPF = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", fischer_anc_6HPF = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", fischer_anc_7HPF = " + splitLine[16]
		if splitLine[17] != 'NA' :
			query += ", fischer_anc_8HPF = " + splitLine[17]
		if splitLine[18] != 'NA' :
			query += ", fischer_anc_9HPF = " + splitLine[18]
		if splitLine[19] != 'NA' :
			query += ", fischer_anc_10HPF = " + splitLine[19]
		if splitLine[20] != 'NA' :
			query += ", fischer_anc_11HPF = " + splitLine[20]
		if splitLine[21] != 'NA' :
			query += ", fischer_anc_12HPF = " + splitLine[21]
		if splitLine[22] != 'NA' :
			query += ", fischer_anc_13HPF = " + splitLine[22]
		if splitLine[23] != 'NA' :
			query += ", fischer_anc_14HPF = " + splitLine[23]
		if splitLine[24] != 'NA' :
			query += ", fischer_anc_15HPF = " + splitLine[24]
		if splitLine[25] != 'NA' :
			query += ", fischer_anc_16HPF = " + splitLine[25]
		if splitLine[26] != 'NA' :
			query += ", fischer_anc_17HPF = " + splitLine[26]
		if splitLine[27] != 'NA' :
			query += ", fischer_anc_18HPF = " + splitLine[27]
		if splitLine[28] != 'NA' :
			query += ", fischer_anc_19HPF = " + splitLine[28]
		if splitLine[29] != 'NA' :
			query += ", helm_anc_2HPF = " + splitLine[29]
		if splitLine[30] != 'NA' :
			query += ", helm_anc_7HPF = " + splitLine[30]
		if splitLine[31] != 'NA' :
			query += ", helm_anc_12HPF = " + splitLine[31]
		if splitLine[32] != 'NA' :
			query += ", helm_anc_24HPF = " + splitLine[32]
		if splitLine[33] != 'NA' :
			query += ", helm_anc_120HPF = " + splitLine[33]
		if splitLine[34] != 'NA' :
			query += ", helm_anc_240HPF = " + splitLine[34]
		if splitLine[35] != 'NA' :
			query += ", mean_2HPF = " + splitLine[35]
		if splitLine[36] != 'NA' :
			query += ", mean_7HPF = " + splitLine[36]
		if splitLine[37] != 'NA' :
			query += ", mean_12HPF = " + splitLine[37]
		if splitLine[38] != 'NA' :
			query += ", mean_24HPF = " + splitLine[38]
		if splitLine[39] != 'NA' :
			query += ", mean_120HPF = " + splitLine[39]
		query += ").save()"

	elif DBTableName == "Annotation" :
		query = 'Annotation(nvertx_id = "' + splitLine[0] +'"'
		query += ', Nemve1_tophit = "' + splitLine[2] +'"'
		if splitLine[3] != "NA" :
			query += ', Nemve1_e_val = ' + splitLine[3]
		else :
			query += ', Nemve1_e_val = 99'
		query += ', Mfuzz_R_Clust = "' + splitLine[4] +'"'
		if splitLine[5] != "NA" :
			query += ', Mfuzz_R_Score = ' + splitLine[5]
		else :
			query += ', Mfuzz_R_Score = -1'
		query += ', Mfuzz_E_Clust = "' + splitLine[6] +'"'
		if splitLine[7] != "NA" :
			query += ', Mfuzz_E_Score = ' + splitLine[7]
		else :
			query += ', Mfuzz_E_Score = -1'
		query += ', Uniprot_ID = "' + splitLine[8] +'"'
		query += ', Uniprot_Description = "' + splitLine[9] +'"'
		query += ', Top_nr_hit_eval = "' + splitLine[10] +'"'
		query += ', Other_nr_hits = "' + splitLine[11] +'"'
		query += ').save()'
	
	elif DBTableName == 'Regen_SE' :
		query = "Regen_SE(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", regen_se_UC = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", regen_se_0HPA = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", regen_se_2HPA = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", regen_se_4HPA = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_se_8HPA = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_se_12HPA = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_se_16HPA = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_se_20HPA = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", regen_se_24HPA = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_se_36HPA = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_se_48HPA = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_se_60HPA = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_se_72HPA = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", regen_se_96HPA = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_se_120HPA = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_se_144HPA = " + splitLine[16]
		query += ").save()"
	
	elif DBTableName == 'Regen_log_SE' :
		query = "Regen_log_SE(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", regen_log_se_UC = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", regen_log_se_0HPA = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", regen_log_se_2HPA = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", regen_log_se_4HPA = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_log_se_8HPA = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_log_se_12HPA = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_log_se_16HPA = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_log_se_20HPA = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", regen_log_se_24HPA = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_log_se_36HPA = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_log_se_48HPA = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_log_se_60HPA = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_log_se_72HPA = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", regen_log_se_96HPA = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_log_se_120HPA = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_log_se_144HPA = " + splitLine[16]
		query += ").save()"
	
	elif DBTableName == 'Embryo_SE' :
		query = "Embryo_SE(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", warner_se_24HPF = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", warner_se_48HPF = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", warner_se_72HPF = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", warner_se_96HPF = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", warner_se_120HPF = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", warner_se_144HPF = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", warner_se_168HPF = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", warner_se_192HPF = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", fischer_se_0HPF = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", fischer_se_1HPF = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", fischer_se_2HPF = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", fischer_se_3HPF = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", fischer_se_4HPF = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", fischer_se_5HPF = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", fischer_se_6HPF = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", fischer_se_7HPF = " + splitLine[16]
		if splitLine[17] != 'NA' :
			query += ", fischer_se_8HPF = " + splitLine[17]
		if splitLine[18] != 'NA' :
			query += ", fischer_se_9HPF = " + splitLine[18]
		if splitLine[19] != 'NA' :
			query += ", fischer_se_10HPF = " + splitLine[19]
		if splitLine[20] != 'NA' :
			query += ", fischer_se_11HPF = " + splitLine[20]
		if splitLine[21] != 'NA' :
			query += ", fischer_se_12HPF = " + splitLine[21]
		if splitLine[22] != 'NA' :
			query += ", fischer_se_13HPF = " + splitLine[22]
		if splitLine[23] != 'NA' :
			query += ", fischer_se_14HPF = " + splitLine[23]
		if splitLine[24] != 'NA' :
			query += ", fischer_se_15HPF = " + splitLine[24]
		if splitLine[25] != 'NA' :
			query += ", fischer_se_16HPF = " + splitLine[25]
		if splitLine[26] != 'NA' :
			query += ", fischer_se_17HPF = " + splitLine[26]
		if splitLine[27] != 'NA' :
			query += ", fischer_se_18HPF = " + splitLine[27]
		if splitLine[28] != 'NA' :
			query += ", fischer_se_19HPF = " + splitLine[28]
		if splitLine[29] != 'NA' :
			query += ", helm_se_2HPF = " + splitLine[29]
		if splitLine[30] != 'NA' :
			query += ", helm_se_7HPF = " + splitLine[30]
		if splitLine[31] != 'NA' :
			query += ", helm_se_12HPF = " + splitLine[31]
		if splitLine[32] != 'NA' :
			query += ", helm_se_24HPF = " + splitLine[32]
		if splitLine[33] != 'NA' :
			query += ", helm_se_120HPF = " + splitLine[33]
		if splitLine[34] != 'NA' :
			query += ", helm_se_240HPF = " + splitLine[34]
		query += ").save()"

	else :
		print "error in DBTableName. Valid inputs are : Regen_cpm ; Fasta ; Embryo_cpm ; Annotation ; Regen_SE ; Embryo_SE"
	return query

#execute the code
def save(DBQueryList) :
	for queries in DBQueryList :
		exec(queries)

##################
### Main
##################

#read the file and process :
def DBFill(fileName,table):
	#file to handle ; temp : to be replaced by command-line args
	dataFile = fileName
	DBQuery = []
	readDataFile = open(dataFile,'r')
	readDataFile.readline() #this only applies if there is a header to skip
	if table == 'Fasta' :
		splitLine = [readDataFile.readline()[1:-1],'']
		for line in readDataFile :
			if line[0] == '>' :
				DBQuery.append(queryCreate(splitLine,table))
				splitLine = [line[1:-1]]
				splitLine.append('')
			else :
				splitLine[1] += line.strip()
		DBQuery.append(queryCreate(splitLine,table))
	else :
		for lines in readDataFile :
			splitLine = dataSplit(lines)
			DBQuery.append(queryCreate(splitLine,table))
	readDataFile.close()
	#save the modifications
	save(DBQuery)
	return 'done filling the DB'





