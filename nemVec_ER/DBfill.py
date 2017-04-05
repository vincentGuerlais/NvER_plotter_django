#DBfill allows to fill the NvERTx Database ;
#how to use it : 
#	1. run the shell via python (in the /dER_plotter/nemVec_ER dir) : python manage.py shell
#	2. import this script in python : import DBfill
#	3. execute with the name of the file and the name of the table : DBfill.DBFill('fileName','table')
#	the tables are : Regen_cpm ; Fasta ; ; ; ; you can chek them in the admin panel
#	4. It may take a while depending on the size of you file 'done filling the DB' should print when done
#For now, it works for :
#SQLite DB ; cpm (regen) data
#Doesn't work if datas are missing !!!

##################
### Var
##################

import os
from ER_plotter.models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Standard_Error

##################
### Functions
##################

#Split the data on the tabs
def dataSplit(dataLine) :
	return dataLine[:-1].split('\t')

#Create the query to fill the database
def queryCreate(splitLine, DBTableName):
	if DBTableName == 'Regen_cpm' :
		query = "%s(nvertx_id = %s, regen_anc_UC= %3.15f, regen_anc_0HPA= %3.15f, regen_anc_2HPA= %3.15f, regen_anc_4HPA= %3.15f, regen_anc_8HPA= %3.15f, regen_anc_12HPA= %3.15f, regen_anc_16HPA= %3.15f, regen_anc_20HPA= %3.15f, regen_anc_24HPA= %3.15f, regen_anc_36HPA= %3.15f, regen_anc_48HPA= %3.15f, regen_anc_60HPA= %3.15f, regen_anc_72HPA= %3.15f, regen_anc_96HPA= %3.15f, regen_anc_120HPA= %3.15f, regen_anc_144HPA= %3.15f).save()" % (DBTableName,splitLine[0],float(splitLine[1]),float(splitLine[2]),float(splitLine[3]),float(splitLine[4]),float(splitLine[5]),float(splitLine[6]),float(splitLine[7]),float(splitLine[8]),float(splitLine[9]),float(splitLine[10]),float(splitLine[11]),float(splitLine[12]),float(splitLine[13]),float(splitLine[14]),float(splitLine[15]),float(splitLine[16]))
	
	elif DBTableName == 'Fasta' :
		query = "%s(nvertx_id = %s, fasta_sequence= %s).save()" % (DBTableName,splitLine[0],splitLine[1])
	
	elif DBTableName == 'Embryo_cpm' :
		query = "Embryo_cpm(nvertx_id = " + splitLine[0]
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
	
	elif DBTableName == 'Annotation' :
		query = "Annotation(nvertx_id = " + splitLine[0]
		if splitLine[2] != 'NA' :
			query += ", nve_hit = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", nve_eval = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", mfuzz_clust = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", mfuzz_score = " + splitLine[5]
		if splitLine[6] != 'NA' and splitLine[6] != 'No_Uniprotmatch' :
			query += ", uniprot_id = " + splitLine[6]
		if splitLine[7] != 'NA' and splitLine[7] != 'No_Description' :
			query += ", uniprot_description = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", top_nr_hit_eval = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", other_nr_hits = " + splitLine[9]
		query += ").save()"
	elif DBTableName == 'Standard_Error' :
		query = "%s(nvertx_id = %s,regen_se_UC = %3.15f,regen_se_0HPA = %3.15f,regen_se_2HPA = %3.15f,regen_se_4HPA = %3.15f,regen_se_8HPA = %3.15f,regen_se_12HPA = %3.15f,regen_se_16HPA = %3.15f,regen_se_20HPA = %3.15f,regen_se_24HPA = %3.15f,regen_se_36HPA = %3.15f,regen_se_48HPA = %3.15f,regen_se_60HPA = %3.15f,regen_se_72HPA = %3.15f,regen_se_96HPA = %3.15f,regen_se_120HPA = %3.15f,regen_se_144HPA = %3.15f,regen_log_se_UC = %3.15f,regen_log_se_0HPA = %3.15f,regen_log_se_2HPA = %3.15f,regen_log_se_4HPA = %3.15f,regen_log_se_8HPA = %3.15f,regen_log_se_12HPA = %3.15f,regen_log_se_16HPA = %3.15f,regen_log_se_20HPA = %3.15f,regen_log_se_24HPA = %3.15f,regen_log_se_36HPA = %3.15f,regen_log_se_48HPA = %3.15f,regen_log_se_60HPA = %3.15f,regen_log_se_72HPA = %3.15f,regen_log_se_96HPA = %3.15f,regen_log_se_120HPA = %3.15f,regen_log_se_144HPA = %3.15f).save()" % (DBTableName,splitLine[1],float(splitLine[2]),float(splitLine[3]),float(splitLine[4]),float(splitLine[5]),float(splitLine[6]),float(splitLine[7]),float(splitLine[8]),float(splitLine[9]),float(splitLine[10]),float(splitLine[11]),float(splitLine[12]),float(splitLine[13]),float(splitLine[14]),float(splitLine[15]),float(splitLine[16]),float(splitLine[17]),float(splitLine[18]),float(splitLine[19]),float(splitLine[20]),float(splitLine[21]),float(splitLine[22]),float(splitLine[23]),float(splitLine[24]),float(splitLine[25]),float(splitLine[26]),float(splitLine[27]),float(splitLine[28]),float(splitLine[29]),float(splitLine[30]),float(splitLine[31]),float(splitLine[32]),float(splitLine[33]))
	else :
		print 'error in DBTableName'
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
	readDataFile.readline()
	for lines in readDataFile :
		splitLine = dataSplit(lines)
		DBQuery.append(queryCreate(splitLine,table))
	readDataFile.close()
	#save the modifications
	save(DBQuery)
	return 'done filling the DB'





