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
		query = query = "%s(nvertx_id= %s, warner_anc_24HPF= %3.15f, warner_anc_48HPF= %3.15f, warner_anc_72HPF= %3.15f, warner_anc_96HPF= %3.15f, warner_anc_120HPF= %3.15f, warner_anc_144HPF= %3.15f, warner_anc_168HPF= %3.15f, warner_anc_192HPF= %3.15f, fischer_anc_0HPF= %3.15f, fischer_anc_1HPF= %3.15f, fischer_anc_2HPF= %3.15f, fischer_anc_3HPF= %3.15f, fischer_anc_4HPF= %3.15f, fischer_anc_5HPF= %3.15f, fischer_anc_6HPF= %3.15f, fischer_anc_7HPF= %3.15f, fischer_anc_8HPF= %3.15f, fischer_anc_9HPF= %3.15f, fischer_anc_10HPF= %3.15f, fischer_anc_11HPF= %3.15f, fischer_anc_12HPF= %3.15f, fischer_anc_13HPF= %3.15f, fischer_anc_14HPF= %3.15f, fischer_anc_15HPF= %3.15f, fischer_anc_16HPF= %3.15f, fischer_anc_17HPF= %3.15f, fischer_anc_18HPF= %3.15f, fischer_anc_19HPF= %3.15f, helm_anc_2HPF= %3.15f, helm_anc_7HPF= %3.15f, helm_anc_12HPF= %3.15f, helm_anc_24HPF= %3.15f, helm_anc_120HPF= %3.15f, helm_anc_240HPF= %3.15f, mean_2HPF= %3.15f, mean_7HPF= %3.15f, mean_12HPF= %3.15f, mean_24HPF= %3.15f, mean_120HPF= %3.15f).save()" % (DBTableName,splitLine[0],float(splitLine[1]),float(splitLine[2]),float(splitLine[3]),float(splitLine[4]),float(splitLine[5]),float(splitLine[6]),float(splitLine[7]),float(splitLine[8]),float(splitLine[9]),float(splitLine[10]),float(splitLine[11]),float(splitLine[12]),float(splitLine[13]),float(splitLine[14]),float(splitLine[15]),float(splitLine[16]),float(splitLine[17]),float(splitLine[18]),float(splitLine[19]),float(splitLine[20]),float(splitLine[21]),float(splitLine[22]),float(splitLine[23]),float(splitLine[24]),float(splitLine[25]),float(splitLine[26]),float(splitLine[27]),float(splitLine[28]),float(splitLine[29]),float(splitLine[30]),float(splitLine[31]),float(splitLine[32]),float(splitLine[33]),float(splitLine[34]),float(splitLine[35]),float(splitLine[36]),float(splitLine[37]),float(splitLine[38]),float(splitLine[39]))
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





