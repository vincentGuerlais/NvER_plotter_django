from __future__ import unicode_literals

from django.db import models

class Fasta(models.Model) :
	nvertx_id = models.CharField(max_length=20, primary_key=True)
	fasta_sequence = models.CharField(max_length=20000,null=True)
	def __str__(self) :
		return self.nvertx_id

class Regen_cpm(models.Model) :
	nvertx_id = models.CharField(max_length=20,primary_key=True)
	regen_anc_UC = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_0HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_2HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_4HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_8HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_12HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_16HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_20HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_24HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_36HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_48HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_60HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_72HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_96HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_120HPA = models.FloatField(max_length=17,null=True,blank=True)
	regen_anc_144HPA = models.FloatField(max_length=17,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id


class Embryo_cpm(models.Model) :
	nvertx_id = models.CharField(max_length=20,primary_key=True)
	warner_anc_24HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_48HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_72HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_96HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_120HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_144HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_168HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_anc_192HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_0HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_1HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_2HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_3HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_4HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_5HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_6HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_7HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_8HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_9HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_10HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_11HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_12HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_13HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_14HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_15HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_16HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_17HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_18HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_anc_19HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_anc_2HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_anc_7HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_anc_12HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_anc_24HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_anc_120HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_anc_240HPF = models.FloatField(max_length=10,null=True,blank=True)
	mean_2HPF = models.FloatField(max_length=10,null=True,blank=True)
	mean_7HPF = models.FloatField(max_length=10,null=True,blank=True)
	mean_12HPF = models.FloatField(max_length=10,null=True,blank=True)
	mean_24HPF = models.FloatField(max_length=10,null=True,blank=True)
	mean_120HPF = models.FloatField(max_length=10,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id
		
class Embryo_SE(models.Model) :
	nvertx_id = models.CharField(max_length=20,primary_key=True)
	warner_se_24HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_48HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_72HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_96HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_120HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_144HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_168HPF = models.FloatField(max_length=10,null=True,blank=True)
	warner_se_192HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_0HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_1HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_2HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_3HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_4HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_5HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_6HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_7HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_8HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_9HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_10HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_11HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_12HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_13HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_14HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_15HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_16HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_17HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_18HPF = models.FloatField(max_length=10,null=True,blank=True)
	fischer_se_19HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_se_2HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_se_7HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_se_12HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_se_24HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_se_120HPF = models.FloatField(max_length=10,null=True,blank=True)
	helm_se_240HPF = models.FloatField(max_length=10,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id

class Annotation(models.Model) :
	nvertx_id = models.CharField(max_length=17,primary_key=True)
	Nemve1_tophit = models.CharField(max_length=60,null=True,blank=True)
	Nemve1_e_val = models.FloatField(max_length=150,null=True,blank=True)
	Mfuzz_R_Clust = models.CharField(max_length=4,null=True,blank=True)
	Mfuzz_R_Score = models.FloatField(max_length=10,null=True,blank=True)
	Mfuzz_E_Clust = models.CharField(max_length=4,null=True,blank=True)
	Mfuzz_E_Score = models.FloatField(max_length=10,null=True,blank=True)
	Uniprot_ID = models.CharField(max_length=10,null=True,blank=True)
	Uniprot_Description = models.CharField(max_length=100,null=True,blank=True)
	Top_nr_hit_eval = models.CharField(max_length=1000,null=True,blank=True)
	Other_nr_hits = models.CharField(max_length=10000,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id		

class Regen_SE(models.Model) :
	nvertx_id = models.CharField(max_length=20,primary_key=True)
	regen_se_UC = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_0HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_2HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_4HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_8HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_12HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_16HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_20HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_24HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_36HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_48HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_60HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_72HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_96HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_120HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_se_144HPA = models.FloatField(max_length=10,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id
		
class Regen_log_SE(models.Model) :
	nvertx_id = models.CharField(max_length=20,primary_key=True)
	regen_log_se_UC = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_0HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_2HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_4HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_8HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_12HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_16HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_20HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_24HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_36HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_48HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_60HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_72HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_96HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_120HPA = models.FloatField(max_length=10,null=True,blank=True)
	regen_log_se_144HPA = models.FloatField(max_length=10,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id

class Mfuzz(models.Model) :
	#name = models.CharField(max_length=20,primary_key=True)
	mfuzz_cluster_nb = models.CharField(max_length=4,primary_key=True)
	cluster_image = models.CharField(max_length=50)
	bp_plot_image = models.CharField(max_length=50)
	#def __str__(self) :
	#	return self.name















