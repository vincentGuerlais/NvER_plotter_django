from __future__ import unicode_literals

from django.db import models

class Fasta(models.Model) :
	nvertx_id = models.CharField(max_length=20, primary_key=True)
	fasta_sequence = models.CharField(max_length=10000,null=True)
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

class Annotation(models.Model) :
	nvertx_id = models.CharField(max_length=20,primary_key=True)
	nve_hit = models.CharField(max_length=30,null=True,blank=True)
	nve_eval = models.FloatField(max_length=10,null=True,blank=True)
	mfuzz_clust = models.IntegerField(null=True,blank=True)
	mfuzz_score = models.FloatField(max_length=10,null=True,blank=True)
	uniprot_id = models.CharField(max_length=10,null=True,blank=True)
	uniprot_description = models.CharField(max_length=100,null=True,blank=True)
	top_nr_hit_eval = models.CharField(max_length=1000,null=True,blank=True)
	other_nr_hits = models.CharField(max_length=10000,null=True,blank=True)
	def __str__(self) :
		return self.nvertx_id		

class Standard_Error(models.Model) :
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
	name = models.CharField(max_length=20,primary_key=True)
	mfuzz_cluster_nb = models.IntegerField()
	cluster_image = models.CharField(max_length=50)
	bp_plot_image = models.CharField(max_length=50)
	def __str__(self) :
		return self.name















