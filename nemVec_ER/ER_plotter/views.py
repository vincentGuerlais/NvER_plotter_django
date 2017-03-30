from django.http import HttpResponse
from django.shortcuts import render
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from ER_plotter.models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Standard_Error, Mfuzz
from .forms import Gene_searchForm, NvERTxForm
import math

def home(request):
	regen_point = [0,2,4,8,12,16,20,24,36,48,60,72,96,120,144]
	gene_search_form = Gene_searchForm(request.POST or None)
	if gene_search_form.is_valid():
		gene_name = gene_search_form.cleaned_data['gene_name']
		gene_search = True
	nvertx_form = NvERTxForm(request.POST or None)
	if nvertx_form.is_valid():
		nvertx_1 = nvertx_form.cleaned_data['nvertx_1']
		if nvertx_1[0] != 'N' :
			nvertx_1 = 'NvERTx.2.' + nvertx_1
		nvertx_2 = nvertx_form.cleaned_data['nvertx_2']
		if nvertx_2 and nvertx_2[0] != 'N' :
			nvertx_2 = 'NvERTx.2.' + nvertx_2
		nvertx_3 = nvertx_form.cleaned_data['nvertx_3']
		if nvertx_3 and nvertx_3[0] != 'N' :
			nvertx_3 = 'NvERTx.2.' + nvertx_3
		nvertx_4 = nvertx_form.cleaned_data['nvertx_4']
		if nvertx_4 and nvertx_4[0] != 'N' :
			nvertx_4 = 'NvERTx.2.' + nvertx_4
		nvertx_5 = nvertx_form.cleaned_data['nvertx_5']
		if nvertx_5 and nvertx_5[0] != 'N' :
			nvertx_5 = 'NvERTx.2.' + nvertx_5
		log2 = nvertx_form.cleaned_data['log2']
		nvertx_search = True
		try :
			sequence_fasta_1 = Fasta.objects.get(nvertx_id=nvertx_1).fasta_sequence
		except :
			sequence_fasta_1 = 'No fasta sequence'
		try :
			if log2 :
				regen_UC_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_UC+1,2),2)
				regen_0_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_0HPA+1,2),2)
				regen_2_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_2HPA+1,2),2)
				regen_4_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_4HPA+1,2),2)
				regen_8_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_8HPA+1,2),2)
				regen_12_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_12HPA+1,2),2)
				regen_16_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_16HPA+1,2),2)
				regen_20_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_20HPA+1,2),2)
				regen_24_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_24HPA+1,2),2)
				regen_36_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_36HPA+1,2),2)
				regen_48_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_48HPA+1,2),2)
				regen_60_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_60HPA+1,2),2)
				regen_72_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_72HPA+1,2),2)
				regen_96_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_96HPA+1,2),2)
				regen_120_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_120HPA+1,2),2)
				regen_144_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_144HPA+1,2),2)
				regen_se_UC_1 = [round(regen_UC_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_UC,2),round(regen_UC_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_UC,2)]
				regen_se_0_1 = [round(regen_0_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_0HPA,2),round(regen_0_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_0HPA,2)]
				regen_se_2_1 = [round(regen_2_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_2HPA,2),round(regen_2_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_2HPA,2)]
				regen_se_4_1 = [round(regen_4_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_4HPA,2),round(regen_4_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_4HPA,2)]
				regen_se_8_1 = [round(regen_8_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_8HPA,2),round(regen_8_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_8HPA,2)]
				regen_se_12_1 = [round(regen_12_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_12HPA,2),round(regen_12_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_12HPA,2)]
				regen_se_16_1 = [round(regen_16_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_16HPA,2),round(regen_16_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_16HPA,2)]
				regen_se_20_1 = [round(regen_20_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_20HPA,2),round(regen_20_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_20HPA,2)]
				regen_se_24_1 = [round(regen_24_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_24HPA,2),round(regen_24_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_24HPA,2)]
				regen_se_36_1 = [round(regen_36_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_36HPA,2),round(regen_36_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_36HPA,2)]
				regen_se_48_1 = [round(regen_48_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_48HPA,2),round(regen_48_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_48HPA,2)]
				regen_se_60_1 = [round(regen_60_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_60HPA,2),round(regen_60_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_60HPA,2)]
				regen_se_72_1 = [round(regen_72_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_72HPA,2),round(regen_72_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_72HPA,2)]
				regen_se_96_1 = [round(regen_96_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_96HPA,2),round(regen_96_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_96HPA,2)]
				regen_se_120_1 = [round(regen_120_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_120HPA,2),round(regen_120_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_120HPA,2)]
				regen_se_144_1 = [round(regen_144_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_144HPA,2),round(regen_144_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_log_se_144HPA,2)]
			else :
				regen_UC_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_UC,2)
				regen_0_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_0HPA,2)
				regen_2_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_2HPA,2)
				regen_4_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_4HPA,2)
				regen_8_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_8HPA,2)
				regen_12_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_12HPA,2)
				regen_16_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_16HPA,2)
				regen_20_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_20HPA,2)
				regen_24_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_24HPA,2)
				regen_36_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_36HPA,2)
				regen_48_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_48HPA,2)
				regen_60_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_60HPA,2)
				regen_72_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_72HPA,2)
				regen_96_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_96HPA,2)
				regen_120_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_120HPA,2)
				regen_144_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_144HPA,2)
				regen_se_UC_1 = [round(regen_UC_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_UC,2),round(regen_UC_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_UC,2)]
				regen_se_0_1 = [round(regen_0_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_0HPA,2),round(regen_0_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_0HPA,2)]
				regen_se_2_1 = [round(regen_2_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_2HPA,2),round(regen_2_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_2HPA,2)]
				regen_se_4_1 = [round(regen_4_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_4HPA,2),round(regen_4_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_4HPA,2)]
				regen_se_8_1 = [round(regen_8_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_8HPA,2),round(regen_8_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_8HPA,2)]
				regen_se_12_1 = [round(regen_12_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_12HPA,2),round(regen_12_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_12HPA,2)]
				regen_se_16_1 = [round(regen_16_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_16HPA,2),round(regen_16_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_16HPA,2)]
				regen_se_20_1 = [round(regen_20_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_20HPA,2),round(regen_20_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_20HPA,2)]
				regen_se_24_1 = [round(regen_24_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_24HPA,2),round(regen_24_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_24HPA,2)]
				regen_se_36_1 = [round(regen_36_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_36HPA,2),round(regen_36_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_36HPA,2)]
				regen_se_48_1 = [round(regen_48_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_48HPA,2),round(regen_48_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_48HPA,2)]
				regen_se_60_1 = [round(regen_60_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_60HPA,2),round(regen_60_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_60HPA,2)]
				regen_se_72_1 = [round(regen_72_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_72HPA,2),round(regen_72_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_72HPA,2)]
				regen_se_96_1 = [round(regen_96_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_96HPA,2),round(regen_96_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_96HPA,2)]
				regen_se_120_1 = [round(regen_120_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_120HPA,2),round(regen_120_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_120HPA,2)]
				regen_se_144_1 = [round(regen_144_1 - Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_144HPA,2),round(regen_144_1 + Standard_Error.objects.get(nvertx_id=nvertx_1).regen_se_144HPA,2)]
		except :
			nvertx_1_regen_invalid = True
		try :
			embryo_warner_24_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_24HPF,2)
			embryo_warner_48_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_48HPF,2)
			embryo_warner_72_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_72HPF,2)
			embryo_warner_96_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_96HPF,2)
			embryo_warner_120_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_120HPF,2)
			embryo_warner_144_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_144HPF,2)
			embryo_warner_168_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_168HPF,2)
			embryo_warner_192_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_192HPF,2)
		except :
			nvertx_1_embryo_warner_invalid = True
		try :
			embryo_fischer_0_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_0HPF,2)
			embryo_fischer_1_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_1HPF,2)
			embryo_fischer_2_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_2HPF,2)
			embryo_fischer_3_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_3HPF,2)
			embryo_fischer_4_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_4HPF,2)
			embryo_fischer_5_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_5HPF,2)
			embryo_fischer_6_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_6HPF,2)
			embryo_fischer_7_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_7HPF,2)
			embryo_fischer_8_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_8HPF,2)
			embryo_fischer_9_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_9HPF,2)	
			embryo_fischer_10_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_10HPF,2)
			embryo_fischer_11_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_11HPF,2)
			embryo_fischer_12_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_12HPF,2)
			embryo_fischer_13_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_13HPF,2)
			embryo_fischer_14_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_14HPF,2)
			embryo_fischer_15_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_15HPF,2)
			embryo_fischer_16_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_16HPF,2)		
			embryo_fischer_17_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_17HPF,2)
			embryo_fischer_18_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_18HPF,2)
			embryo_fischer_19_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_19HPF,2)
		except :
			nvertx_1_embryo_fischer_invalid = True
		try :
			embryo_helm_2_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_2HPF,2)
			embryo_helm_7_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_7HPF,2)
			embryo_helm_12_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_12HPF,2)
			embryo_helm_24_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_24HPF,2)
			embryo_helm_120_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_120HPF,2)
			embryo_helm_240_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_240HPF,2)
		except :
			nvertx_1_embryo_helm_invalid = True
		try :
			embryo_mean_2_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_2HPF,2)
			embryo_mean_7_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_7HPF,2)
			embryo_mean_12_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_12HPF,2)
			embryo_mean_24_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_24HPF,2)
			embryo_mean_120_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_120HPF,2)
		except :
			nvertx_1_embryo_mean_invalid = True
		try :
			annot_nve_hit_1 = Annotation.objects.get(nvertx_id=nvertx_1).nve_hit
			annot_nve_eval_1 = round(Annotation.objects.get(nvertx_id=nvertx_1).nve_eval,2)
			annot_mfuzz_clust_1 = Annotation.objects.get(nvertx_id=nvertx_1).mfuzz_clust
			annot_mfuzz_score_1 = round(Annotation.objects.get(nvertx_id=nvertx_1).mfuzz_score,2)
			annot_uniprot_id_1 = Annotation.objects.get(nvertx_id=nvertx_1).uniprot_id
			annot_uniprot_description_1 = Annotation.objects.get(nvertx_id=nvertx_1).uniprot_description
			annot_top_nr_hit_eval_1 = Annotation.objects.get(nvertx_id=nvertx_1).top_nr_hit_eval
			annot_other_nr_hits_1 = Annotation.objects.get(nvertx_id=nvertx_1).other_nr_hits
		except :
			nvertx_1_annot_invalid = True
		if nvertx_2 :
			try :
				sequence_fasta_2 = Fasta.objects.get(nvertx_id=nvertx_2).fasta_sequence
			except :
				sequence_fasta_2 = 'No fasta sequence'
			try :
				if log2 :
					regen_UC_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_UC+1,2),2)
					regen_0_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_0HPA+1,2),2)
					regen_2_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_2HPA+1,2),2)
					regen_4_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_4HPA+1,2),2)
					regen_8_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_8HPA+1,2),2)
					regen_12_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_12HPA+1,2),2)
					regen_16_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_16HPA+1,2),2)
					regen_20_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_20HPA+1,2),2)
					regen_24_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_24HPA+1,2),2)
					regen_36_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_36HPA+1,2),2)
					regen_48_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_48HPA+1,2),2)
					regen_60_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_60HPA+1,2),2)
					regen_72_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_72HPA+1,2),2)
					regen_96_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_96HPA+1,2),2)
					regen_120_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_120HPA+1,2),2)
					regen_144_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_144HPA+1,2),2)
					regen_se_UC_2 = [round(regen_UC_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_UC,2),round(regen_UC_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_UC,2)]
					regen_se_0_2 = [round(regen_0_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_0HPA,2),round(regen_0_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_0HPA,2)]
					regen_se_2_2 = [round(regen_2_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_2HPA,2),round(regen_2_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_2HPA,2)]
					regen_se_4_2 = [round(regen_4_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_4HPA,2),round(regen_4_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_4HPA,2)]
					regen_se_8_2 = [round(regen_8_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_8HPA,2),round(regen_8_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_8HPA,2)]
					regen_se_12_2 = [round(regen_12_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_12HPA,2),round(regen_12_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_12HPA,2)]
					regen_se_16_2 = [round(regen_16_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_16HPA,2),round(regen_16_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_16HPA,2)]
					regen_se_20_2 = [round(regen_20_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_20HPA,2),round(regen_20_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_20HPA,2)]
					regen_se_24_2 = [round(regen_24_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_24HPA,2),round(regen_24_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_24HPA,2)]
					regen_se_36_2 = [round(regen_36_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_36HPA,2),round(regen_36_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_36HPA,2)]
					regen_se_48_2 = [round(regen_48_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_48HPA,2),round(regen_48_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_48HPA,2)]
					regen_se_60_2 = [round(regen_60_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_60HPA,2),round(regen_60_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_60HPA,2)]
					regen_se_72_2 = [round(regen_72_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_72HPA,2),round(regen_72_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_72HPA,2)]
					regen_se_96_2 = [round(regen_96_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_96HPA,2),round(regen_96_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_96HPA,2)]
					regen_se_120_2 = [round(regen_120_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_120HPA,2),round(regen_120_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_120HPA,2)]
					regen_se_144_2 = [round(regen_144_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_144HPA,2),round(regen_144_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_log_se_144HPA,2)]
				else :
					regen_UC_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_UC,2)
					regen_0_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_0HPA,2)
					regen_2_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_2HPA,2)
					regen_4_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_4HPA,2)
					regen_8_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_8HPA,2)
					regen_12_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_12HPA,2)
					regen_16_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_16HPA,2)
					regen_20_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_20HPA,2)
					regen_24_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_24HPA,2)
					regen_36_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_36HPA,2)
					regen_48_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_48HPA,2)
					regen_60_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_60HPA,2)
					regen_72_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_72HPA,2)
					regen_96_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_96HPA,2)
					regen_120_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_120HPA,2)
					regen_144_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_144HPA,2)
					regen_se_UC_2 = [round(regen_UC_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_UC,2),round(regen_UC_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_UC,2)]
					regen_se_0_2 = [round(regen_0_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_0HPA,2),round(regen_0_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_0HPA,2)]
					regen_se_2_2 = [round(regen_2_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_2HPA,2),round(regen_2_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_2HPA,2)]
					regen_se_4_2 = [round(regen_4_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_4HPA,2),round(regen_4_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_4HPA,2)]
					regen_se_8_2 = [round(regen_8_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_8HPA,2),round(regen_8_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_8HPA,2)]
					regen_se_12_2 = [round(regen_12_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_12HPA,2),round(regen_12_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_12HPA,2)]
					regen_se_16_2 = [round(regen_16_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_16HPA,2),round(regen_16_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_16HPA,2)]
					regen_se_20_2 = [round(regen_20_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_20HPA,2),round(regen_20_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_20HPA,2)]
					regen_se_24_2 = [round(regen_24_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_24HPA,2),round(regen_24_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_24HPA,2)]
					regen_se_36_2 = [round(regen_36_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_36HPA,2),round(regen_36_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_36HPA,2)]
					regen_se_48_2 = [round(regen_48_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_48HPA,2),round(regen_48_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_48HPA,2)]
					regen_se_60_2 = [round(regen_60_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_60HPA,2),round(regen_60_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_60HPA,2)]
					regen_se_72_2 = [round(regen_72_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_72HPA,2),round(regen_72_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_72HPA,2)]
					regen_se_96_2 = [round(regen_96_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_96HPA,2),round(regen_96_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_96HPA,2)]
					regen_se_120_2 = [round(regen_120_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_120HPA,2),round(regen_120_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_120HPA,2)]
					regen_se_144_2 = [round(regen_144_2 - Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_144HPA,2),round(regen_144_2 + Standard_Error.objects.get(nvertx_id=nvertx_2).regen_se_144HPA,2)]
				embryo_warner_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_24HPF,2)
				embryo_warner_48_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_48HPF,2)
				embryo_warner_72_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_72HPF,2)
				embryo_warner_96_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_96HPF,2)
				embryo_warner_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_120HPF,2)
				embryo_warner_144_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_144HPF,2)
				embryo_warner_168_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_168HPF,2)
				embryo_warner_192_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_192HPF,2)
				embryo_fischer_0_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_0HPF,2)
				embryo_fischer_1_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_1HPF,2)
				embryo_fischer_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_2HPF,2)
				embryo_fischer_3_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_3HPF,2)
				embryo_fischer_4_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_4HPF,2)
				embryo_fischer_5_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_5HPF,2)
				embryo_fischer_6_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_6HPF,2)
				embryo_fischer_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_7HPF,2)
				embryo_fischer_8_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_8HPF,2)
				embryo_fischer_9_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_9HPF,2)	
				embryo_fischer_10_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_10HPF,2)
				embryo_fischer_11_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_11HPF,2)
				embryo_fischer_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_12HPF,2)
				embryo_fischer_13_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_13HPF,2)
				embryo_fischer_14_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_14HPF,2)
				embryo_fischer_15_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_15HPF,2)
				embryo_fischer_16_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_16HPF,2)		
				embryo_fischer_17_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_17HPF,2)
				embryo_fischer_18_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_18HPF,2)
				embryo_fischer_19_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_19HPF,2)
				embryo_helm_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_2HPF,2)
				embryo_helm_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_7HPF,2)
				embryo_helm_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_12HPF,2)
				embryo_helm_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_24HPF,2)
				embryo_helm_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_120HPF,2)
				embryo_helm_240_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_240HPF,2)
				embryo_mean_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_2HPF,2)
				embryo_mean_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_7HPF,2)
				embryo_mean_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_12HPF,2)
				embryo_mean_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_24HPF,2)
				embryo_mean_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_120HPF,2)
				annot_nve_hit_2 = Annotation.objects.get(nvertx_id=nvertx_2).nve_hit
				annot_nve_eval_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).nve_eval,2)
				annot_mfuzz_clust_2 = Annotation.objects.get(nvertx_id=nvertx_2).mfuzz_clust
				annot_mfuzz_score_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).mfuzz_score,2)
				annot_uniprot_id_2 = Annotation.objects.get(nvertx_id=nvertx_2).uniprot_id
				annot_uniprot_description_2 = Annotation.objects.get(nvertx_id=nvertx_2).uniprot_description
				annot_top_nr_hit_eval_2 = Annotation.objects.get(nvertx_id=nvertx_2).top_nr_hit_eval
				annot_other_nr_hits_2 = Annotation.objects.get(nvertx_id=nvertx_2).other_nr_hits
			except :
				nvertx_search_2_invalid = True
	return render(request, 'ER_plotter/home.html', locals())

def resultats(request,nvertx_id):
	return render(request, 'ER_plotter/resultats.html',
		{'sequence_fasta': Fasta.objects.get(nvertx_id=nvertx_id).fasta_sequence,
		'nvertx_id':nvertx_id,
		'regen_UC':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_UC,2),
		'regen_0':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_0HPA,2),
		'regen_2':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_2HPA,2),
		'regen_4':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_4HPA,2),
		'regen_8':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_8HPA,2),
		'regen_12':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_12HPA,2),
		'regen_16':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_16HPA,2),
		'regen_20':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_20HPA,2),
		'regen_24':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_24HPA,2),
		'regen_36':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_36HPA,2),
		'regen_48':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_48HPA,2),
		'regen_60':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_60HPA,2),
		'regen_72':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_72HPA,2),
		'regen_96':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_96HPA,2),
		'regen_120':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_120HPA,2),
		'regen_144':round(Regen_cpm.objects.get(nvertx_id=nvertx_id).regen_anc_144HPA,2),
		'embryo_warner_24':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_24HPF,2),
		'embryo_warner_48':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_48HPF,2),
		'embryo_warner_72':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_72HPF,2),
		'embryo_warner_96':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_96HPF,2),
		'embryo_warner_120':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_120HPF,2),
		'embryo_warner_144':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_144HPF,2),
		'embryo_warner_168':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_168HPF,2),
		'embryo_warner_192':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).warner_anc_192HPF,2),
		'embryo_fischer_0':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_0HPF,2),
		'embryo_fischer_1':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_1HPF,2),
		'embryo_fischer_2':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_2HPF,2),
		'embryo_fischer_3':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_3HPF,2),
		'embryo_fischer_4':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_4HPF,2),
		'embryo_fischer_5':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_5HPF,2),
		'embryo_fischer_6':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_6HPF,2),
		'embryo_fischer_7':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_7HPF,2),
		'embryo_fischer_8':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_8HPF,2),
		'embryo_fischer_9':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_9HPF,2),	
		'embryo_fischer_10':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_10HPF,2),
		'embryo_fischer_11':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_11HPF,2),
		'embryo_fischer_12':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_12HPF,2),
		'embryo_fischer_13':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_13HPF,2),
		'embryo_fischer_14':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_14HPF,2),
		'embryo_fischer_15':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_15HPF,2),
		'embryo_fischer_16':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_16HPF,2),		
		'embryo_fischer_17':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_17HPF,2),
		'embryo_fischer_18':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_18HPF,2),
		'embryo_fischer_19':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).fischer_anc_19HPF,2),
		'embryo_helm_2':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).helm_anc_2HPF,2),
		'embryo_helm_7':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).helm_anc_7HPF,2),
		'embryo_helm_12':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).helm_anc_12HPF,2),
		'embryo_helm_24':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).helm_anc_24HPF,2),
		'embryo_helm_120':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).helm_anc_120HPF,2),
		'embryo_helm_240':round(Embryo_cpm.objects.get(nvertx_id=nvertx_id).helm_anc_240HPF,2),
		'annot_nve_hit':Annotation.objects.get(nvertx_id=nvertx_id).nve_hit,
		'annot_nve_eval':round(Annotation.objects.get(nvertx_id=nvertx_id).nve_eval,2),
		'annot_mfuzz_clust':Annotation.objects.get(nvertx_id=nvertx_id).mfuzz_clust,
		'annot_mfuzz_score':round(Annotation.objects.get(nvertx_id=nvertx_id).mfuzz_score,2),
		'annot_uniprot_id':Annotation.objects.get(nvertx_id=nvertx_id).uniprot_id,
		'annot_uniprot_description':Annotation.objects.get(nvertx_id=nvertx_id).uniprot_description,
		'annot_top_nr_hit_eval':Annotation.objects.get(nvertx_id=nvertx_id).top_nr_hit_eval,
		'annot_other_nr_hits':Annotation.objects.get(nvertx_id=nvertx_id).other_nr_hits
		}
		)

def mfuzz(request):
	gene_search_form = Gene_searchForm(request.POST or None)
	if gene_search_form.is_valid():
		gene_name = gene_search_form.cleaned_data['gene_name']
		gene_search = True
	nvertx_form = NvERTxForm(request.POST or None)
	if nvertx_form.is_valid():
		nvertx_1 = nvertx_form.cleaned_data['nvertx_1']
		if nvertx_1[0] != 'N' :
			nvertx_1 = 'NvERTx.2.' + nvertx_1
		nvertx_2 = nvertx_form.cleaned_data['nvertx_2']
		if nvertx_2 and nvertx_2[0] != 'N' :
			nvertx_2 = 'NvERTx.2.' + nvertx_2
		nvertx_3 = nvertx_form.cleaned_data['nvertx_3']
		if nvertx_3 and nvertx_3[0] != 'N' :
			nvertx_3 = 'NvERTx.2.' + nvertx_3
		nvertx_4 = nvertx_form.cleaned_data['nvertx_4']
		if nvertx_4 and nvertx_4[0] != 'N' :
			nvertx_4 = 'NvERTx.2.' + nvertx_4
		nvertx_5 = nvertx_form.cleaned_data['nvertx_5']
		if nvertx_5 and nvertx_5[0] != 'N' :
			nvertx_5 = 'NvERTx.2.' + nvertx_5
		log2 = nvertx_form.cleaned_data['log2']
		nvertx_search = True
	clusters_list = Mfuzz.objects.all()
	mfuzz_all = Annotation.objects.all()
	mfuzz1_all = mfuzz_all.filter(mfuzz_clust=1)
	mfuzz2_all = mfuzz_all.filter(mfuzz_clust=2)
	mfuzz3_all = mfuzz_all.filter(mfuzz_clust=3)
	mfuzz4_all = mfuzz_all.filter(mfuzz_clust=4)
	mfuzz5_all = mfuzz_all.filter(mfuzz_clust=5)
	mfuzz6_all = mfuzz_all.filter(mfuzz_clust=6)
	mfuzz7_all = mfuzz_all.filter(mfuzz_clust=7)
	mfuzz8_all = mfuzz_all.filter(mfuzz_clust=8)
	mfuzz9_all = mfuzz_all.filter(mfuzz_clust=9)
	mfuzz10_all = mfuzz_all.filter(mfuzz_clust=10)
	mfuzz11_all = mfuzz_all.filter(mfuzz_clust=11)
	mfuzz12_all = mfuzz_all.filter(mfuzz_clust=12)
	return render(request, 'ER_plotter/mfuzz.html', locals())

def test(request):
	mfuzz1_all = Annotation.objects.filter(mfuzz_clust=1)
	return render(request, 'ER_plotter/test.html', locals())



