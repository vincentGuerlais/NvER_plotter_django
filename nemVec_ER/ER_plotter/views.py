from django.http import HttpResponse
from django.shortcuts import render
from django.forms import formset_factory
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from ER_plotter.models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Standard_Error, Mfuzz
from .forms import Gene_searchForm, NvERTxForm, ConvertForm
import math, re, diggPaginator

def results(request):
	regen_point = [0,2,4,8,12,16,20,24,36,48,60,72,96,120,144]
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	#get id
	nvid = request.GET.getlist('Nvid', '')

	nvertx_1 = False
	nvertx_2 = False
	nvertx_3 = False
	nvertx_4 = False
	nvertx_5 = False
	log2 = False

	if nvid :
		nvertx_1 = nvid[0]
		if len(nvid) >= 2 :
			nvertx_2 = nvid[1]
		if len(nvid) >= 3 :
			nvertx_3 = nvid[2]
		if len(nvid) >= 4 :
			nvertx_4 = nvid[3]
		if len(nvid) == 5 :
			nvertx_5 = nvid[4]
		nvertx_search = True


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
		
	if nvertx_1 :
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
			nr_hit_graph_1 = re.search('\| (.+),', annot_top_nr_hit_eval_1).group(1) + ']'
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
			except :
				nvertx_2_regen_invalid = True
			try :
				embryo_warner_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_24HPF,2)
				embryo_warner_48_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_48HPF,2)
				embryo_warner_72_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_72HPF,2)
				embryo_warner_96_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_96HPF,2)
				embryo_warner_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_120HPF,2)
				embryo_warner_144_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_144HPF,2)
				embryo_warner_168_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_168HPF,2)
				embryo_warner_192_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_192HPF,2)
			except :
				nvertx_2_embryo_warner_invalid = True
			try :
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
			except :
				nvertx_2_embryo_fischer_invalid = True
			try :
				embryo_helm_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_2HPF,2)
				embryo_helm_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_7HPF,2)
				embryo_helm_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_12HPF,2)
				embryo_helm_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_24HPF,2)
				embryo_helm_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_120HPF,2)
				embryo_helm_240_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_240HPF,2)
			except :
				nvertx_2_embryo_helm_invalid = True
			try :
				embryo_mean_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_2HPF,2)
				embryo_mean_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_7HPF,2)
				embryo_mean_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_12HPF,2)
				embryo_mean_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_24HPF,2)
				embryo_mean_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_120HPF,2)
			except :
				nvertx_2_embryo_mean_invalid = True
			try :
				annot_nve_hit_2 = Annotation.objects.get(nvertx_id=nvertx_2).nve_hit
				annot_nve_eval_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).nve_eval,2)
				annot_mfuzz_clust_2 = Annotation.objects.get(nvertx_id=nvertx_2).mfuzz_clust
				annot_mfuzz_score_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).mfuzz_score,2)
				annot_uniprot_id_2 = Annotation.objects.get(nvertx_id=nvertx_2).uniprot_id
				annot_uniprot_description_2 = Annotation.objects.get(nvertx_id=nvertx_2).uniprot_description
				annot_top_nr_hit_eval_2 = Annotation.objects.get(nvertx_id=nvertx_2).top_nr_hit_eval
				annot_other_nr_hits_2 = Annotation.objects.get(nvertx_id=nvertx_2).other_nr_hits
				nr_hit_graph_2 = re.search('\| (.+),', annot_top_nr_hit_eval_2).group(1) + ']'
			except :
				nvertx_2_annot_invalid = True

		if nvertx_3 :
			try :
				sequence_fasta_3 = Fasta.objects.get(nvertx_id=nvertx_3).fasta_sequence
			except :
				sequence_fasta_3 = 'No fasta sequence'
			try :
				if log2 :
					regen_UC_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_UC+1,2),2)
					regen_0_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_0HPA+1,2),2)
					regen_2_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_2HPA+1,2),2)
					regen_4_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_4HPA+1,2),2)
					regen_8_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_8HPA+1,2),2)
					regen_12_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_12HPA+1,2),2)
					regen_16_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_16HPA+1,2),2)
					regen_20_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_20HPA+1,2),2)
					regen_24_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_24HPA+1,2),2)
					regen_36_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_36HPA+1,2),2)
					regen_48_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_48HPA+1,2),2)
					regen_60_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_60HPA+1,2),2)
					regen_72_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_72HPA+1,2),2)
					regen_96_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_96HPA+1,2),2)
					regen_120_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_120HPA+1,2),2)
					regen_144_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_144HPA+1,2),2)
					regen_se_UC_3 = [round(regen_UC_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_UC,2),round(regen_UC_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_UC,2)]
					regen_se_0_3 = [round(regen_0_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_0HPA,2),round(regen_0_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_0HPA,2)]
					regen_se_2_3 = [round(regen_2_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_2HPA,2),round(regen_2_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_2HPA,2)]
					regen_se_4_3 = [round(regen_4_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_4HPA,2),round(regen_4_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_4HPA,2)]
					regen_se_8_3 = [round(regen_8_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_8HPA,2),round(regen_8_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_8HPA,2)]
					regen_se_12_3 = [round(regen_12_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_12HPA,2),round(regen_12_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_12HPA,2)]
					regen_se_16_3 = [round(regen_16_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_16HPA,2),round(regen_16_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_16HPA,2)]
					regen_se_20_3 = [round(regen_20_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_20HPA,2),round(regen_20_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_20HPA,2)]
					regen_se_24_3 = [round(regen_24_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_24HPA,2),round(regen_24_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_24HPA,2)]
					regen_se_36_3 = [round(regen_36_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_36HPA,2),round(regen_36_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_36HPA,2)]
					regen_se_48_3 = [round(regen_48_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_48HPA,2),round(regen_48_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_48HPA,2)]
					regen_se_60_3 = [round(regen_60_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_60HPA,2),round(regen_60_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_60HPA,2)]
					regen_se_72_3 = [round(regen_72_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_72HPA,2),round(regen_72_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_72HPA,2)]
					regen_se_96_3 = [round(regen_96_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_96HPA,2),round(regen_96_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_96HPA,2)]
					regen_se_120_3 = [round(regen_120_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_120HPA,2),round(regen_120_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_120HPA,2)]
					regen_se_144_3 = [round(regen_144_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_144HPA,2),round(regen_144_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_log_se_144HPA,2)]
				else :
					regen_UC_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_UC,2)
					regen_0_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_0HPA,2)
					regen_2_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_2HPA,2)
					regen_4_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_4HPA,2)
					regen_8_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_8HPA,2)
					regen_12_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_12HPA,2)
					regen_16_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_16HPA,2)
					regen_20_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_20HPA,2)
					regen_24_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_24HPA,2)
					regen_36_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_36HPA,2)
					regen_48_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_48HPA,2)
					regen_60_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_60HPA,2)
					regen_72_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_72HPA,2)
					regen_96_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_96HPA,2)
					regen_120_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_120HPA,2)
					regen_144_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_144HPA,2)
					regen_se_UC_3 = [round(regen_UC_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_UC,2),round(regen_UC_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_UC,2)]
					regen_se_0_3 = [round(regen_0_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_0HPA,2),round(regen_0_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_0HPA,2)]
					regen_se_2_3 = [round(regen_2_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_2HPA,2),round(regen_2_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_2HPA,2)]
					regen_se_4_3 = [round(regen_4_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_4HPA,2),round(regen_4_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_4HPA,2)]
					regen_se_8_3 = [round(regen_8_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_8HPA,2),round(regen_8_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_8HPA,2)]
					regen_se_12_3 = [round(regen_12_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_12HPA,2),round(regen_12_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_12HPA,2)]
					regen_se_16_3 = [round(regen_16_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_16HPA,2),round(regen_16_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_16HPA,2)]
					regen_se_20_3 = [round(regen_20_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_20HPA,2),round(regen_20_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_20HPA,2)]
					regen_se_24_3 = [round(regen_24_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_24HPA,2),round(regen_24_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_24HPA,2)]
					regen_se_36_3 = [round(regen_36_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_36HPA,2),round(regen_36_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_36HPA,2)]
					regen_se_48_3 = [round(regen_48_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_48HPA,2),round(regen_48_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_48HPA,2)]
					regen_se_60_3 = [round(regen_60_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_60HPA,2),round(regen_60_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_60HPA,2)]
					regen_se_72_3 = [round(regen_72_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_72HPA,2),round(regen_72_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_72HPA,2)]
					regen_se_96_3 = [round(regen_96_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_96HPA,2),round(regen_96_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_96HPA,2)]
					regen_se_120_3 = [round(regen_120_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_120HPA,2),round(regen_120_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_120HPA,2)]
					regen_se_144_3 = [round(regen_144_3 - Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_144HPA,2),round(regen_144_3 + Standard_Error.objects.get(nvertx_id=nvertx_3).regen_se_144HPA,2)]
			except :
				nvertx_3_regen_invalid = True
			try :
				embryo_warner_24_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_24HPF,2)
				embryo_warner_48_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_48HPF,2)
				embryo_warner_72_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_72HPF,2)
				embryo_warner_96_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_96HPF,2)
				embryo_warner_120_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_120HPF,2)
				embryo_warner_144_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_144HPF,2)
				embryo_warner_168_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_168HPF,2)
				embryo_warner_192_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_192HPF,2)
			except :
				nvertx_3_embryo_warner_invalid = True
			try :
				embryo_fischer_0_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_0HPF,2)
				embryo_fischer_1_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_1HPF,2)
				embryo_fischer_2_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_2HPF,2)
				embryo_fischer_3_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_3HPF,2)
				embryo_fischer_4_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_4HPF,2)
				embryo_fischer_5_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_5HPF,2)
				embryo_fischer_6_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_6HPF,2)
				embryo_fischer_7_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_7HPF,2)
				embryo_fischer_8_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_8HPF,2)
				embryo_fischer_9_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_9HPF,2)	
				embryo_fischer_10_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_10HPF,2)
				embryo_fischer_11_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_11HPF,2)
				embryo_fischer_12_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_12HPF,2)
				embryo_fischer_13_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_13HPF,2)
				embryo_fischer_14_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_14HPF,2)
				embryo_fischer_15_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_15HPF,2)
				embryo_fischer_16_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_16HPF,2)		
				embryo_fischer_17_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_17HPF,2)
				embryo_fischer_18_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_18HPF,2)
				embryo_fischer_19_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_19HPF,2)
			except :
				nvertx_3_embryo_fischer_invalid = True
			try :
				embryo_helm_2_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_2HPF,2)
				embryo_helm_7_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_7HPF,2)
				embryo_helm_12_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_12HPF,2)
				embryo_helm_24_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_24HPF,2)
				embryo_helm_120_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_120HPF,2)
				embryo_helm_240_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_240HPF,2)
			except :
				nvertx_3_embryo_helm_invalid = True
			try :
				embryo_mean_2_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_2HPF,2)
				embryo_mean_7_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_7HPF,2)
				embryo_mean_12_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_12HPF,2)
				embryo_mean_24_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_24HPF,2)
				embryo_mean_120_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_120HPF,2)
			except :
				nvertx_3_embryo_mean_invalid = True
			try :
				annot_nve_hit_3 = Annotation.objects.get(nvertx_id=nvertx_3).nve_hit
				annot_nve_eval_3 = round(Annotation.objects.get(nvertx_id=nvertx_3).nve_eval,2)
				annot_mfuzz_clust_3 = Annotation.objects.get(nvertx_id=nvertx_3).mfuzz_clust
				annot_mfuzz_score_3 = round(Annotation.objects.get(nvertx_id=nvertx_3).mfuzz_score,2)
				annot_uniprot_id_3 = Annotation.objects.get(nvertx_id=nvertx_3).uniprot_id
				annot_uniprot_description_3 = Annotation.objects.get(nvertx_id=nvertx_3).uniprot_description
				annot_top_nr_hit_eval_3 = Annotation.objects.get(nvertx_id=nvertx_3).top_nr_hit_eval
				annot_other_nr_hits_3 = Annotation.objects.get(nvertx_id=nvertx_3).other_nr_hits
				nr_hit_graph_3 = re.search('\| (.+),', annot_top_nr_hit_eval_3).group(1) + ']'
			except :
				nvertx_3_annot_invalid = True

		if nvertx_4 :
			try :
				sequence_fasta_4 = Fasta.objects.get(nvertx_id=nvertx_4).fasta_sequence
			except :
				sequence_fasta_4 = 'No fasta sequence'
			try :
				if log2 :
					regen_UC_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_UC+1,2),2)
					regen_0_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_0HPA+1,2),2)
					regen_2_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_2HPA+1,2),2)
					regen_4_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_4HPA+1,2),2)
					regen_8_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_8HPA+1,2),2)
					regen_12_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_12HPA+1,2),2)
					regen_16_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_16HPA+1,2),2)
					regen_20_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_20HPA+1,2),2)
					regen_24_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_24HPA+1,2),2)
					regen_36_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_36HPA+1,2),2)
					regen_48_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_48HPA+1,2),2)
					regen_60_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_60HPA+1,2),2)
					regen_72_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_72HPA+1,2),2)
					regen_96_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_96HPA+1,2),2)
					regen_120_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_120HPA+1,2),2)
					regen_144_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_144HPA+1,2),2)
					regen_se_UC_4 = [round(regen_UC_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_UC,2),round(regen_UC_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_UC,2)]
					regen_se_0_4 = [round(regen_0_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_0HPA,2),round(regen_0_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_0HPA,2)]
					regen_se_2_4 = [round(regen_2_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_2HPA,2),round(regen_2_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_2HPA,2)]
					regen_se_4_4 = [round(regen_4_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_4HPA,2),round(regen_4_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_4HPA,2)]
					regen_se_8_4 = [round(regen_8_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_8HPA,2),round(regen_8_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_8HPA,2)]
					regen_se_12_4 = [round(regen_12_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_12HPA,2),round(regen_12_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_12HPA,2)]
					regen_se_16_4 = [round(regen_16_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_16HPA,2),round(regen_16_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_16HPA,2)]
					regen_se_20_4 = [round(regen_20_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_20HPA,2),round(regen_20_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_20HPA,2)]
					regen_se_24_4 = [round(regen_24_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_24HPA,2),round(regen_24_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_24HPA,2)]
					regen_se_36_4 = [round(regen_36_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_36HPA,2),round(regen_36_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_36HPA,2)]
					regen_se_48_4 = [round(regen_48_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_48HPA,2),round(regen_48_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_48HPA,2)]
					regen_se_60_4 = [round(regen_60_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_60HPA,2),round(regen_60_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_60HPA,2)]
					regen_se_72_4 = [round(regen_72_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_72HPA,2),round(regen_72_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_72HPA,2)]
					regen_se_96_4 = [round(regen_96_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_96HPA,2),round(regen_96_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_96HPA,2)]
					regen_se_120_4 = [round(regen_120_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_120HPA,2),round(regen_120_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_120HPA,2)]
					regen_se_144_4 = [round(regen_144_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_144HPA,2),round(regen_144_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_log_se_144HPA,2)]
				else :
					regen_UC_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_UC,2)
					regen_0_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_0HPA,2)
					regen_2_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_2HPA,2)
					regen_4_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_4HPA,2)
					regen_8_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_8HPA,2)
					regen_12_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_12HPA,2)
					regen_16_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_16HPA,2)
					regen_20_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_20HPA,2)
					regen_24_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_24HPA,2)
					regen_36_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_36HPA,2)
					regen_48_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_48HPA,2)
					regen_60_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_60HPA,2)
					regen_72_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_72HPA,2)
					regen_96_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_96HPA,2)
					regen_120_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_120HPA,2)
					regen_144_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_144HPA,2)
					regen_se_UC_4 = [round(regen_UC_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_UC,2),round(regen_UC_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_UC,2)]
					regen_se_0_4 = [round(regen_0_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_0HPA,2),round(regen_0_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_0HPA,2)]
					regen_se_2_4 = [round(regen_2_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_2HPA,2),round(regen_2_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_2HPA,2)]
					regen_se_4_4 = [round(regen_4_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_4HPA,2),round(regen_4_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_4HPA,2)]
					regen_se_8_4 = [round(regen_8_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_8HPA,2),round(regen_8_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_8HPA,2)]
					regen_se_12_4 = [round(regen_12_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_12HPA,2),round(regen_12_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_12HPA,2)]
					regen_se_16_4 = [round(regen_16_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_16HPA,2),round(regen_16_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_16HPA,2)]
					regen_se_20_4 = [round(regen_20_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_20HPA,2),round(regen_20_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_20HPA,2)]
					regen_se_24_4 = [round(regen_24_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_24HPA,2),round(regen_24_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_24HPA,2)]
					regen_se_36_4 = [round(regen_36_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_36HPA,2),round(regen_36_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_36HPA,2)]
					regen_se_48_4 = [round(regen_48_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_48HPA,2),round(regen_48_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_48HPA,2)]
					regen_se_60_4 = [round(regen_60_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_60HPA,2),round(regen_60_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_60HPA,2)]
					regen_se_72_4 = [round(regen_72_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_72HPA,2),round(regen_72_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_72HPA,2)]
					regen_se_96_4 = [round(regen_96_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_96HPA,2),round(regen_96_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_96HPA,2)]
					regen_se_120_4 = [round(regen_120_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_120HPA,2),round(regen_120_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_120HPA,2)]
					regen_se_144_4 = [round(regen_144_4 - Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_144HPA,2),round(regen_144_4 + Standard_Error.objects.get(nvertx_id=nvertx_4).regen_se_144HPA,2)]
			except :
				nvertx_4_regen_invalid = True
			try :
				embryo_warner_24_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_24HPF,2)
				embryo_warner_48_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_48HPF,2)
				embryo_warner_72_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_72HPF,2)
				embryo_warner_96_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_96HPF,2)
				embryo_warner_120_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_120HPF,2)
				embryo_warner_144_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_144HPF,2)
				embryo_warner_168_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_168HPF,2)
				embryo_warner_192_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_192HPF,2)
			except :
				nvertx_4_embryo_warner_invalid = True
			try :
				embryo_fischer_0_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_0HPF,2)
				embryo_fischer_1_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_1HPF,2)
				embryo_fischer_2_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_2HPF,2)
				embryo_fischer_3_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_3HPF,2)
				embryo_fischer_4_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_4HPF,2)
				embryo_fischer_5_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_5HPF,2)
				embryo_fischer_6_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_6HPF,2)
				embryo_fischer_7_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_7HPF,2)
				embryo_fischer_8_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_8HPF,2)
				embryo_fischer_9_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_9HPF,2)	
				embryo_fischer_10_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_10HPF,2)
				embryo_fischer_11_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_11HPF,2)
				embryo_fischer_12_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_12HPF,2)
				embryo_fischer_13_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_13HPF,2)
				embryo_fischer_14_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_14HPF,2)
				embryo_fischer_15_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_15HPF,2)
				embryo_fischer_16_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_16HPF,2)		
				embryo_fischer_17_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_17HPF,2)
				embryo_fischer_18_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_18HPF,2)
				embryo_fischer_19_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_19HPF,2)
			except :
				nvertx_4_embryo_fischer_invalid = True
			try :
				embryo_helm_2_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_2HPF,2)
				embryo_helm_7_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_7HPF,2)
				embryo_helm_12_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_12HPF,2)
				embryo_helm_24_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_24HPF,2)
				embryo_helm_120_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_120HPF,2)
				embryo_helm_240_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_240HPF,2)
			except :
				nvertx_4_embryo_helm_invalid = True
			try :
				embryo_mean_2_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_2HPF,2)
				embryo_mean_7_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_7HPF,2)
				embryo_mean_12_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_12HPF,2)
				embryo_mean_24_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_24HPF,2)
				embryo_mean_120_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_120HPF,2)
			except :
				nvertx_4_embryo_mean_invalid = True
			try :
				annot_nve_hit_4 = Annotation.objects.get(nvertx_id=nvertx_4).nve_hit
				annot_nve_eval_4 = round(Annotation.objects.get(nvertx_id=nvertx_4).nve_eval,2)
				annot_mfuzz_clust_4 = Annotation.objects.get(nvertx_id=nvertx_4).mfuzz_clust
				annot_mfuzz_score_4 = round(Annotation.objects.get(nvertx_id=nvertx_4).mfuzz_score,2)
				annot_uniprot_id_4 = Annotation.objects.get(nvertx_id=nvertx_4).uniprot_id
				annot_uniprot_description_4 = Annotation.objects.get(nvertx_id=nvertx_4).uniprot_description
				annot_top_nr_hit_eval_4 = Annotation.objects.get(nvertx_id=nvertx_4).top_nr_hit_eval
				annot_other_nr_hits_4 = Annotation.objects.get(nvertx_id=nvertx_4).other_nr_hits
				nr_hit_graph_4 = re.search('\| (.+),', annot_top_nr_hit_eval_4).group(1) + ']'
			except :
				nvertx_4_annot_invalid = True

		if nvertx_5 :
			try :
				sequence_fasta_5 = Fasta.objects.get(nvertx_id=nvertx_5).fasta_sequence
			except :
				sequence_fasta_5 = 'No fasta sequence'
			try :
				if log2 :
					regen_UC_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_UC+1,2),2)
					regen_0_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_0HPA+1,2),2)
					regen_2_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_2HPA+1,2),2)
					regen_4_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_4HPA+1,2),2)
					regen_8_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_8HPA+1,2),2)
					regen_12_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_12HPA+1,2),2)
					regen_16_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_16HPA+1,2),2)
					regen_20_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_20HPA+1,2),2)
					regen_24_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_24HPA+1,2),2)
					regen_36_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_36HPA+1,2),2)
					regen_48_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_48HPA+1,2),2)
					regen_60_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_60HPA+1,2),2)
					regen_72_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_72HPA+1,2),2)
					regen_96_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_96HPA+1,2),2)
					regen_120_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_120HPA+1,2),2)
					regen_144_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_144HPA+1,2),2)
					regen_se_UC_5 = [round(regen_UC_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_UC,2),round(regen_UC_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_UC,2)]
					regen_se_0_5 = [round(regen_0_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_0HPA,2),round(regen_0_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_0HPA,2)]
					regen_se_2_5 = [round(regen_2_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_2HPA,2),round(regen_2_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_2HPA,2)]
					regen_se_4_5 = [round(regen_4_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_4HPA,2),round(regen_4_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_4HPA,2)]
					regen_se_8_5 = [round(regen_8_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_8HPA,2),round(regen_8_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_8HPA,2)]
					regen_se_12_5 = [round(regen_12_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_12HPA,2),round(regen_12_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_12HPA,2)]
					regen_se_16_5 = [round(regen_16_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_16HPA,2),round(regen_16_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_16HPA,2)]
					regen_se_20_5 = [round(regen_20_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_20HPA,2),round(regen_20_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_20HPA,2)]
					regen_se_24_5 = [round(regen_24_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_24HPA,2),round(regen_24_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_24HPA,2)]
					regen_se_36_5 = [round(regen_36_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_36HPA,2),round(regen_36_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_36HPA,2)]
					regen_se_48_5 = [round(regen_48_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_48HPA,2),round(regen_48_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_48HPA,2)]
					regen_se_60_5 = [round(regen_60_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_60HPA,2),round(regen_60_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_60HPA,2)]
					regen_se_72_5 = [round(regen_72_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_72HPA,2),round(regen_72_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_72HPA,2)]
					regen_se_96_5 = [round(regen_96_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_96HPA,2),round(regen_96_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_96HPA,2)]
					regen_se_120_5 = [round(regen_120_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_120HPA,2),round(regen_120_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_120HPA,2)]
					regen_se_144_5 = [round(regen_144_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_144HPA,2),round(regen_144_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_log_se_144HPA,2)]
				else :
					regen_UC_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_UC,2)
					regen_0_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_0HPA,2)
					regen_2_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_2HPA,2)
					regen_4_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_4HPA,2)
					regen_8_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_8HPA,2)
					regen_12_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_12HPA,2)
					regen_16_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_16HPA,2)
					regen_20_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_20HPA,2)
					regen_24_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_24HPA,2)
					regen_36_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_36HPA,2)
					regen_48_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_48HPA,2)
					regen_60_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_60HPA,2)
					regen_72_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_72HPA,2)
					regen_96_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_96HPA,2)
					regen_120_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_120HPA,2)
					regen_144_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_144HPA,2)
					regen_se_UC_5 = [round(regen_UC_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_UC,2),round(regen_UC_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_UC,2)]
					regen_se_0_5 = [round(regen_0_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_0HPA,2),round(regen_0_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_0HPA,2)]
					regen_se_2_5 = [round(regen_2_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_2HPA,2),round(regen_2_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_2HPA,2)]
					regen_se_4_5 = [round(regen_4_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_4HPA,2),round(regen_4_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_4HPA,2)]
					regen_se_8_5 = [round(regen_8_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_8HPA,2),round(regen_8_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_8HPA,2)]
					regen_se_12_5 = [round(regen_12_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_12HPA,2),round(regen_12_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_12HPA,2)]
					regen_se_16_5 = [round(regen_16_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_16HPA,2),round(regen_16_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_16HPA,2)]
					regen_se_20_5 = [round(regen_20_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_20HPA,2),round(regen_20_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_20HPA,2)]
					regen_se_24_5 = [round(regen_24_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_24HPA,2),round(regen_24_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_24HPA,2)]
					regen_se_36_5 = [round(regen_36_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_36HPA,2),round(regen_36_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_36HPA,2)]
					regen_se_48_5 = [round(regen_48_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_48HPA,2),round(regen_48_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_48HPA,2)]
					regen_se_60_5 = [round(regen_60_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_60HPA,2),round(regen_60_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_60HPA,2)]
					regen_se_72_5 = [round(regen_72_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_72HPA,2),round(regen_72_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_72HPA,2)]
					regen_se_96_5 = [round(regen_96_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_96HPA,2),round(regen_96_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_96HPA,2)]
					regen_se_120_5 = [round(regen_120_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_120HPA,2),round(regen_120_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_120HPA,2)]
					regen_se_144_5 = [round(regen_144_5 - Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_144HPA,2),round(regen_144_5 + Standard_Error.objects.get(nvertx_id=nvertx_5).regen_se_144HPA,2)]
			except :
				nvertx_5_regen_invalid = True
			try :
				embryo_warner_24_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_24HPF,2)
				embryo_warner_48_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_48HPF,2)
				embryo_warner_72_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_72HPF,2)
				embryo_warner_96_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_96HPF,2)
				embryo_warner_120_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_120HPF,2)
				embryo_warner_144_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_144HPF,2)
				embryo_warner_168_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_168HPF,2)
				embryo_warner_192_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_192HPF,2)
			except :
				nvertx_5_embryo_warner_invalid = True
			try :
				embryo_fischer_0_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_0HPF,2)
				embryo_fischer_1_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_1HPF,2)
				embryo_fischer_2_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_2HPF,2)
				embryo_fischer_3_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_3HPF,2)
				embryo_fischer_4_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_4HPF,2)
				embryo_fischer_5_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_5HPF,2)
				embryo_fischer_6_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_6HPF,2)
				embryo_fischer_7_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_7HPF,2)
				embryo_fischer_8_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_8HPF,2)
				embryo_fischer_9_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_9HPF,2)	
				embryo_fischer_10_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_10HPF,2)
				embryo_fischer_11_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_11HPF,2)
				embryo_fischer_12_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_12HPF,2)
				embryo_fischer_13_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_13HPF,2)
				embryo_fischer_14_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_14HPF,2)
				embryo_fischer_15_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_15HPF,2)
				embryo_fischer_16_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_16HPF,2)		
				embryo_fischer_17_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_17HPF,2)
				embryo_fischer_18_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_18HPF,2)
				embryo_fischer_19_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_19HPF,2)
			except :
				nvertx_5_embryo_fischer_invalid = True
			try :
				embryo_helm_2_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_2HPF,2)
				embryo_helm_7_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_7HPF,2)
				embryo_helm_12_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_12HPF,2)
				embryo_helm_24_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_24HPF,2)
				embryo_helm_120_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_120HPF,2)
				embryo_helm_240_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_240HPF,2)
			except :
				nvertx_5_embryo_helm_invalid = True
			try :
				embryo_mean_2_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_2HPF,2)
				embryo_mean_7_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_7HPF,2)
				embryo_mean_12_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_12HPF,2)
				embryo_mean_24_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_24HPF,2)
				embryo_mean_120_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_120HPF,2)
			except :
				nvertx_5_embryo_mean_invalid = True
			try :
				annot_nve_hit_5 = Annotation.objects.get(nvertx_id=nvertx_5).nve_hit
				annot_nve_eval_5 = round(Annotation.objects.get(nvertx_id=nvertx_5).nve_eval,2)
				annot_mfuzz_clust_5 = Annotation.objects.get(nvertx_id=nvertx_5).mfuzz_clust
				annot_mfuzz_score_5 = round(Annotation.objects.get(nvertx_id=nvertx_5).mfuzz_score,2)
				annot_uniprot_id_5 = Annotation.objects.get(nvertx_id=nvertx_5).uniprot_id
				annot_uniprot_description_5 = Annotation.objects.get(nvertx_id=nvertx_5).uniprot_description
				annot_top_nr_hit_eval_5 = Annotation.objects.get(nvertx_id=nvertx_5).top_nr_hit_eval
				annot_other_nr_hits_5 = Annotation.objects.get(nvertx_id=nvertx_5).other_nr_hits
				nr_hit_graph_5 = re.search('\| (.+),', annot_top_nr_hit_eval_5).group(1) + ']'
			except :
				nvertx_5_annot_invalid = True

	return render(request, 'ER_plotter/nvertxResults.html', locals())

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
	return render(request, 'ER_plotter/mfuzz.html', locals())

### Mfuzz home page
def mfuzzHome(request):
	#Side panel forms
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	#list of clusters. Used to create the buttons
	clusters_list = Mfuzz.objects.all()

	return render(request, 'ER_plotter/mfuzzHome.html', locals())


### Mfuzz page when a cluster has been selected (button clicked)
def mfuzzResults(request,mfuzz_nb):
	#Side panel forms
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	#list of clusters. Used to create the buttons
	clusters_list = Mfuzz.objects.all()

	#Title and plots for the selected cluster
	mfuzz_cluster_nb = clusters_list.get(name='Mfuzz'+str(mfuzz_nb)).mfuzz_cluster_nb
	mfuzz_graph = clusters_list.get(name='Mfuzz'+str(mfuzz_nb)).cluster_image
	mfuzz_bp_plot = clusters_list.get(name='Mfuzz'+str(mfuzz_nb)).bp_plot_image

	#Details of the selected cluster
	cluster_list = Annotation.objects.filter(mfuzz_clust=mfuzz_nb)

	#pagination for the details
	page = request.GET.get('page', 1)
	paginator = diggPaginator.DiggPaginator(cluster_list, 50, body=5)
	try:
		cluster_table = paginator.page(page)
	except PageNotAnInteger:
		cluster_table = paginator.page(1)
	except EmptyPage:
		cluster_table = paginator.page(paginator.num_pages)

	return render(request, 'ER_plotter/mfuzzResults.html', locals())

def home(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	return render(request, 'ER_plotter/home.html', locals())

def searchResults(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	if gene_search_form.is_valid():
		search_query = gene_search_form.cleaned_data['gene_name']

	search_result_all = Annotation.objects.filter(nvertx_id__icontains=search_query) | Annotation.objects.filter(nve_hit__icontains=search_query) | Annotation.objects.filter(uniprot_id__icontains=search_query) | Annotation.objects.filter(uniprot_description__icontains=search_query) | Annotation.objects.filter(top_nr_hit_eval__icontains=search_query) | Annotation.objects.filter(other_nr_hits__icontains=search_query)

	#pagination for the details
	page = request.GET.get('page', 1)
	paginator = diggPaginator.DiggPaginator(search_result_all, 50, body=5)
	try:
		search_result = paginator.page(page)
	except PageNotAnInteger:
		search_result = paginator.page(1)
	except EmptyPage:
		search_result = paginator.page(paginator.num_pages)

	return render(request, 'ER_plotter/searchResults.html', locals())


def test(request):
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)

	cluster_list = Annotation.objects.filter(mfuzz_clust=1)
	paginator = diggPaginator.DiggPaginator(cluster_list, 10, body=6, padding=2)

	return render(request, 'ER_plotter/test.html', {'page': paginator.page(1)})

	#return render(request, 'ER_plotter/test.html', locals())



