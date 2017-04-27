from django import forms

class Gene_searchForm(forms.Form):
    gene_name = forms.CharField(max_length=20)

class NvERTxForm(forms.Form):
    nvertx_1 = forms.CharField(max_length=25)
    nvertx_2 = forms.CharField(max_length=25, required=False)
    nvertx_3 = forms.CharField(max_length=25, required=False)
    nvertx_4 = forms.CharField(max_length=25, required=False)
    nvertx_5 = forms.CharField(max_length=25, required=False)
    log2 = forms.BooleanField(initial=False, required=False)

class ConvertForm(forms.Form):
	nvertx = forms.CharField(max_length=25)