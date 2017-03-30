from django.contrib import admin

from .models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Standard_Error, Mfuzz

admin.site.register(Fasta)
admin.site.register(Regen_cpm)
admin.site.register(Embryo_cpm)
admin.site.register(Annotation)
admin.site.register(Standard_Error)
admin.site.register(Mfuzz)