from django.contrib import admin

from .models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Regen_SE, Mfuzz, Embryo_SE, Regen_log_SE

admin.site.register(Fasta)
admin.site.register(Regen_cpm)
admin.site.register(Embryo_cpm)
admin.site.register(Annotation)
admin.site.register(Regen_SE)
admin.site.register(Mfuzz)
admin.site.register(Embryo_SE)
admin.site.register(Regen_log_SE)