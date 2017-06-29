from django.conf.urls import url

from blastplus import views

urlpatterns = [
    url(r'^blastn/', views.blastn),
    url(r'^tblastn/', views.tblastn, name='tblastn'),
    url(r'^blast/$', views.blastn, name='blast'),
    url(r'^blastp/$', views.blastp),
    url(r'^blastx/$', views.blastx),

]
