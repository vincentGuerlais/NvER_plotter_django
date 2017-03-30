from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^home$', views.home, name='home'),
    url(r'^res/(NvERTx.2.\d+)$', views.resultats, name='resultats'),
    url(r'^mfuzz$', views.mfuzz, name='mfuzz'),
    url(r'^test$', views.test, name='test'),
]