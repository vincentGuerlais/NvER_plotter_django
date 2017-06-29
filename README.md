# NvER_plotter_django
django-based NvER_plotter

NvERTx - A comparative embryogenesis & regeneration temporal gene expression plotter
Understanding the relationship between embryogenesis and regeneration is a long lasting question in regenerative biology as both developmental strategies lead to fully functional organisms. Modern functional genomics enables us today to re-address this question using integrative approaches that enable us to identify the similarities as well as importantly the regeneration specific elements that are potential candidates for regenerative medicine.

The starlet sea anemone Nematostella vectensis (Anthozoa, Cnidaria) is a unique model to study embryonic development and whole body regeneration in the same organism. In complement to an existing spatial gene expression database, we present here a novel tool to mine simultaneously (re)-analyzed/normalized temporal RNAseq datasets for embryonic development and oral regeneration of Nematostella vectensis. Data can be accessed by searching for your gene of interest using a gene name or the NvERTx ID of your gene of interest, scrolling trough the gene list corresponding to a given expression clusters, or by using the integrated BLAST ® function.

# Getting Started with Django

## system requirements
numpy==1.8.0rc1
Unix-like OS (we have run this on MacOSX, and SUSELinux)

## Other dependencies 
*Python 2.7 or higher 
*Apache/2.4 or higher (Other servers could work with minor tweaks.) 
*mod_wgsi 
*[BLAST executables 2.4 or higher](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 
### python packages
django 
biopython 


## Install Dependencies

First we install pip to make our life easier:

``` sh
sudo easy_install pip ;
pip --version ;
sudo pip install --upgrade pip
```

```
sudo pip install django ;
sudo pip install django-widget-tweaks ;
#might not need the blastplus
sudo pip install django-blastplus ;
sudo pip install biopython
```

Need NCBI dependencies

```
ftp ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ;
get ncbi-blast-2.6.0+-x64-macosx.tar.gz ;
quit ;
tar -zxvf ncbi-blast-2.6.0+-x64-macosx.tar.gz ;
sudo mv ncbi-blast-2.6.0+/bin/* /usr/local/bin/ ;
#test that they work
blastx -h
```

```
#might not need this
#python-dev -version
```

```
pip freeze
cd /Library/Python/2.7/site-packages/
#see if I can do this without cd-ing to the directory
sudo git clone https://github.com/vincentGuerlais/NvERtx_blastplus.git
sudo mv NvERtx_blastplus/blastplus/ ./

#modify this to point at the blast db!!!
vi /Library/Python/2.7/site-packages/blastplus/settings.py

```

# run the server

python manage.py runserver


##
##
## modifying the database
## 


python manage.py makemigrations ER_plotter
python manage.py migrate

## empty the database
python manage.py flush
rm db.sqlite3
rm ER_plotter/migrations/*

## repopulate
python manage.py makemigrations ER_plotter
python manage.py migrate

## fill the database
python manage.py shell
import DBfill
DBfill.DBFill('fileName','table')

## create admin
python manage.py createsuperuser
