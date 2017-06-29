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
django-widget-tweaks  
biopython  

## Installation  
#### This tutorial is for MacOSX. For SUSELinux scroll down.

First we install pip to make our life easier:  
``` sh
sudo easy_install pip ;
pip --version ;
sudo pip install --upgrade pip
```
Next we install the python packages
```
sudo pip install django ;
sudo pip install django-widget-tweaks ;
sudo pip install biopython
```
Next we get the NCBI dependencies
```
ftp ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ;
get ncbi-blast-2.6.0+-x64-macosx.tar.gz ;
quit ;
tar -zxvf ncbi-blast-2.6.0+-x64-macosx.tar.gz ;
sudo mv ncbi-blast-2.6.0+/bin/* /usr/local/bin/ ;
#test that they work
blastx -h
```
Now we install the Django app
```
sudo git clone https://github.com/vincentGuerlais/NvER_plotter_django
#run it
cd NvER_plotter_django/nemVec_ER/
python manage.py
```
Get the databases one:
```
scp -v db.sqlite3 user@host:NvER_plotter_django/nemVec_ER/db.sqlite3
scp -v -r nemVec_ER/blast_db user@host:/NvER_plotter_django/nemVec_ER/
```
Configure the server:  
```
# MacOSX has apache2 installed already
# Check the version
sudo httpd -v
```
Install wgsi
```
pip install apache2-mod_wsgi
```
Edit the httpd file:
```
cp /etc/apache2/httpd.conf /etc/apache2/httpd.bak
sudo vi /etc/apache2/httpd.conf
#######add these lines to the bottom of ~/../../etc/apache2/httpd.conf
#this serves the static content
Alias /static/ /home/ec2-user/NvER_plotter_django/nemVec_ER/ER_plotter/static/
<Directory /home/ec2-user/NvER_plotter_django/nemVec_ER/ER_plotter/static>
Require all granted
</Directory>
#this serves the site
WSGIScriptAlias / /home/ec2-user/NvER_plotter_django/nemVec_ER/nemVec_ER/wsgi.py 
WSGIPythonPath /home/ec2-user/NvER_plotter_django/nemVec_ER/
<Directory /home/ec2-user/NvER_plotter_django/nemVec_ER/>
<Files wsgi.py>
Require all granted
</Files>
</Directory> 
######## end edit to ~/../../etc/apache2/httpd.conf
```
Start the server
```
apachectl start
```

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
