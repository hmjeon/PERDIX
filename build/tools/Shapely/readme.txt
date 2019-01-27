Anaconda 2(64bit)
Add path
Integrate python 2.7

Requirements:
•	shapely
•	numpy
•	pyinstaller

Environment:
Anaconda 64bit + Python 2.7

Install:
conda install shapely
pyinstaller.exe -F Spec_Conv.spec Shapely.py

Etc:
Unofficial Windows Binaries for Python Extension Packages – wheel files
https://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely
conda config --add channels conda-forge


############################################

Linux
Install pip for python
wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py --user
pip install shapely