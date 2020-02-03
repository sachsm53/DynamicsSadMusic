#!/usr/bin/python

import os
import sys
import csv
import pandas as pd
import numpy as np

reload(sys)
sys.setdefaultencoding('utf8')

#set paths
pathbase = "/Volumes/MusicProject/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-1/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-2/NaturalisticNetwork"

csvfile = "%s/survey_post_fMRI_raw-02-22-19.csv" %pathbase

#BLOCKs DESIGN
# 1. 3 - Additional Sad Music Questions
# 2. 4 - Aesthetic Response Music
# 3. 5 - iri
# 4. 6 - PHQ
# 5. 7 - GAD (anxiety)
# 5. 8 - Gold-MSI
# 7. 9 - Gold MSI part 2
# 8. 10 - Demographics 1 
# 9. 11 - Demographics 2 (11.2 - age)


print "Opening Survey %s" %csvfile
reader = pd.read_csv(csvfile, skiprows = [1])


#Background Info
reader = reader.rename(columns={'Q11.2':'age', 'Q11.3':'sex', 'Q11.8':'language'})


#Block 3 Like Sad music
lsm = reader[[elem for elem in reader.columns if 'Q3.' in elem]]
lsm = lsm.drop(lsm.columns[[range(0,4)]], axis=1)


#Block 5: Interpersonal Reactivity Index
iri = reader[[elem for elem in reader.columns if 'Q5' in elem]] - 1
iri = iri.drop(iri.columns[[range(0,4)]], axis=1)

reverse_items = ['Q5.5','Q5.6','Q5.9','Q5.15','Q5.16','Q5.17','Q5.19','Q5.22','Q5.23']
for items in reverse_items:
	iri[items] = (4 - iri[items])

iri["fantasy"] = iri[['Q5.3','Q5.7','Q5.9','Q5.15','Q5.20','Q5.28','Q5.32']].sum(axis=1)
iri["perspective"] = iri[['Q5.5','Q5.11','Q5.14','Q5.19','Q5.25','Q5.31','Q5.34']].sum(axis=1)
iri["empathic"] = iri[['Q5.4','Q5.6','Q5.12','Q5.17','Q5.22','Q5.24','Q5.27']].sum(axis=1)
iri["pers_distress"] = iri[['Q5.8','Q5.13','Q5.16','Q5.21','Q5.23','Q5.30','Q5.33']].sum(axis=1)

iri_final = iri[["fantasy","perspective","empathic","pers_distress"]]
iri_final["global"] = iri[["fantasy","perspective","empathic","pers_distress"]].mean(axis=1)

#GOLDSMITH
# 5. 8 - Gold-MSI
# 7. 9 - Gold MSI part 2


#Output File
outcsv = os.path.join(pathbase,'naturalistic_network_survey_results.csv')
reader.to_csv(path_or_buf = outcsv)


