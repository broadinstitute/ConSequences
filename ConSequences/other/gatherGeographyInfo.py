import os
import sys
from Bio import Entrez
import time
import argparse
from collections import defaultdict
from geopy.geocoders import Nominatim
geolocator = Nominatim(user_agent="gaherGeographyInfo")

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Use Entrez functionalities in Biopython and geopy to gather geographic origins of hits to NCBI's nt database.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--nt_hits_file', help='File listing Nucleotide sample IDs (hits from BLASTing to nt).', required=True)
	parser.add_argument('-m', '--mail', help='email to use for Entrez.', required=True)
	args = parser.parse_args()
	return args

# This program was made more general than the scripts used for Pironti and Salamzade et al. 2020

myargs = create_parser()
nt_hits_file = myargs.nt_hits_file
mail = myargs.mail

Entrez.email = mail

try:
	assert(os.path.isfile(nt_hits_file))
except:
	sys.stderr.write("ERROR: validating file exists!")

with open(nt_hits_file) as onhf:
	for line in onhf:
		nt_hit = line.strip()
		biosample = None
		try:
			handle = Entrez.elink(dbfrom="nucleotide", db='biosample', id=nt_hit)
			for l in handle.readlines():
				l = l.strip().decode('utf-8')
				if l.startswith("<Id>"):
					biosample = l.split('<Id>')[1].split('</Id>')[0]
			handle.close()
		except:
			sys.stderr.write("Issue with getting nucleotide db information for nt hit %s\n" % nt_hit)
		if biosample:
			lat = None;
			lon = None;
			geo_attributes = []
			try:
				handle = Entrez.efetch(db='biosample', id=biosample)
				l = handle.readlines()[1].decode('utf-8')
				handle.close()
				for it in l.split('</'):
					if '_name="latitude">' in it:
						lat = it.split('>')[-1].strip()
					if '_name="longitude">' in it:
						lon = it.split('>')[-1].strip()
					if '_name="geographic location">' in it:
						geo_attributes.append(it.split('>')[-1].strip())
					if '_name="latitude and longitude">' in it:
						lat, da, lon, do  = it.split('>')[-1].strip().split()
						if da == 'S': lat = '-' + lat
						if do == 'W': lon = '-' + lon
			except:
				sys.stderr.write("Issue with getting biosample db information for nt hit %s\n" % nt_hit)

			long_lat_support = []
			flag = False
			try:
				float(lat)
				float(lon)
				flag = True
			except:
				pass
			if flag:
				location = geolocator.reverse(lat + ', ' + lon)
				long_lat_support = []
				for desc in location.raw.values():
					long_lat_support.append(str(desc).lower())

			print(nt_hit + '\t' + '; '.join(long_lat_support) + '\t' + '; '.join(geo_attributes))
		time.sleep(1)
