# --------------------------------------------------------------------------------------------- 
# Script to make overlay cutouts:
#
#   -- SDSS:          greyscale, red circles for positions, optionally plot index numbers
#   -- Hodge/FIRST:   blue contour, blue + for positions
#   -- JVLA Stripe82: green contour, green + for positions
#   -- DR7 Quasars:   cyan circles for positions, optionally plot index numbers
#
# --------------------------------------------------------------------------------------------- 

# use backend Agg so we can run this script without X
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from astropy.io.votable import parse
from astropy.io import fits
import aplpy

import urllib, urllib2
import os

# to disable the aplpy logging - its too much! - take this out if you want to debug
import logging
logging.disable(logging.CRITICAL)
# supress annoying unicode warnings - take this out if you want to debug
logging.captureWarnings(True)

# --------------------------------------------------------------------------------------------- 
# Parameters

# Numer of processes to parallelise over
NPROCESSES = 0             # integer. set this to zero if you don't want to use multiprocessing

# base directory
BASEDIR = '/data2/matt/Overlaycode'

# source positions
POSNSFILE = 'JVLAsdssfilenameWEST.dat'
posfile = '{0}/{1}'.format(BASEDIR,POSNSFILE)

# FITS file locations
SDSS_DIR = '/data2/matt/attempt_2_unnested'
JVLA_DIR = '/data2/matt/JVLAcutouts/Westcutoutslarger'

# catalogue files
CATALOGUE1 = 'Hodge.vot'
CATALOGUE2 = 'West.vot'
CATALOGUE3 = 'MATCHESWITHIN40.vot'
CATALOGUE4 = 'DR7QuasarsBEST.vot'

# suffix for the files names (if desired)
SUFFIX = ''

# cutout sizes
CUTOUTSIZE = 1.0           # arcminutes
CUTOUTSIZEVLA = 3.0        # Needs to be about 2 arcminutes bigger than CUTOUTSIZE                            
CUTOUTSIZE_DEG = CUTOUTSIZE/60.

# SDSS circle and index name parameters
SDSS_CIRCLE = 1.0          # arcseconds 
SDSS_INDEX = True          # display index names?
SDSS_OFFSET = 4.3          # offset of index names from positions  
SDSS_CIRCLE_DEG = SDSS_CIRCLE/60./60.
SDSS_OFFSET_DEG = SDSS_OFFSET/60./60.

# Quasar circle parameters
QUASAR_CIRCLE = 1.4        # arcseconds 
QUASAR_INDEX = True        # display index names?
QUASAR_OFFSET = 3.8        # offset of index names from positions  
QUASAR_CIRCLE_DEG = QUASAR_CIRCLE/60./60.
QUASAR_OFFSET_DEG = QUASAR_OFFSET/60./60.

# --------------------------------------------------------------------------------------------- 
# Functions

def get_url_fits(query_dict, url, filename, err_string='No data is available'):
    # get FITS cutout from given survey url

    # query url 
    params = urllib.urlencode(query_dict)
    req = urllib2.Request(url,params)
    response = urllib2.urlopen(req)

    # put result into file, line by line
    with open(filename, 'wb') as response_file:
        line = response.read()
        if not err_string in line:
            response_file.write(line)
        else:
            # if there is no data available, return False
            os.system('rm {0}'.format(filename))
            return False
    return True

def get_hodge(query_dict, url, filename):
    return get_url_fits(query_dict, url, filename, err_string='No Stripe 82 data is available')

def get_first(query_dict, url, filename):
    return get_url_fits(query_dict, url, filename, err_string='No FIRST data is available')

def get_background_variance(data,sigma_clip=5.0,tolerance=0.01):
    """Compute the variance by iteratively removing outliers greater than a given sigma
    until the mean changes by no more than tolerance.
    (function stolen from TomM)

    Inputs
    ------
    data - 1d numpy array of data to compute variance
    sigma_clip - the amount of sigma to clip the data before the next iteration
    tolerance - the fractional change in the mean to stop iterating

    Outputs
    -------
    variance - the final background variance in the sigma clipped image
    """
    #Initialise diff and data_clip and mean and std
    diff = 1
    mean = np.nanmean(data)
    data_clip = data
    while diff > tolerance:
        data_clip = data_clip[np.abs(data_clip)<mean+sigma_clip*np.nanstd(data_clip)]
        newmean = np.nanmean(data_clip)
        diff = np.abs(mean-newmean)/(mean+newmean)
        mean = newmean
    return np.nanvar(data_clip)

# --------------------------------------------------------------------------------------------- 
# main cutout function
def make_cutout(target_parameters):
	"""
	Make a png cutout of the various images, in greyscale and contour, and catalogue position markers.

	Inputs
	------
	target parameters - list containing parameters for the cutout
	   target_nname - name of the target source, string
	   ra_string - Right Ascension of the target, string
	   dec_string - Declination of the target, string
	   field - SDSS field to use at the target position
	"""

	target_name, ra_string, dec_string, field = target_parameters

	ra = float(ra_string)
	dec = float(dec_string)

	print ' -- ', target_name, ra, dec, field

	# ---------------------------------------------------------------------------------------------
	# set up figure
	fig = plt.figure(figsize=(14, 14))

	# ---------------------------------------------------------------------------------------------
	# sdss greyscale
	sdss_fitsfile = '{0}/{1}'.format(SDSS_DIR,field)

	f = aplpy.FITSFigure(sdss_fitsfile, figure=fig, dimensions=[1,0])
	f.show_grayscale() 

	f.set_title(target_name)

	# ---------------------------------------------------------------------------------------------
	# Hodge/FIRST contour 

	# download fits file  
	#  - first try Hodge, then FIRST if Hodges is not available
	hodge_file, first_file = False, False

	# quirk of the interface is that the full position string goes into the RA field
	ra_hours = 24.*ra/360.
	pos_string = '{0} {1}'.format(ra_hours, dec)
	query_dict['RA'] = pos_string

	# cutout name
	fitsfile = '{0}_stripe82_cutout.fits'.format(target_name,) 

	# query for Hodge cutout
	hodge_file = get_hodge(query_dict, HODGE_URL, fitsfile)
	# if Hodge cutout unavailable, query for FIRST cutout
	if not hodge_file:
	    # search for FIRST cutout
	    first_file = get_first(query_dict, FIRST_URL, fitsfile)
	    if first_file:
	        print 'Using FIRST rather than Hodges for Stripe82 cutout'
	    else: 
	        print 'Neither Hodge nor FIRST data available for this field'

	if (first_file or hodge_file):
	    # if we have either Hodge of FIRST cutout downloaded, plot contours for it

	    # get rid of vestigial image header parameters
	    HF_fits = fits.open(fitsfile, mode='update')
	    HF_header = HF_fits[0].header
	    remove_list = ['CTYPE3','CRVAL3','CDELT3','CRPIX3','CROTA3','CTYPE4','CRVAL4','CDELT4','CRPIX4','CROTA4']
	    for key in remove_list:
	        HF_header.remove(key)
	    HF_fits.close()

	    # get data variance for calculating contour levels
	    fitsimage = fits.open(fitsfile)
	    data = fitsimage[0].data
	    sigma = np.sqrt(get_background_variance(data.flatten()))
	    HF_levels = np.logspace(np.log10(4.0*sigma),np.log10(0.9*np.nanmax(data)),num=5)  

	    # plot Hodge/FIRST contours
	    f.show_contour(fitsfile,levels=HF_levels,colors='b',overlap=True,smooth=1)

	# delete the temporary cutout file
	os.system('rm {0}'.format(fitsfile,))

	# ---------------------------------------------------------------------------------------------
	# JVLA contour
	vla_fitsfile = '{0}/{1}.fits'.format(JVLA_DIR,target_name)

	f.recenter(ra, dec, radius=CUTOUTSIZE/60.)  # radius is in degrees

	# get data varience for calculating contour levels
	fitsimage = fits.open(vla_fitsfile)
	data = fitsimage[0].data
	sigma = np.sqrt(get_background_variance(data.flatten()))
	JVLA_levels = np.logspace(np.log10(4.0*sigma),np.log10(0.9*np.nanmax(data)),num=10)  

	f.show_contour(vla_fitsfile,levels=JVLA_levels,colors='g',overlap=True,smooth=1)

	# ---------------------------------------------------------------------------------------------
	# catalogue 1 - Hodge source locations

	f.show_markers(r_hodge,d_hodge,edgecolor='b',marker='+',s=300,linewidths=3)

	# ---------------------------------------------------------------------------------------------
	# catalogue 2 - JVLA source locations

	f.show_markers(r_vla,d_vla,edgecolor='g',marker='+',s=300,linewidths=3)

	# ---------------------------------------------------------------------------------------------
	# catalogue 3 - SDSS source circles

	# only extract data in the cutout area 
	#  takes too long otherwise - the SDSS vot is huge!
	cutout_mask = (np.abs(r_sdss-ra)<CUTOUTSIZE_DEG)&(np.abs(d_sdss-dec)<CUTOUTSIZE_DEG)

	r_cutout = r_sdss[cutout_mask]
	d_cutout = d_sdss[cutout_mask]

	s = np.ones(len(r_cutout))*SDSS_CIRCLE_DEG
	f.show_circles(r_cutout,d_cutout,s,edgecolor='r')

	if SDSS_INDEX:
	    indices = index_sdss[cutout_mask]
	    for ri, di, index in zip(r_cutout, d_cutout,indices):
	        f.add_label(ri+SDSS_OFFSET_DEG,di,index,color='r',size=12)

	# ---------------------------------------------------------------------------------------------
	# catalogue 4 - Quasar source circles

	cutout_mask = (np.abs(r_quasar-ra)<CUTOUTSIZE_DEG)&(np.abs(d_quasar-dec)<CUTOUTSIZE_DEG)

	r_cutout = r_quasar[cutout_mask]
	d_cutout = d_quasar[cutout_mask]

	# if there are quasars in this cutout, plot them
	if len(r_cutout) > 0:
	    s = np.ones(len(r_cutout))*QUASAR_CIRCLE_DEG
	    f.show_circles(r_cutout,d_cutout,s,edgecolor='c')

	    if QUASAR_INDEX:
	        indices = index_quasar[cutout_mask]
	        for ri, di, index in zip(r_cutout, d_cutout,indices):
	            f.add_label(ri,di+QUASAR_OFFSET_DEG,index,color='c',size=12)

	# ---------------------------------------------------------------------------------------------

	plt.savefig('{0}{1}.png'.format(target_name,SUFFIX))
	# closing the plot is necessary to relase the memory
	#  (this is a pylab issue)
	plt.close()

# ---------------------------------------------------------------------------------------------                            
# parse vot files

# catalogue 1 - Hodge source locations
cat = parse('{0}/{1}'.format(BASEDIR,CATALOGUE1))
table = cat.get_first_table()
r_hodge = table.array['RA']
d_hodge = table.array['DEC']

# catalogue 2 - JVLA source locations
cat = parse('{0}/{1}'.format(BASEDIR,CATALOGUE2))
table = cat.get_first_table()
r_vla = table.array['RA']
d_vla = table.array['DEC']

# catalogue 3 - SDSS source circles
cat = parse('{0}/{1}'.format(BASEDIR,CATALOGUE3))
table = cat.get_first_table()
r_sdss = table.array['RA']
d_sdss = table.array['DEC']
index_sdss = table.array['INDEX']

# catalogue 4 - Quasar source circles
cat = parse('{0}/{1}'.format(BASEDIR,CATALOGUE4))
table = cat.get_first_table()
r_quasar = table.array['RA']
d_quasar = table.array['DEC']
index_quasar = table.array['INDEX']

# ---------------------------------------------------------------------------------------------
# set image download parameters
HODGE_URL = 'http://third.ucllnl.org/cgi-bin/stripe82image'
FIRST_URL = 'http://third.ucllnl.org/cgi-bin/firstimage'

query_dict={}
# quirk with this is that the Ra and DEC request are in the RA field
query_dict['DEC'] = ''
query_dict['Equinox'] = 'J2000'
query_dict['ImageSize'] = CUTOUTSIZEVLA
query_dict['MaxInt'] = 3
query_dict['Survey'] = 'stripe82'
query_dict['FITS'] = 1

# ---------------------------------------------------------------------------------------------

# get array of targets
target_array = np.loadtxt(posfile,dtype='string')

if NPROCESSES > 1:
	# divy out the list to multiple processes and make cutouts
	print 'Using multiprocessing: {0} cores'.format(NPROCESSES,)
	from multiprocessing import Pool
	pool = Pool(processes=NPROCESSES)
	pool.map(make_cutout, target_array)
else:
	# or if you don't want to use multiprocessing, just loop through the list and make cutouts 
	print 'Not using multiprocessing'
	for line in target_array:
		make_cutout(line)


