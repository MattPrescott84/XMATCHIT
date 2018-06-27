# --------------------------------------------------------------------------------------------- 

# This is a script that generates a smaller .fits cutout image from a larger mosciac file.
# Change the .fits files to the name of your choice and give it a txt file with a list of coordiniates of your sources.
# This makes small .fits cutouts that can be used to produce overlays using MAKEOVERLAYS.py

# Requires pyfits and worldpos 

# --------------------------------------------------------------------------------------------- 

import pyfits
import os
import copy
import numpy as np 
import pywcs
import worldpos # Contains definition of RA,Dec to x,y conversion


fits_file = 'YOURFITS.fits'
out_file= 'cutout1.fits'


fits_data = pyfits.getdata(fits_file,0)
header = pyfits.getheader(fits_file)
image_data = fits_data[0][0]

ID,RA,Dec = np.loadtxt('/INSERT/YOUR/DIRECTORY/WITH/coordinates.txt', unpack=True, usecols=[0,1,2], dtype=[('ID','|S30'),('RA','f8'),('dec','f8')])
# Convert RA and Dec to x and y

wcs = worldpos.wcs(fits_file, ext=0, rot='CROTA2', cd=False) # crota2 assumed to be zero (i.e. image has no rotation wrt sky)
x,y = wcs.sky2xy(RA, Dec) # POS_X, POS_Y

# Pixel numbers
pixel_x = np.rint(x).astype('i') # Rounding to the nearest integer, then ensuring type is changed from float to integer
pixel_y = np.rint(y).astype('i')

def cutout(file,xc,yc,name,xw=100,yw=100,units='pixels',JPEGFolder = "WESTTEST",clobber=True):
    """
    Inputs:
        file  - .fits filename or pyfits HDUList (must be 2D)
        xc,yc - x and y coordinates in the fits files' coordinate system
        xw,yw - x and y width
        units - specify units to use: either pixels or wcs
        outfile - optional output file
    """

    if isinstance(file,str):
        file = pyfits.open(file)
        opened=True
    elif isinstance(file,pyfits.HDUList):
        opened=False
    else:
        raise Exception("cutout: Input file is wrong type (string or HDUList are acceptable).")

    head = file[0].header
    if head['NAXIS'] > 2:
        raise Exception("Too many (%i) dimensions!" % head['NAXIS'])
    try:
        cd1 = head['CDELT1']
        cd2 = head['CDELT2']
    except:
        try:
            cd1 = head['CD1_1']
            cd2 = head['CD2_2']
        except:
            raise Exception("No CD or CDELT keywords in header")

    lonarr = ((np.arange(head['NAXIS1'])-head['CRPIX1'])*cd1 + head['CRVAL1'] )
    latarr = ((np.arange(head['NAXIS2'])-head['CRPIX2'])*cd2 + head['CRVAL2'] )

#wcs = pywcs.WCS(head)

    #xx = np.argmin(np.abs(xc-lonarr))
    #yy = np.argmin(np.abs(yc-latarr))
    #xx,yy = wcs.wcs_sky2pix(xc,yc,0)


    if os.path.exists(JPEGFolder) == False:
            os.makedirs(JPEGFolder)
    
    
    for n,x,y in zip(name,xc,yc):
        naxis1=head['NAXIS1']
        naxis2=head['NAXIS2']       
        crpix1=head['CRPIX1']
        crpix2=head['CRPIX2']
 
        if units=='pixels':
            xmin,xmax = np.max([0,x-xw]),np.min([naxis1,x+xw])
            ymin,ymax = np.max([0,y-yw]),np.min([naxis2,y+yw])
        elif units=='wcs':
            xmin,xmax = np.max([0,x-xw/np.abs(cd1)]),np.min([naxis1,x+xw/np.abs(cd1)])
            ymin,ymax = np.max([0,y-yw/np.abs(cd2)]),np.min([naxis2,y+yw/np.abs(cd2)])
        else:
            raise Exception("Can't use units %s." % units)

        if xmax < 0 or ymax < 0:
            raise Exception("Coordinate is outside of map.")

        img = file[0].data[ymin:ymax,xmin:xmax]

        crpix1-=xmin
        crpix2-=ymin
        naxis1=img.shape[1]
        naxis1=img.shape[0]
	
        newhead=copy.deepcopy(head)
        newhead['CRPIX1']=crpix1
        newhead['CRPIX2']=crpix2
        newhead['NAXIS1']=naxis1
        newhead['NAXIS2']=naxis2
        newfile = pyfits.PrimaryHDU(data=img,header=newhead)
        outfile= JPEGFolder+os.path.sep+n+'.fits'
        if isinstance(outfile,str):
            newfile.writeto(outfile,clobber=clobber)

       
       
    #if opened:
     #   file.close()

cutout1 = cutout(file=fits_file,xc=pixel_x,yc=pixel_y,xw=1000,yw=1000,units='pixels',name=ID,clobber=True)

