# XMATCHIT
This repository contains code to make overlays for radio/optical surveys, and code to cross-match and classify radio sources with their optical counterparts.   

 
# GETCUTOUT.py  
 This is a script that generates a small .fits cutout images from a larger radio mosiac file.
 Change the .fits files to the name of your choice and give it a .dat file with a list of coordiniates of your radio sources.
 These small cutouts can then be used in MAKEOVERLAYS.py to produce overlays. 

 This requires pyfits and worldpos.py 

# MAKEOVERLAY.py

This script can be used to make overlays of radio sources in .png format, like the image below.  

In this example radio contours from two radio surveys: The VLA Stripe 82 Snapshot Survey, Heywood et al. 2016, http://adsabs.harvard.edu/abs/2016MNRAS.460.4433H (green contour) and Hodge et al. 2011 (blue contours) are overlayed on top of greyscale SDSS images. 

![alt text](https://github.com/MattPrescottAstro/XMATCHIT/blob/master/EXAMPLE.jpg)

This requires matplotlib and the APLpy package of functions which can be found here: https://aplpy.github.io

# XMATCHIT.py

 Script for the visual matching and classification of radio sources using overlays produced by MAKEOVERLAY.py.
 In this example it needs two sets of overlays of different size (eg. 0.5 x 0.5 arcminutes and 3.0 x 3.0 arcminutes) placed in a directory.  
 The script prompts the user to classify whether or not radio source components have an optical counterpart, and then prompts the user to classify the radio source morphology (compact, extended, FRI/II etc.) 
 
 To run type:
 python XMATCHIT.py /dir/containing/overlays
 
This produces a file called eyeball_out.dat which contains a list of source IDs and their classsifications.

