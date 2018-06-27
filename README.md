# XMATCHIT
This repository contains code to make overlays for radio/optical surveys, and code to cross-match and classify radio sources with their optical counterparts.   

 
# GETOVERLAYS.py  
 This is a script that generates a small .fits cutout images from a larger radio mosiac file.
 Change the .fits files to the name of your choice and give it a .dat file with a list of coordiniates of your radio sources.
 These small cutouts can then be used in MAKEOVERLAYS.py to produce overlays. 

 This requires pyfits and worldpos.py 

# MAKEOVERLAYS.py



This requires matplotlib and the APLpy https://aplpy.github.io


