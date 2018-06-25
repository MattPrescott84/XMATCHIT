
#============================================
# Script for eyeballing overlays

#To run type:
#  python XMATCHIT.py <<dir with overlays>>
#or
#  python XMATCHIT.py -l <<list of object names (J...)>> <<dir with overlays>>
#Script assumes overlay images files are named J.....png and J....ZOOMED.png

#To recover from a crash use:
#  python XMATCHIT.py -r
#This will make an eyeball_out.dat file from the eyeball_log.log file
#This is needed since eyeball_out.dat isn't written on the fly in the code
#and eyeball_log.log isn't always the same as _out.dat

#If you stop and start again can annotate onto the end of the 'eyeball_out.dat' data file
#rather than starting again. Default is to append onto end of file.

#Added back option- skips to previous object- in case you mess up
#Added option to recover from crash


#============================================
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage
import sys
import glob
import os

#============================================
print '###################################################'
print ' Hi, and welcome to XMATCHIT\n version 2\n!'
print '###################################################\n'
#Command line stuff
if len(sys.argv)<=1:
    print 'to run type: python XMATCHIT.py <<dir with overlays>>\n or python  XMATCHIT.py -l <<list of objects>> <<dir with overlays>>'
    print ' to recover from a crash use\n python  XMATCHIT -r \n witch will create an eyeball_out.dat file from eyeball_log.log ready'
    print ' for running  XMATCHIT.py again, dont just run  XMATCHIT.py again after a crash since it will overwrite the log file,'
    print ' at least save the log under a different name'
    sys.exit()
if sys.argv[1]=='-l':
    file_dir = sys.argv[3]
    objlist = np.loadtxt(sys.argv[2],usecols=[0],dtype=str)
elif sys.argv[1]=='-r':
    print 'recovering from eyeball_log.log to eyeball_out.dat..........'
    logfile = open('eyeball_log.log','r')
    logdat = [[line.split()[0],line] for line in logfile]
    logfile.close()
    outnames,outdat = [],[]
    for ldi in logdat[::-1]:
        if ldi[0] not in outnames:
            outdat.append(ldi[1])
            outnames.append(ldi[0])
    outfile = open('eyeball_out.dat','w')
    outfile.write('#name match FR_class match_ID flags user comment\n')
    for odi in outdat[::-1]: outfile.write(odi)
    outfile.close()
    print 'done........ hopefully'
    sys.exit()
else:
    file_dir = sys.argv[1]
    objlist = [pathi.split('/')[-1][:20] for pathi in glob.glob(file_dir+'/*ZOOMED.png')]

############################################################################
#                              FUNCTIONS
#-----------------------------------------
#Silly formating func
def resize(instr,mlen):
    if len(instr)>mlen: print 'string too long!!', instr
    return instr + (mlen - len(instr))*' '

#-------------------------------------------
#Function to read and display images
def readAndDisp(fdir,oname):
    gs1 = gridspec.GridSpec(1,2)
    gs1.update(wspace=0.,hspace=0.,left=0.,right=1.,bottom=0.,top=1.)
    plt.subplot(gs1[0,0])
    plt.imshow(plt.imread(fdir+'/'+oname+'ZOOMED.png'))
    plt.xticks([],[])
    plt.yticks([],[])

    plt.subplot(gs1[0,1])
    plt.imshow(plt.imread(fdir+'/'+oname+'.png'))
    plt.xticks([],[])
    plt.yticks([],[])
#--------------------------------------------
#Way to get input that must be certain values
def question(q_text,possible_vals_list,helptxt='not helping'):
    answer = 'wtf'
    while answer == 'wtf':
        answer = raw_input(q_text)
        if not answer in possible_vals_list:
            print '\nNot a valid input!\n'
            answer = 'wtf'
        elif answer == '?':
            print helptxt
            answer = 'wtf'
    return answer

####################################################################
#Get initials of user
mlen_in = 5
user_initials = raw_input(' Please enter your initials:\n')
if (len(user_initials)>mlen_in):
    print ' you have too many initials and are too posh to play '
    sys.exit()
user_initials = resize(user_initials,mlen_in)

#-----------------------------------------------------
#Does the output file already exist,
#do you want to just add onto the end of this...
#Default is to append onto end of file
prevFileAp='n'
done_objects = []
if os.path.isfile('./eyeball_out.dat'):
    prevFileAp = raw_input(' Found previous results file eyeball_out.dat\n annotate from the end of this file (y/n)\n n will destroy previous file, q to quit now\n')
    if prevFileAp=='q': sys.exit()
logfile = open('eyeball_log.log','w')
if prevFileAp=='n':
    ofile = open('eyeball_out.dat','w')
    ofile.write('#name match FR_class match_ID flags user comment\n')
else:
    ofile = open('eyeball_out.dat','r+')
    for line in ofile:
        done_objects.append(line.split()[0])
        logfile.write(line)

#--------------------------------------------
plt.ion()
plt.figure(figsize=[18,8]) # ********************* LR ******************
nobj    = len(objlist)
out_arr = []
mlen_id = 12
mlen_m  = 2
mlen_Fid= 3
mlen_flg= 1

#--------------------------------------------
#Loop over files, display images, ask for user input
#Love that fortran 66 loop
i=-1
while i < nobj-1:
    i += 1
    objname = objlist[i]
    print i, objname
    if objname in done_objects: continue

    readAndDisp(file_dir,objname)
    plt.show()

#Now human input stuff
    match    = 'wtf'
    match_id = '0'
    FR_class = '9'
    flags    = 0
    comment  = '\'\''

    match = question(' Is this a match? If so red photo z or cyan spectro z \n g/s/n or ? for more options\n',['g','s','g*','s*','n','n?','g?','s?','back','quit','?'],\
                         helptxt='  options are:\n g  = red galaxy\n s  = cyan specZ\n n  = no match\n add a \'*\' for a near match where radio centroid isnt quite on optical\n add a \'?\' to flag sources to be looked at again\n \'back\' to go to previous object\n \'quit\' to exit and save progress\n')

    if match =='back':
        out_arr = out_arr[:-1]
        i = i-2

    if match=='quit': break

    if len(match)==2:
        if match[1]=='?': flags += 1
        if match[1]=='*': flags += 2
        match = match[0]

    if match in ['n','g','s']:
        if match !='n':
#Faf to make sure match_id is an integer, someone show me how to do better.
            loop = True
            while loop:
                match_id = raw_input(' What is the match ID number?\n')
                try:
                    crap = int(match_id)
                    loop = False
                except:
                    print '\nnot an integer!\n'
        else:
            if question(' Part of a larger radio structure (y/n)\n',['y','n']) == 'y': flags += 4
        FR_class = question(' FR I/II/compact/unclear? (1/2/0/5 also 1?/2?/1+2/5 or ? for help)\n',\
                                ['1','2','0','1?','2?','1+2','5','?'],\
                                helptxt=' 1/2 for FR 1s and 2s\n 0 for compact\n 5 for extended but unclear\n can also have 1?/2? for probable identifications\n and 1+2 for something that clearly has attributed of both\n')

        match    = resize(match,mlen_m)
        match_id = resize(match_id,mlen_id)
        FR_class = resize(FR_class,mlen_Fid)
        flags    = resize(str(flags),mlen_flg)
        comment  = '\'' + raw_input(' Comments (not required!)\n') +'\''
        outline  = '  '.join([objname,match,FR_class,match_id,flags,user_initials,comment,'\n'])
        out_arr.append(outline)
        logfile.write(outline)


    plt.clf()


for line in out_arr: ofile.write(line)
ofile.close()
logfile.close()


