   
'''
fits_WCS_Seq_updater.py

By: Roger Boulanger
Date:  06/18/2025   Alpha6.18

        ***** THIS PROGRAM IS EXPERIMENTAL CODE NOT INTENDED FOR PRODUCTION! *****

This Windows program will attempt to plate solve fits image files and, if successful, update the fits header with the
WCS conversion data. The updated image file is then moved to a subfolder defined by wcs_dest variable. Unsuccessful plate
solves will result in the fits image being moved to a subfolder defined by the nowcs_dest variable.

In addition, NINA generated images have an RFO Sequence Number applied to the filename. At program
startup the last sequence number used is read from the sequence number file and then is incremented as NINA files are processed. The
last sequence number is written to the sequence number file once all files have been processed.

And finally, LIGHT files are calibrated using master Flat calibration files located in the CalFlats folder as well as  master Dark
in the CalDark folder. The master Dark is scaled for exposure so that a single master Dark can be used to calibrate Light file
shat at different exposures. Successfully calibrated files are stored in the wcs_dest folder but identified by "MNc ...." vs 
noncalibrated as "MN ...."

Python3 must be available on the computer and Astap plate solver with an appropriate star database must be 
installed in the "Program Files/astap" folder.

The program can be launched from the command line or from a script inserted into a NINA Advanced Sequence.
The path to the folder containing the images to be processed must be provided as an argument.
'''

import os
import io
import sys
import pathlib
import shutil
import subprocess
from datetime import datetime, timedelta
from astropy.io import fits
import numpy as np
from numpy import float32
import logging




'''
The getMstrDrk function is called by add_WCS_Coordinates function and loads the master dark file and populates several variables 
and returns these so that this operation is only executed once and the returned variables are used over and over again as 
individual LIGHT files are calibrated.
'''


def getMstrDrk(darkpath, biaspath):
    try:
        with fits.open(darkpath+'\\Dark_Master-c.fits') as fdrk:
            dark_hdr = fdrk[0].header
            dark_exp = float32(dark_hdr['EXPTIME'])
            dark_binning = int(dark_hdr['XBINNING'])
            # scale the data 0.000 to 1.000
            global dark_buf
            dark_buf = (float32(fdrk[0].data))/65535.0

        with fits.open(biaspath+'\\Bias_Master-c.fits') as fbias:
            global bias_buf
            bias_buf = (float32(fbias[0].data))/65535.0
	
        return [dark_exp, dark_binning]
    except:
        logger.error(" Can't open " + darkpath +' or ' + biaspath)
        return 2

'''
 The calibrateImage function is called by the process_NINA_images function. "ip_dir' is the working directory and "calinpath" is 
 the path to the image file to be calibrated. This folder will be searched for a FLAT file with matching Filter and Binning settings. 
 The calibration formula used is: (LIGHT - DARK) / (FLAT - BIAS). A synthetic value of BIAS is derived from the ASI2600MM Specifications 
 and serves until bias files are collected.
''' 
def calibrate(ip_dir, calinpath, dark_exp, dark_binning):
    global bias_buf        
    global dark_buf
    try:
        try:
            fimg = fits.open (calinpath, 'update')
        except:
            logger.error(" Can't open " + calinpath)
            return 2
            
        hdr = fimg[0].header
        img_binning = int(hdr['XBINNING'])
        img_type = hdr['IMAGETYP'].upper()
        img_filter = hdr['FILTER']
     
        imgbuf = (float32(fimg[0].data))/65535.0
        img_exp = float32(hdr['EXPTIME'])

        # If more than 20s exposure difference..
        if (abs(img_exp - dark_exp) > 20):
            exp_scaleFactor = float32(img_exp / dark_exp)
            logger.info('expScaleFactor = '+str(exp_scaleFactor))
            hotpix_threshold = float32(2024.0/65535.0)

            dark_buf -= bias_buf
            for x in dark_buf:
                chk = (x < hotpix_threshold)
                x[chk] *= exp_scaleFactor
            dark_buf += bias_buf     
        
        diff_img= np.subtract(imgbuf, dark_buf)
        # Eliminate unrealistic numbers...
        diff_median = np.median(diff_img)
        logger.info("DIFFIMG-MEDIAN="+str(diff_median))
        y = diff_median * 0.25
        for x in diff_img:
            chk = (x <= 0)
            x[chk] = y            
 
        try:
            got_flat = False
            flatpath = ip_dir + '\\FLATS'
            for flat in os.scandir(flatpath):
               #logger.info(flat.name)
                if flat.name.startswith ('F') and flat.name.endswith('.fits'):
                    with fits.open (flatpath + '\\' + flat.name) as f :
                        flt_hdr = f[0].header
                        calType = flt_hdr['IMAGETYP'].upper()

                        if (calType.startswith('FLAT') and flt_hdr['FILTER'].startswith(img_filter) and flt_hdr['XBINNING'] == img_binning):
                            flat_data = float32(f[0].data) / 65535.0
                            got_flat = True
                            logger.info('Flat file found')
                            break
                                            
        except Exception as e :
            logger.error("ERROR " + str(e) +'\n')
            logger.error( "Can't load Flat for "+ calinpath)
            fimg.close()
            return 3

        if not (got_flat) :
            logger.error("No matching Flat for " + calinpath)
            fimg.close()
            return 4

        flat_data -= bias_buf
        maxFlat = np.max(flat_data)
        flat_data /= maxFlat
        diff_img /= flat_data
        fimg[0].data = np.int16(diff_img * 65535 - 32768)
        hdr['CALSTAT'] = 'BDF'
        hdr['BZERO'] = 32768
        hdr['BSCALE'] = 1
        fimg.close()
        return 0

    except Exception as e :
        print('ERROR= '+str(e)+'\n')
        logger.error('ERROR= '+str(e)+'\n') 
        return 1


'''
The gen_seqNbr function rplaces the NINA generated file sequence number with and RFO standard  continuously incrementing 
eight digit sequence number. This overcomes the problem of NINA's numbers being reset to zero at startup. It is called by 
the add_WCS_coordinates() function.
'''
def gen_seqNbr(img, seqNbr):
    # replace NINA seq nbr with master sequential number...
    try:
        img0, ext = img.rsplit('.',1)
        imgmain = img0.rsplit (' ',1)
        fn = imgmain[0]
        newname = fn + ' ' + str(seqNbr).zfill(8)+'.'+ext
        return newname
    except  Exception as e:
        logger.error("ERROR Numbering" + str(e) +'\n')
        return img
	    
'''              
The process_NINA_images function first calls the gen_seqNbr function to re-number all NINA generated image files. Then with 
LIGHT files only it plate-solves the image using ASTAP plate solver in order to add WCS coordinates to the image header and 
if successful, designate it as having been plate solved and written to the wcs_dest folder. An unsuccessful plate-solve 
results in the file being written to the nowcs_dest folder and no further action is taken on this file. if successful the file 
is then copied but with a leading "MNc..." to the filename and then sent to the calibrateImage function to be calibrated. If 
the "MNc...." file has been successfully calibrated its header is modified to indicate its calibrated status and then it is 
written to the wcs_dest folder. If not successful the file is deleted.
'''
def process_NINA_images(ip_dir, wcs_dest, nowcs_dest):
    try:
        with open (ip_dir+'\\SEQ_NBR\\nina_seqNbr.txt', 'r') as fs:
            seqNbr = int (fs.readline() )
    except:
        seqNbr = 0
    try:
        darkpath = ip_dir+'\\DARK'
        biaspath = ip_dir+'\\BIAS'
        [dark_exp, dark_binning] = getMstrDrk(darkpath, biaspath)
        
        
    except:
        logger.error("Cannot get: "+darkpath+'\\xxxx')
        return 5
       
    try:
        for entry in os.scandir (ip_dir):
            imgname = entry.name
            inpath = str(ip_dir) + '\\'+imgname

            if imgname.startswith ('MN ') and imgname.endswith('.fits'):
                
                #logger.info ("INITIAL INPATH: "+inpath)
                try:
                    hdul = fits.open (inpath)
                    # Get image file header to determine type of image (ie, LIGHT, DARK, FLAT, etc)
                    hdr = hdul[0].header    
                    imageType = hdr['IMAGETYP'].upper()
                    try:
                        swcreate = hdr['SWCREATE']
                    except :
                        swcreate = 'xxxx'
                    hdul.close()

                    if swcreate.startswith ('N.I.N.A.'):
                        seqNbr += 1
                        imagename = gen_seqNbr (imgname, seqNbr)
                        newpath = str(ip_dir)+'\\'+ imagename 
                        try:
                            os.rename (inpath, newpath)
                            inpath = newpath
                        except:
                            logger.error("RENAME ERROR " + imgname)                            
                    
                    # make sure this is a "LIGHT" file needing plate-solving
                    if (imageType.startswith ('LIGHT')):                        
                        try:
                            # Filenames with embedded spaces must be enclosed in " marks. The "update" option causes 
                            # the image header to be updated with the WCS information.
                            res = subprocess.run ("C:\\Program Files\\astap\\astap.exe -f "+'\"'+inpath+'\"'+" -update")
                            code = str(res).split('=')
                            logger.info('ASTAP Result= '+ code[2])
                        except Exception as e:
                            logger.error("ASTAP EXECUTION ERROR " + str(e) +'\n')                  
                    
                        if code[2] == '0)':
                            logger.info('PLATESOLVE SUCCESS- '+ inpath)
                            try:
                                # perform image calibration...
                                calinpath = inpath.replace('MN ', 'MNc ')

                                shutil.copy2(inpath, calinpath) 

                                res = calibrate (ip_dir, calinpath, dark_exp, dark_binning)
                                logger.info('Calibrate Result= ' +str(res))
                            except:
                                res = 1
                                logger.error("Calibration Error, File: "+ calinpath)

                            logger.info("CalSTATUS = "+str(res))
                  
                            if (res == 0):
                                shutil.move (calinpath, wcs_dest)
                            else:
                                os.remove(calinpath)
                                
                            shutil.move (inpath, wcs_dest)
                            # This will point us to the corresponding extra result files and delete them
                            rslt_ini = inpath.replace (".fits", ".ini")
                            os.remove (rslt_ini)
                            wcsfile = rslt_ini.replace (".ini", ".wcs")
                            os.remove (wcsfile)

                        else:
                            # ASTAP may have created this error file
                            rslt_ini = inpath.replace (".fits", ".ini")
                            if os.path.isfile(rslt_ini):
                                os.remove (rslt_ini)
                            logger.error('PLATESOLVE FAILED!!!! for ' + inpath)
                            shutil.move (inpath, nowcs_dest)
                    else:
                        try:
                            shutil.move (inpath, nowcs_dest)
                        except:
                            logger.info(inpath + ' already exists')
                            os.remove(inpath)

                except Exception as e:
                    logger.error("ERROR 000 " + inpath + " " + str(e) +'\n')
                    shutil.move (inpath, nowcs_dest)
            else:
                if os.path.isfile(inpath):
                    os.remove(inpath)

        with open (ip_dir+'\\SEQ_NBR\\nina_seqNbr.txt', 'w') as fs:
            fs.write(str(seqNbr))
            logger.info('FINAL SeqNbr = '+str(seqNbr))
        return 0
        
    except Exception as e:
        logger.error ("MAIN BODY ERROR *** " + str(e) +'\n')
        with open (ip_dir+'\\SEQ_NBR\\nina_seqNbr.txt', 'w') as fs:
            fs.write(str(seqNbr))
        return 1
#================================================================================================================
#================================================================================================================

logger = logging.getLogger('NINA_OROCESSING')   #(__name__)
logging.basicConfig(filename="C:\\PLT-SLV\\Plt-Slv_log.txt",format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %H:%M:%S', level=logging.DEBUG)

logger.info("\nSTARTING\n")
print("\nPROCESSING IMAGES\n")
try:
    if (len(sys.argv) != 2) :
        logger.error("\n**** Add Path for the input files to be processed! ****\n")
        print("\n**** Add Path for the input files to be processed! ****\n")
        sys.exit (1)
    
    today = datetime.today() - timedelta(hours=12)
    foldername = today.strftime("%Y-%m-%d")

    #wcs_dest = "R:\\Eagle\\SkyX\\images\\"+foldername+'\\'
    wcs_dest = sys.argv[1]+'\\WCS\\'                 # **** For Test Purposes!

    #nowcs_dest = "R:\\Eagle\\SkyX\\images\\"+foldername+'\\'
    nowcs_dest = sys.argv[1]+'\\ERR\\'               # **** For Test Purposes!
    if not os.path.isdir(wcs_dest) :
         pathlib.Path(wcs_dest).mkdir(parents=True, exist_ok=True)
         
    if not os.path.isdir(nowcs_dest) :
         pathlib.Path(nowsc_dest).mkdir(parents=True, exist_ok=True)

    res = process_NINA_images (sys.argv[1], wcs_dest, nowcs_dest)
    logger.info ("End Session with res = " + str(res))

except Exception as e:
    logger.error("ERROR *** " + str(e) +'\n')
        
