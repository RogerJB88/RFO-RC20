    
'''
fits_WCS_updater.py

By: Roger Boulanger
Date:  04/29/2025   xxxx03/16/2025   xxxx03/14/2025  xxxx02/24.2025     xxxxx02/04/2025

        ***** THIS PROGRAM IS EXPERIMENTAL CODE NOT INTENDED FOR PRODUCTION! *****

This Windows program will attempt to plate solve fits image files and, if successful, update the fits header with the
WCS conversion data. The updated image file is then moved to a subfolder defined by wcs_dest variable. Unsuccessful plate
solves will result in the fits image being moved to a subfolder defined by the nowcs_dest variable.

In addition, NINA generated images have an RFO Sequence Number applied to the filename. At program
startup the last sequence  used is read from the sequence number file and then is incremented as NINA files are processed. The
last sequence number is written to the sequence number file once all liles have been processed.

And finally, LIGHT files are calibrated using Flat calibration files located in the CalFlats folder. Successfully
calibrated files are stored in the wcs_dest folder but identified by "MNc ...." vs noncalibrated as "MN ...."

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
from datetime import datetime
from astropy.io import fits
import numpy as np
from numpy import float32
import logging


def calibrateImage (inpath, flatpath):
    try:
        hdul = fits.open (inpath, 'update')
        # Get image file header to determine type of image (ie, LIGHT, DARK, FLAT, etc)
        hdr = hdul[0].header
        imageType = hdr['IMAGETYP'].upper()

        if not (imageType.startswith ('LIGHT') ):
            logger.info('NOT a LIGHT')
            hdul.close()
            return 1
        else:
            imgFilter = hdr['FILTER']
            binning = int(hdr['XBINNING'])
            gotFlat = False
            # Try to find an appropriate Flat file...
            try:
                for entry in os.scandir(flatpath):
                    logger.info(entry.name)
                    if not entry.name.startswith ('.') and entry.name.endswith('fits'):
                        with fits.open (flatpath + '/' + entry.name) as f :
                            flt = f[0].header
                            calType = flt['IMAGETYP'].upper()


                            if (calType.startswith('FLAT') and flt['FILTER'] == imgFilter and int(flt['XBINNING']) == binning):
                            #if (calType.startswith('FLAT') and int(flt['XBINNING']) == binning):
                               
                                flat_data = float32(f[0].data) / 65535.0
                                gotFlat = True
                                logger.info('Flat file found')
                                break
                            else:
                                logger.info('looping')
                            
            except Exception as e :
                logger.error("ERROR " + str(e) +'\n')
                logger.error( "Can't load Flat for "+ inpath)
                hdul.close()
                return 2
            if not (gotFlat) :
                logger.error("No matching Flat for " + inpath)
                hdul.close()
                return 3
                
            # Normalize the integer image data resulting in values 0.0000 to 1.0000
            image_data = float32(hdul[0].data) / 65535.0
            medImg = np.mean(image_data)
            logger.info('IMG MEAN= '+str(medImg))

            x, y = (image_data.shape)
            xf, yf = (flat_data.shape)
                        
            logger.info ('yf = ' + str(yf))
            # if flat file came from SkyX, trim it to 1044 x 1560
            if (yf > y) :
                flat_data = flat_data[0:1044, 1:1561]

            # Sythetic bias is derived from the camera adu offset
            bias = float32( hdr['OFFSET'] * 10 + 2.0) / 65535.0
            for x in flat_data:
                x -= bias            # subtract bias from flat

            # Synthetic dark is derived from camera adu offset + noise floor estimate at diff exposures
            exp = hdr['EXPTIME']            # image exposure time
            if (exp <= 120): exp = 0         # no factor at 120s or less exposure
            expFactor = 2 + exp * 0.0040

            darkBias = (hdr['OFFSET'] * 10 + expFactor)/65535.0       # 10 is camera adu mult factor
            for x in image_data:
                x -= darkBias             # subtract  dark from image
            
            medFlat = np.mean(flat_data)
            #logger.info('medImg= '+str(medImg))

            cal_data = float32(image_data / flat_data)
            #logger.info (cal_data)
            maxCal = np.max(cal_data)
            maxImg = np.max(image_data)
            medCal = np.mean(cal_data)

            #normalize the data to match original image mean adu
            #norm_data = cal_data  * medImg / medCal
            norm_data = cal_data  * maxImg / maxCal
            #print(final_data)
            hdul[0].data = np.int16((norm_data * 65535.0) - 32768)
            hdr['CALSTAT'] = 'BDF'
            hdr['BZERO'] = 32768
            hdr['BSCALE'] = 1
            hdul.close()
        return 0
        
    except Exception as e :
        logger.error("ERROR " + str(e) +'\n')
        logger.error("Failed calibration " + inpath)
        hdul.close()
        return 4


def gen_seqNbr(img, seqNbr):
    try:
        img0, ext = img.rsplit('.',1)
        imgmain = img0.rsplit (' ',1)
        fn = imgmain[0]
        newname = fn + ' ' + str(seqNbr).zfill(8)+'.'+ext
        return newname
    except  Exception as e:
        logger.error("ERROR Numbering" + str(e) +'\n')
        return img
        
def add_WCS_Coordinates (ip_dir, wcs_dest, nowcs_dest, flatpath):
    try:
        with open (ip_dir+'\\SEQ_NBR\\nina_seqNbr.txt', 'r') as fs:
            seqNbr = int (fs.readline() )
    except:
        seqNbr = 0       
    try:
        for entry in os.scandir (ip_dir):
            imgname = entry.name
            inpath = str(ip_dir) + '\\'+imgname
            #print('0: '+inpath)
            if not imgname.startswith ('.') and imgname.endswith('fits'):
                
                logger.info ("INITIAL INPATH: "+inpath)
                try:
                    hdul = fits.open (inpath)
                    # Get image file header to determine type of image (ie, LIGHT, DARK, FLAT, etc)
                    hdr = hdul[0].header    
                    imageType = hdr['IMAGETYP'].upper()
                    #print('1: '+imageType)
                    try:
                        swcreate = hdr['SWCREATE']
                    except :
                        swcreate = 'xxxx'
                    hdul.close()
                    if swcreate.startswith ('N.I.N.A.'):
                        seqNbr += 1
                        imagename = gen_seqNbr (imgname, seqNbr)
                        newpath = str(ip_dir)+'\\'+ imagename
                        #print('NEWPATH= '+newpath)
                        try:
                            os.rename (inpath, newpath)
                            inpath = newpath
                        except:
                            logger.error("RENAME ERROR " + imgname)
                            
                        '''
                        #print('B4:  CalibrateImage ')
                        # perform image calibration...
                        res = calibrateImage (inpath, flatpath)
                        logger.info('Calibrate Result= ' +str(res))
                        
                        #print('RES= '+str(res))
                        
                        if (res == 0):
                            calImagename = imagename.replace ('MN ', 'MNc ')
                            newpath = ip_dir + '\\' + calImagename
                            try:
                                os.rename (inpath, newpath)
                                inpath = newpath
                                logger.info('RENAMED= '+inpath)
                            except:
                                logger.error("RENAME ERROR " + newpath)
                        '''
                    
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
                        
                        calinpath = inpath.replace('MN ', 'MNc ')
                        print('calinpath= '+calinpath)
                        shutil.copy2(inpath, calinpath)
                    
                        if code[2] == '0)':
                            logger.info('PLATESOLVE SUCCESS- '+ inpath)
                            try:
                                # perform image calibration...
                                res = calibrateImage (calinpath, flatpath)
                                logger.info('Calibrate Result= ' +str(res))
                            except:
                                res = 1
                                logger.error("Calibration Error, File: "+ calinpath)
                            print("CalERROR = "+str(res))
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
                            logger.error('PLATESOLVE FAILED!!!! for ' + inpath)
                            shutil.move (inpath, nowcs_dest)
                    else:
                        try:
                            shutil.move (inpath, nowcs_dest)
                        except:
                            logger.info(inpath + ' alreary exists')
                            
                except Exception as e:
                    logger.error("ERROR 000 " + inpath + " " + str(e) +'\n')
                    shutil.move (inpath, nowcs_dest)
            else:
                if os.path.isfile(inpath):
                    os.remove(inpath)

        with open (ip_dir+'\\SEQ_NBR\\nina_seqNbr.txt', 'w') as fs:
            fs.write(str(seqNbr))
            logger.info('FINAL SeqNbr = '+str(seqNbr))
        
        
    except Exception as e:
        logger.error ("MAIN BODY ERROR *** " + str(e) +'\n')
        with open (ip_dir+'\\SEQ_NBR\\nina_seqNbr.txt', 'w') as fs:
            fs.write(str(seqNbr))

#==========================================================================================

logger = logging.getLogger('NINA_OROCESSING')   #(__name__)
logging.basicConfig(filename="C:\\PLT-SLV\\WCS_SeqNbr.log",format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %H:%M:%S', level=logging.DEBUG)

logger.info("\nSTARTING\n")

try:
    if (len(sys.argv) != 2) :
        logger.error("\n**** Add Path for the input files to be processed! ****\n")
        sys.exit (1)
    
    today = datetime.now()
    foldername = today.strftime("%Y-%m-%d")

 #   wcs_dest = "R:\\Eagle\\SkyX\\images\\"+foldername+'\\'
    wcs_dest = sys.argv[1]+'\\WCS\\'                 # **** For Test Purposes!

 #   nowcs_dest = "R:\\Eagle\\SkyX\\images\\"+foldername+'\\'
    nowcs_dest = sys.argv[1]+'\\ERR\\'               # **** For Test Purposes!
    if not os.path.isdir(wcs_dest) :
         pathlib.Path(wcs_dest).mkdir(parents=True, exist_ok=True)
         
    if not os.path.isdir(nowcs_dest) :
         pathlib.Path(nowsc_dest).mkdir(parents=True, exist_ok=True)
    flatpath = sys.argv[1] + '\\' + 'CalFlats'
    #print('FlatDir= '+flatpath)
    add_WCS_Coordinates (sys.argv[1], wcs_dest, nowcs_dest, flatpath)
    
except Exception as e:
        logger.error("ERROR *** " + str(e) +'\n')
        

