
"launchWCS.bat" is a single line Windows bat file to run the fits_WCS_updater.py script to 
re-number images in a continuous sequence, plate-solve and update the header of all LIGHT image 
files to WCS coordinates. Successfully platesolved files are moved to the "WCS subfolder while 
image files not successfully converted are moved to the "ERR" subfolder.

In addition, a copy of each LIGHT file is made and is Calibrated automatically using the 
appropriate FLAT file shot through the matching filter). This file is then renamed with a 
leading "MNc ..." vs "MN ..." to designate it as a calibrated file.

Here is a sample script to kick off the WCS image conversion from a NINA Advanced Sequence.
At a minimum it should be inserted towards the end of the sequence (i.e., in the Shutdown 
section) as a single line (with quotes):

"C:\PLT-SLV\launchWCS.bat"   "C:\PLT-SLV\TEMP_IMAGES" 

The first group needs to point to the PLT-SLV folder (here it's at the root of drive C:\)
and the second group (the argument) needs to contain the path and folder name of the configured 
"Image file path" in the NINA Options->Imaging configuration window.

Preferably multiple identical scripts should also be placed so that they are periodically 
executed during a long imaging session.
 
# RFO-RC20
