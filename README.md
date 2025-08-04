This program is intended to be used alongside the NINA astro imaging application.

The NINA Advanced Sequencer supports user scripts to invoke external programs but it is limited to Windows 
bat scripts. "launchWCS.bat" is a very simple Windows bat file who's sole function is to launch the 
fits_WCS_updater.py Python script. 

This re-numbers images in a continuous sequence, plate-solves and updates the header of all LIGHT image 
files to WCS coordinates. Successfully platesolved files are moved to the "wcs_dest" subfolder while 
image files not successfully converted are moved to the "nowcs_dest" subfolder.

In addition, a copy of each LIGHT file is made and is Calibrated automatically using master BIAS, DARK, and 
appropriate FLAT files shot through the matching filter). This LIGHT file header is updated to Calstate = 'BFD'
and the file renamed with a leading "MNc ..." vs the normal "MN ..." to designate it as a calibrated file.

Here is a sample script to kick off the WCS image conversion from a NINA Advanced Sequence. This script is 
placed at the end of each target template so that images are processed as each group is complete.

Here is an example of the NINA script to launch the processing functions:

"C:\PLT-SLV\.\launchWCS.bat C:\PLT-SLV\TEMP_IMAGES" (without quotes)

The first group needs to point to the PLT-SLV folder (here it's at the root of drive C:\) and the second 
group (the argument) needs to contain the path and folder name of the configured "Image file path" in the 
NINA Options->Imaging configuration window. This is a temporary holding buffer and this example shows it
on the local drive but it should be located on a remote harddrive or on a ram disk.

# RFO-RC20
