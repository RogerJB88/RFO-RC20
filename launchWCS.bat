@ECHO OFF
SET mypath=%~dp0
python  %mypath:~0,-1%\fits_WCS_Seq_updater.py  %1

