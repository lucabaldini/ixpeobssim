@ECHO OFF
REM See this stackoverflow question
REM http:\\stackoverflow.com\questions\3827567\how-to-get-the-path-of-the-batch-script-in-windows
REM for the magic in this command
SET WIN_SETUP_DIR=%~dp0
SET SETUP_DIR=%WIN_SETUP_DIR:\=/%

REM
REM Base package root. All the other releavant folders are relative to this
REM location.
REM
SET IXPEOBSSIM_ROOT=%SETUP_DIR:~0,-1%
ECHO "IXPEOBSSIM_ROOT set to " %IXPEOBSSIM_ROOTROOT%

REM
REM Setup the PYTHONPATH.
REM
set PYTHONPATH=%IXPEOBSSIM_ROOT%;%PYTHONPATH%
echo "PYTHONPATH set to " %PYTHONPATH%

REM
REM Add the bin folder to the PATH environmental variable.
REM 
set PATH=%IXPEOBSSIM_ROOT%\ixpeobssim\bin;%PATH%
echo "PATH set to " %PATH%
