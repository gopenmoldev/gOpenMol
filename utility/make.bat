REM 
REM make executables of the utility programs
REM Leif Laaksonen/2001
REM Eero Häkkinen/2003
REM

REM make in utility

f90  ambera2b.f
copy ambera2b.exe ..\bin
f90  autodock2plt.f
copy autodock2plt.exe ..\bin
f90  charmmtrj.f
copy charmmtrj.exe ..\bin
cl   /MD /O2 contman.c
copy contman.exe ..\bin
f90  convert.f
copy convert.exe ..\bin
f90  delphi2plt.f
copy delphi2plt.exe ..\bin
cl   /MD /O2 g94cub2pl.c
copy g94cub2pl.exe ..\bin
cl   /MD /O2 gamess2plt.c
copy gamess2plt.exe ..\bin
cl   /MD /O2 gcube2plt.c
copy gcube2plt.exe ..\bin
f90  ins2plt.f
copy ins2plt.exe ..\bin
cl   /MD /O2 jaguar2plt.c
copy jaguar2plt.exe ..\bin
f90  kont2plt.f
copy kont2plt.exe ..\bin
cl   /MD /O2 mpcube2plt.cpp
copy mpcube2plt.exe ..\bin
cl   /MD /O2 pltfile.c
copy pltfile.exe ..\bin
f90  sybyl2amber.f
copy sybyl2amber.exe ..\bin
f90  trajmerge.f
copy trajmerge.exe ..\bin
f90  xmol2bamber.f
copy xmol2bamber.exe ..\bin
f90  xplor2charmm.f
copy xplor2charmm.exe ..\bin

REM make in utility/turbotools

cd turbotools
cl   /MD /O2 tmole2plt.c
copy tmole2plt.exe ..\..\bin
cd ..

REM make in utility/UHBD

cd UHBD
f90  gridasc2plt.f
copy gridasc2plt.exe ..\..\bin
f90  gridbin2plt.f
copy gridbin2plt.exe ..\..\bin
cd ..

cd ..

REM make in density

cd density
cl   /MD /O2 density.c
copy density.exe ..\bin
cd ..

REM make in icon8

cd icon8
f90  icon8.f rsp.f output.f
copy icon8.exe ..\bin
cd ..

REM make in probsurf

cd probsurf
cl   /MD /O2 probsurf2.0.c
copy probsurf2.0.exe ..\bin\probsurf.exe
cd ..

REM make in vss

cd vss
cl   /O2 vssmod.c
copy vssmod.exe ..\bin\vss.exe
cd ..

cd utility
pause
