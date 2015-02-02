Install gOpenMol
================

Steps to do to install gOpenMol on UNIX

I'm only distributing gOpenMol (currently) as a binary
distribution, if you need the source code you have to get in
touch with me. The program is distributed free of charge for
academic installations and academic use.

Most likely you have just got a compressed tar file
version of gOpenMol (gopenmol.tar.Z). You have to
uncompress the file first with the command:


    uncompress gopenmol.XXX.tar.Z

Copy first the gopenmol.tar file to directory
where you want to install gOpenMol. The tar file will
create a directory gopenmol containg the
needed subdirectories, programs and data. You can untar
the tar file with the command:

    tar xf gopenmol.xxx.tar

Change directory down to the gopenmol directory
Run the installation script install with the command: 

    ./install

You are now ready to use gOpenMol. To run gOpenMol write:
    

     bin/rungOpenMol

====================================================

Steps to do to install gOpenMol on Windows 95/98 and NT
====================================================


gOpenMol is now distributed with source code!
The program is distributed free of charge for
academic installations and academic use.


    1.  You have got the gopenmol.zip file downloaded from the
        net. Take a program (WinZip is a good candidate) to unzip
        the file.
    2.  Choose the directory where you want to download the
        program including the files. The unzip program creates a
        directory gopenmol and places the file in the
        respective directories.
    3.  Change directory to the gopenmol directory and run the install.bat
        file. This procedure will call the wish shell and run the
        installation procedure. After the installation procedure
        you should have a new bat file in your gopenmol/bin
        directory called rungOpenMol.bat. This file starts
        gOpenMol and makes all the needed environment
        assignments.
    4.  Run the script rungOpenMol.bat always when you
        want to run the gOpenMol program.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Very often the Windows 95/98 DOS shell has, by default, a very small amount
of memory served for the environment values. This will prevent gOpenMol 
from reserving all the needed environment values through the set command. 
What you see is that the DOS shell complains about running out of 
environment space when you run the bin\rungOpenMol.bat script.

To cure this problem click on the upper left corner of your DOS window
and go to the Property option to give more space (memory) to your DOS 
shell.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Uninstall gOpenMol

Both the Unix and the Windows versions can be removed by just deleting
the root directory and the directories under the root.

The Windows version uses no register entries and places no files in the
system directory.
 
Maintained by (http://www.csc.fi/lul/leif/leif.laaksonen.html) Leif
Laaksonen (Leif.Laaksonen@csc.fi)

Last modified 2003/08/08
