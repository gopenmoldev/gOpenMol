#
# Instructions to install the program
#
#  Eero Häkkinen 2005
#

To install gOpenMol please do the following:

(1) Make sure you have these libraries and utilities installed:

      1) Tcl/Tk
           1) Tcl programs and libraries (version 8.4 or newer).
                http://www.tcl.tk/
           2) Tk libraries (the same version as Tcl).
                http://www.tcl.tk/
           3) BWidget
                http://sourceforge.net/projects/tcllib/
           4) Tcllib
                http://sourceforge.net/projects/tcllib/
           5) TclXML
                http://sourceforge.net/projects/tclxml/

      2) OpenGL
           1) OpenGL library (libGL).
           2) OpenGL Utility Library (libGLU).
           3) OpenGL Utility Library Toolkit (libglut).

           MESA 3D Graphics Libraries contains these.
             http://www.mesa3d.org
           Both MesaLib and MesaDemos packages are needed
           (the latter contains glut).

           These are standard OpenGL libraries and should come
           with any OpenGL implementation.

      3) The Independent JPEG Group's JPEG runtime library
           http://www.ijg.org/

      4) Python (optional)
           http://www.python.org

      If you install from source code, you have to configure, build and
      install these libraries and utilities.

      If you install precompiled binaries (binary tar packages or
      RPM or DEB packages or similar), make sure you install
      development files (headers and library links), too.
      They are usually in *-dev packages.

(2) Copy the content of the lib/ directory in the Tcl/Tk 
    distribution to the gOpenMol lib/ directory.

(3) Run configure. To get list of possible options run `./configure --help'.
    If you installed libraries to non-standard locations (i.e. using
    --prefix="$HOME/gom"), you probably have to at least
    give --with-dirs option or set CPPFLAGS and LDFLAGS variables.
    If you want to keep object files separated from source files, you can
    build gOpenMol in another directory.
    Examples:
        ./configure

        ../gOpenMol/src/configure \
                --prefix="$HOME/gom" \
                --without-python \
                --with-dirs=$HOME/Mesa:$HOME/tcl

(4) Build gOpenMol, plugins and utility programs using make.
    Command (see the end of config.mk for other possible targets):
        make

(5) Install gOpenMol.
      Command (only one of the followings):
          make altinstall     # Replace files in bin and plugins directories.
                              # --prefix=PREFIX has no effect.
          make install        # Install all files.
          make install-local  # Install only platform dependent files.

(6) If you did not use --with-dirs and some libraries are in non-standard
    locations, you may have to copy or link shared library files
    (libtclX.X.so.*, libtkX.X.so.*, libGL.so.*, libGLU.so.*, libglut.so.*,
    libjpeg.so.* and libpythonX.X.so.*) to the lib/ directory.
