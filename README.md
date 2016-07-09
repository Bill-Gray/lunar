# lunar
Basic astronomical functions for solar system ephemerides,  time systems,
coordinate systems,  etc.  This includes some utilities based on these
functions,  such as a calendar computer and a utility to numerically
integrate asteroid orbits.  The code can be built for Windows,  Linux,
or BSD.  Some documentation (in need of updates) is at

http://www.projectpluto.com/source.htm#astrocalc

On Linux,  run 'make' to build the library and various test executables.
You can then run 'make install' to put libraries in /usr/local/lib and some
include files in /usr/local/include.  (You will probably have to make that
'sudo make install'.)  For BSD,  and probably OS/X,  run 'gmake CLANG=Y'
(GNU make,  with the clang compiler),  then 'sudo gmake install'.
