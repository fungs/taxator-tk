# Compile and install

If you downloaded the source, make sure you installed all needed
software development packages. You need to have a recent version of the Boost
C++ libararies (www.boost.org) installed (>= version 1.34), including the
development header files. You will also need the build management tool
"cmake" >= version 3.1 and the UNIX tools "make" and "uname" which should be
contained in a standard installation.

On a Unix (Linux) system type

    ./build.sh

and all programs will be built in a sub-directory called "bin".

There is no installation procedure. Either copy the executables from the build
directory to your prefered directory in your PATH variable like /usr/local/bin
or execute them by prefixing them with the build directory.
