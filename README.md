OpenFCI by Simen Kvaal
======================

**Quick links:** [documentation][docs], [releases][rels].

This is a fork of [OpenFCI][origin] to fix various bugs.

Dependencies
------------

  - C, C++, and Fortran compilers.

  - [ARPACK][arpack]: can usually be installed via the package manager.  If
    you want to install it manually, consider downloading this
    [version][arpackng] instead.

  - [LAPACK][lapack]: widely available and can be installed via the package
    manager.

  - [LAPACK Plus Plus][lpp]: this is a header-only library.  Simply
    [download lpp.0.2.b][lppdl], extract, and run:

        make HEADERDEST=/usr/local/include

    Substitute `/usr/local` with wherever you want it to be installed.  If you
    choose to install it in an unusual location such as `/home/alice`, make
    sure the `include` directory is present in the `CPLUS_INCLUDE_PATH`
    variable.  If it isn't, run:

        CPLUS_INCLUDE_PATH=/home/alice/include:$CPLUS_INCLUDE_PATH
        export CPLUS_INCLUDE_PATH

Installation
------------

Download a [pre-packaged release][rels] and then run:

    ./configure --prefix=/usr/local
    make
    make install

Substitute `/usr/local` with wherever you want it to be installed.

If you want to build directly from the *source tree* rather than download one
of the pre-packaged releases, you will need to first generate `./configure`
using [autotools][autotools]:

    autoreconf -i

[autotools]: https://en.wikipedia.org/wiki/GNU_build_system
[arpack]:    http://www.caam.rice.edu/software/ARPACK
[arpackng]:  https://github.com/opencollab/arpack-ng
[docs]:      http://folk.uio.no/simenkva/openfci/html/index.html
[lapack]:    http://netlib.org/lapack
[lpp]:       http://sourceforge.net/projects/lpp
[lppdl]:     http://sourceforge.net/projects/lpp/files/lpp/lpp-0.2.b/lpp.0.2.b.tgz
[origin]:    http://folk.uio.no/simenkva/openfci.shtml
[rels]:      https://github.com/xrf/simen-openfci/releases
