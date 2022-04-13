# Configure script for crumble.
#
#    Copyright (C) 2018 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

dnl Process this file with autoconf to produce a configure script
AC_INIT([Crumble], 0.9.0)
AC_PREREQ([2.63])  dnl This version introduced 4-argument AC_CHECK_HEADER

AM_INIT_AUTOMAKE
AM_INIT_AUTOMAKE([foreign])

m4_include([ax_with_htslib.m4])

AC_PROG_CC
AC_SYS_LARGEFILE

AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib development files not found.])
fi

AC_OUTPUT(Makefile)

echo
echo "Using Htslib CPPFLAGS: $HTSLIB_CPPFLAGS"
echo "Using Htslib LDFLAGS:  $HTSLIB_LDFLAGS"
echo
