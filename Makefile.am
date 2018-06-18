bin_PROGRAMS    = crumble
crumble_SOURCES = snp_score.c str_finder.c bed.c
pkginclude_HEADERS = bed.h str_finder.h tree.h utlist.h
crumble_LDADD   = -lm -lhts

AM_CPPFLAGS     = @HTSLIB_CPPFLAGS@
AM_LDFLAGS      = @HTSLIB_LDFLAGS@

man_MANS        = crumble.1
EXTRA_DIST      = README.md $(man_MANS)
