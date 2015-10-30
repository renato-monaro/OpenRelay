
if HAVE_DOXYGEN
directory = $(top_srcdir)/doc/man/man3/

dist_man_MANS = $(directory)
$(directory): doxyfile.stamp

doxyfile.stamp:
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

html_DATA = html/* 

# Truque para o make não reclamar quando não achar os alvos com nomes
# dos arquivos gerados pelo Doxygen
$(html_DATA):

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp
clean-local:
	rm -rf $(top_srcdir)/doc/man $(top_srcdir)/doc/html
endif
