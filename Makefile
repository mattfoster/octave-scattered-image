## Generic Makefile to allow the octave-forge packages to be build and 
## installed using "configure; make all; make install". This Makefile
## includes the capability to install to a temporary location, and install
## an on_uninstall.m file that prevents the user removing this package
## with Octave's package manager. This is useful for for the distribution's
## various package managers and is forced by defining DESTDIR and DISTPKG.

PKGDIR := $(shell pwd | sed -e 's|^.*/||')
TMPDIR ?= /tmp
PACKAGE ?= $(TMPDIR)/$(PKGDIR).tar.gz
PKG := $(shell echo $(PKGDIR) | sed -e 's|^\(.*\)-.*|\1|')

all: build package

build:
	@if [ -e src/Makefile ]; then \
	   $(MAKE) -C src all; \
	fi

package: build
	@if [ -e src/Makefile ]; then \
	  mv src/Makefile src/Makefile.disable; \
	fi; \
	if [ -e src/configure ]; then \
	  mv src/configure src/configure.disable; \
	fi; \
	cd ..; tar -czf $(PACKAGE) $(PKGDIR); \
	echo tar -czf $(PACKAGE) $(PKGDIR);

install:
	@cd ../; \
	if [ "X${DISTPKG}X" != "XX" ]; then \
	  stripcmd="unlink(pkg('local_list'));unlink(pkg('global_list'));"; \
	fi; \
	if [ "X$(DESTDIR)X" = "XX" ]; then \
	  pkgdir=`octave -H -q --no-site-file --eval "warning('off','all');pkg('install','$(PACKAGE)');l=pkg('list');disp(l{cellfun(@(x)strcmp(x.name,'$(PKG)'),l)}.dir);$$stripcmd;"`; \
	else \
	  shareprefix=$(DESTDIR)/`octave -H -q --no-site-file --eval "disp(fullfile(OCTAVE_HOME(),'share','octave'));"`; \
	  libexecprefix=$(DESTDIR)/`octave -H -q --no-site-file --eval "disp(fullfile(octave_config_info('libexecdir'),'octave'));"`; \
	  octprefix=$$shareprefix/packages; \
	  archprefix=$$libexecprefix/packages; \
	  if [ ! -e $$octprefix ]; then \
	    mkdir -p $$octprefix; \
	  fi; \
	  if [ ! -e $$archprefix ]; then \
	    mkdir -p $$archprefix; \
	  fi; \
	  octave -H -q --no-site-file --eval "warning('off','all');pkg('prefix','$$octprefix','$$archprefix');pkg('global_list',fullfile('$$shareprefix','octave_packages'));pkg('local_list',fullfile('$$shareprefix','octave_packages'));pkg('install','-nodeps','-verbose','$(PACKAGE)');"; \
	  pkgdir=`octave -H -q --no-site-file --eval "warning('off','all');pkg('prefix','$$octprefix','$$archprefix');pkg('global_list',fullfile('$$shareprefix','octave_packages'));pkg('local_list',fullfile('$$shareprefix','octave_packages'));l=pkg('list');disp(l{cellfun(@(x)strcmp(x.name,'$(PKG)'),l)}.dir);$$stripcmd;"`; \
	fi; \
	if [ "X${DISTPKG}X" != "XX" ]; then \
	  if [ -e $$pkgdir/packinfo/on_uninstall.m ]; then \
	    mv $$pkgdir/packinfo/on_uninstall.m \
	       $$pkgdir/packinfo/on_uninstall.m.orig; \
	  fi; \
	  echo "function on_uninstall (desc)" > $$pkgdir/packinfo/on_uninstall.m; \
	  echo "  error ('Can not uninstall %s installed by the $(DISTPKG) package manager', desc.name);" >> $$pkgdir/packinfo/on_uninstall.m; \
	  echo "endfunction" >> $$pkgdir/packinfo/on_uninstall.m; \
	  echo "#! /bin/sh -f" > $$pkgdir/packinfo/dist_admin; \
	  echo "if [ \"\$$1\" == \"install\" ]; then" >> $$pkgdir/packinfo/dist_admin; \
	  echo "  octave -H -q --no-site-file --eval \"pkg('rebuild');\"" >> $$pkgdir/packinfo/dist_admin; \
	  echo "else"  >> $$pkgdir/packinfo/dist_admin; \
	  echo "  pkgdir=\`octave -H -q --no-site-file --eval \"pkg('rebuild');l=pkg('list');disp(l{cellfun(@(x)strcmp(x.name,'$(PKG)'),l)}.dir);\"\`" >> $$pkgdir/packinfo/dist_admin; \
	  echo "  rm \$$pkgdir/packinfo/on_uninstall.m" >> $$pkgdir/packinfo/dist_admin; \
	  echo "  if [ -e \$$pkgdir/packinfo/on_uninstall.m.orig ]; then" >> $$pkgdir/packinfo/dist_admin; \
	  echo "    mv \$$pkgdir/packinfo/on_uninstall.m.orig \$$pkgdir/packinfo/on_uninstall.m" >> $$pkgdir/packinfo/dist_admin; \
	  echo "    cd \$$pkgdir/packinfo" >> $$pkgdir/packinfo/dist_admin; \
	  echo "    octave -q -H --no-site-file --eval \"l=pkg('list');on_uninstall(l{cellfun(@(x)strcmp(x.name,'$(PKG)'),l)});\"" >> $$pkgdir/packinfo/dist_admin; \
	  echo "  fi"  >> $$pkgdir/packinfo/dist_admin; \
	  echo "fi"  >> $$pkgdir/packinfo/dist_admin; \
	  chmod a+x $$pkgdir/packinfo/dist_admin; \
	fi;

clean:
	rm $(PACKAGE)
	$(MAKE) -C src clean
