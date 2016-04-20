################################################################################
#
# Makefile grand central
#
################################################################################

TOPDIR   = ../..
BUILDDIR = .

VPATH = $(BUILDDIR)

include $(TOPDIR)/sysdep/Makedef.$(NICAM_SYS)
include $(TOPDIR)/Mkinclude

all:
	$(MAKE) makedir
	$(MAKE) makebin

allclean: clean cleandir

PRJ = dckernel_divdamp.exe

PRG = main

MODS = \
	mod_precision.o      \
	mod_adm.o            \
	mod_debug.o          \
	mod_grd.o            \
	mod_gmtr.o           \
	mod_vmtr.o           \
	mod_oprt.o           \
	mod_oprt3d.o

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(NICAM_BINDIR)

makebin: makedir $(PRJ)



mod_precision.o      : mod_precision.f90
mod_adm.o            : mod_adm.f90
mod_debug.o          : mod_debug.f90
mod_grd.o            : mod_grd.f90
mod_gmtr.o           : mod_gmtr.f90
mod_vmtr.o           : mod_vmtr.f90
mod_oprt.o           : mod_oprt.f90
mod_oprt3d.o         : mod_oprt3d.f90

main.o : main.f90 $(MODS)

$(BUILDDIR)/$(PRJ) : $(PRG).o $(BUILDDIR)/$(MODS)
	$(LD) $(LFLAGS) -o $@ $^
	$(INSTALL) $@ $(NICAM_BINDIR)
	@echo; echo "##### PROGRAM $(notdir $@) is installed Successfully #####"; echo



cleandir:
	-rm $(NICAM_BINDIR)/*

clean:
	rm -f $(PRJ)
	rm -f *~ *.mod *.o *.lst *.L
	cd $(BUILDDIR); rm -f $(PRJ); rm -f *~ *.mod *.o *.lst *.L

.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) $(MODDIROPT) -o $(OBJDIR)/$@ -c $<

%.o: %.mod

$(BUILDDIR)/%.mod: %.f90
	$(make) $(patsubst %.f90,$(BUILDDIR)/%.o,$<)

.PHONY : clean cleandir allclean
