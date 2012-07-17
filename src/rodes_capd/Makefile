# Last modified: Fri Jul 14 2000
#
# A simple Makefile for the programs.

all :
	@-echo ""; exit 0
	@echo "************* RODES - Rigorous ODE Solver *************"
	@echo "";
	@echo "You must specify 'make target', where the 'target' can be:"
	@echo ""
	@echo "    rodes        (the main program in my paper)"
	@echo "";
	@echo "    expansion    (used for estimating the accumulated expansion)"
	@echo "";
	@echo "    smalldiv     (estimates the modulus of the small divisors)"
	@echo "";
	@echo "    coeff        (estimates the size of the coefficients)"
	@echo "";

# -----------------------------------------------------------------------

HERE    = ./
R_EFILE = $(HERE)/rodes
E_EFILE = $(HERE)/expansion
S_EFILE = $(HERE)/smalldiv
C_EFILE = $(HERE)/coeff
LOC     = /usr
INCLS   = -I$(LOC)/include
CXX	= g++

# Use this line for Linux.
#CFLAGS  = -Wall -O2 #-D__I386__

# Use this line for Solaris. 
#CFLAGS  = -w -O2   

# For CAPD
CAPDBINDIR =/Users/jberwald/src/capd/bin/

# setting compiler and linker flags
# CAPDFLAGS = `${CAPDBINDIR}capd-config --cflags`
CAPDFLAGS = -D__USE_FILIB__ -frounding-math -ffloat-store -I/Users/jberwald/src/capd/capdExt/./filibsrc -I/Users/jberwald/src/capd/capdExt/filibsrc -I/Users/jberwald/src/capd/capdExt/./boost/include -I/Users/jberwald/src/capd/capdDynSys/include -I/Users/jberwald/src/capd/capdAux/include -I/Users/jberwald/src/capd/capdAlg/include -I/Users/jberwald/src/capd/capdExt/include -I/Users/jberwald/src/capd/capdRedHom/include

#CAPDLIBS = `${CAPDBINDIR}capd-config --libs`
CAPDLIBS = /Users/jberwald/src/capd/capdExt/boost/src/libcapd_boost_program_options.a /Users/jberwald/src/capd/capdExt/boost/src/libcapd_boost_serialization.a /Users/jberwald/src/capd/capdExt/boost/src/libcapd_boost_filesystem.a /Users/jberwald/src/capd/capdExt/boost/src/libcapd_boost_system.a /Users/jberwald/src/capd/capdDynSys/src/capd/dynsysfacade/libcapddynsysfacade.a /Users/jberwald/src/capd/capdDynSys/src/capd/covrel/libcapdcovrel.a /Users/jberwald/src/capd/capdDynSys/src/capd/diffIncl/libcapddiffIncl.a /Users/jberwald/src/capd/capdDynSys/src/capd/poincare/libcapdpoincare.a /Users/jberwald/src/capd/capdDynSys/src/capd/diffAlgebra/libcapddiffAlgebra.a /Users/jberwald/src/capd/capdDynSys/src/capd/dynsys/libcapddynsys.a /Users/jberwald/src/capd/capdDynSys/src/capd/map/libcapdmap.a /Users/jberwald/src/capd/capdDynSys/src/capd/dynset/libcapddynset.a /Users/jberwald/src/capd/capdDynSys/src/capd/geomset/libcapdgeomset.a /Users/jberwald/src/capd/capdAux/src/capd/auxil/libcapdauxil.a /Users/jberwald/src/capd/capdAlg/src/capd/basicalg/libcapdbasicalg.a /Users/jberwald/src/capd/capdAlg/src/capd/algfacade/libcapdalgfacade.a /Users/jberwald/src/capd/capdAlg/src/capd/matrixAlgorithms/libcapdmatrixAlgorithms.a /Users/jberwald/src/capd/capdAlg/src/capd/rounding/libcapdrounding.a /Users/jberwald/src/capd/capdAlg/src/capd/vectalg/libcapdvectalg.a /Users/jberwald/src/capd/capdAlg/src/capd/intervals/libcapdintervals.a /Users/jberwald/src/capd/capdExt/src/capd/alglib/libcapdalglib.a /Users/jberwald/src/capd/capdExt/src/capd/bzip2/libcapdbzip2.a /Users/jberwald/src/capd/capdExt/src/capd/chom/libcapdchom.a /Users/jberwald/src/capd/capdExt/src/capd/homology/libcapdhomology.a /Users/jberwald/src/capd/capdRedHom/src/capd/bitSet/libcapdbitSet.a /Users/jberwald/src/capd/capdRedHom/src/capd/cubSet/libcapdcubSet.a /Users/jberwald/src/capd/capdRedHom/src/capd/homAlgebra/libcapdhomAlgebra.a /Users/jberwald/src/capd/capdRedHom/src/capd/redAlg/libcapdredAlg.a /Users/jberwald/src/capd/capdRedHom/src/capd/repSet/libcapdrepSet.a -L/Users/jberwald/src/capd/capdExt/filibsrc/libprim/.libs -lprim 

CXXFLAGS = $(CAPDFLAGS) $(INCLS) -I./include -O2 -Wall



# Enter the path for the PROFIL/BIAS directory here.
#PROFIL = /home/warwick/RODES/IA/MY_PROFIL

#P_INCL = -I$(PROFIL)/incl -IBIAS
#P_LIBS = -L$(PROFIL)/lib -lProfil -LBIAS -lBias

# -----------------------------------------------------------------------

R_OBJS   = classes.o  fixed_point.o vector_field.o low_functions.o \
	   flow_functions.o return_map.o convert.o request.o rodes.o

# -----------------------------------------------------------------------

E_OBJS   = classes.o request.o expansion.o 

# -----------------------------------------------------------------------

S_MAIN   = smalldiv.cc

# -----------------------------------------------------------------------

C_MAIN   = coeff.cc

# -----------------------------------------------------------------------

clean:
	rm -rf *~ *.o coeff expansion rodes smalldiv

# -----------------------------------------------------------------------

rodes: $(R_OBJS)
	@echo "Linking to CAPD..."
	@echo "  compiler args: "
	$(CXX) $(CXXFLAGS) -o $(R_EFILE) $(R_OBJS) $(CAPDLIBS)
	@echo "                    ... Pronto!"

# -----------------------------------------------------------------------

expansion: $(E_OBJS)
	@echo "Linking to CAPD..."
	@$(CXX) $(CXXFLAGS) -o $(E_EFILE) $(R_OBJS) $(CAPDLIBS)	
	@echo "                    ... Pronto!"

# -----------------------------------------------------------------------

smalldiv: $(S_MAIN)
	@echo "Linking to CAPD..."
	@$(CXX) $(CXXFLAGS) -o $(S_EFILE) $(R_OBJS) $(CAPDLIBS)
	@echo "                    ... Pronto!"

# -----------------------------------------------------------------------

coeff: $(C_MAIN)
	@echo "Linking to CAPD..."
	@$(CXX) $(CXXFLAGS) -o $(C_EFILE) $(R_OBJS) $(CAPDLIBS)
	@echo "                    ... Pronto!"

# -----------------------------------------------------------------------

classes.o: classes.cc  classes.h 
	@echo "Updating 'classes.o'"
	$(CXX) $(CAPDFLAGS)-MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 	

request.o: request.cc  request.h
	@echo "Updating 'request.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

vector_field.o: vector_field.cc  vector_field.h \
	        classes.cc  classes.h   \
	        error_handler.h 
	@echo "Updating 'vector_field.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

fixed_point.o: fixed_point.cc  fixed_point.h \
	       classes.cc  classes.h  \
	       list.h
	@echo "Updating 'fixed_point.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

low_functions.o: low_functions.cc low_functions.h \
	         classes.cc  classes.h   \
	         error_handler.h \
	         vector_field.cc  vector_field.h
	@echo "Updating 'low_functions.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

flow_functions.o: flow_functions.cc flow_functions.h \
	          list.h  error_handler.h \
	          classes.cc  classes.h   \
	          low_functions.cc low_functions.h \
	          vector_field.cc  vector_field.h
	@echo "Updating 'flow_functions.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

return_map.o: return_map.cc  return_map.h \
	      classes.cc  classes.h \
	       list.h  error_handler.h \
	      flow_functions.cc flow_functions.h \
	      fixed_point.cc  fixed_point.h \
	      low_functions.cc low_functions.h
	@echo "Updating 'return_map.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

convert.o: convert.cc convert.h classes.cc  classes.h \
	   list.h 2d_classes.h \
	   flow_functions.cc flow_functions.h 
	@echo "Updating 'convert.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 


rodes.o:  rodes.cc \
	 2d_classes.h list.h \
	 error_handler.h \
	 classes.cc  classes.h \
	 convert.cc  convert.h \
	 request.cc  request.h \
	 return_map.cc  return_map.h \
	 fixed_point.cc  fixed_point.h
	@echo "Updating 'rodes.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

expansion.o: expansion.cc \
	     classes.cc  classes.h \
	     2d_classes.h \
	     request.cc  request.h  list.h 
	@echo "Updating 'expansion.o'"
	$(CXX) $(CAPDFLAGS) -MT $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

smalldiv.o: smalldiv.cc 
	@echo "Updating 'smalldiv.o'"
	$(CXX) $(CAPDFLAGS) $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

coeff.o: coeff.cc 
	@echo "Updating 'coeff.o'"
	$(CXX) $(CAPDFLAGS) $@ -MD -MP -MF ${@:%=%.d} -c $< -o $@ 

# -----------------------------------------------------------------------

