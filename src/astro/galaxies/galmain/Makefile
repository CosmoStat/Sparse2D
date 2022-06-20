
# @(#)Makefile 1.2 99/02/03 CEA - Jean-Luc Starck

CC = g++-4.0  
# CC = /usr/pubsw/bin/g++
# CC = purify g++

CCC = $(CC) -c
# GCC_FLAG_OPTIM = -/s -O6 -fomit-frame-pointer -funroll-loops
# TIFFDIR=../tiff/sources/libtiff
TIFFDIR=./TIFF

# GCC_FLAG_OPTIM = -g -Wall -DIO_TIFF -I$(TIFFDIR)
GCC_FLAG_OPTIM = -g -Wall
DIR_Cxx_LIB = /network/produits/gnu/lib
# TIFFLIB=$(TIFFDIR)/libtiff.a
DIR_Cxx_LIB = /produits/GNUDev/lib
# LIBCxx      =$(DIR_Cxx_LIB)/libstdc++.a
LIBCxx=-lstdc++ 
LIB_SYS = $(LIBCxx)
# LIB_SYS = $(TIFFLIB) $(LIBCxx) -lm
# LIB_SYS = $(TIFFLIB)  -lm

# GCC_FLAG_OPTIM = -g
CFLAGS=  -DNDEBUG $(GCC_FLAG_OPTIM)
.SUFFIXES: .o .cc .h .o~ .cc~ .h~

OBJ = Genus3D.o Atrou3D.o GetOpt.o IM_IO.o DefMath.o  Memory.o NR_nrutil.o NR_sort.o \
      FFTN_3D.o FFTN.o DefPoint.o CF_Tools.o CUBE_Segment.o CF_LibAna.o


FILES = mr3d_mwt genus coxline simgauss mr3d_genus mr3d_cat2fits cf_ana

current:  mr3d_mwt

all: $(FILES)

.cc.o:
	$(CCC) $(CFLAGS) $(NR_INCLUDE) $(DIR_INCLUDE) $*.cc

.cc~.o:
	$(GET) $*.cc
	$(CCC) $(CFLAGS) $*.cc

.h~.h:
	$(GET) $*.h

genus: genus.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  genus.cc"
	$(CC) $(CFLAGS) $(NR_INCLUDE)  genus.cc -o genus $(OBJ)  $(LIB_SYS)
#	$(CC) $(CFLAGS) $(NR_INCLUDE)  genus.cc -o genus $(OBJ) -ldl -lm $(LIB_SYS)
	@echo "------------------------------"

mr3d_genus: genus.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  genus.cc"
	$(CC) -DWAVELET $(CFLAGS) $(NR_INCLUDE)  genus.cc -o mr3d_genus $(OBJ)  $(LIB_SYS)
#	$(CC) $(CFLAGS) $(NR_INCLUDE)  genus.cc -o genus $(OBJ) -ldl -lm $(LIB_SYS)
	@echo "------------------------------"

mr3d_cat2fits: mr3d_cat2fits.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  mr3d_cat2fits.cc"
	$(CC)  $(CFLAGS) $(NR_INCLUDE)  mr3d_cat2fits.cc -o mr3d_cat2fits $(OBJ)  $(LIB_SYS)
#	$(CC) $(CFLAGS) $(NR_INCLUDE)  genus.cc -o genus $(OBJ) -ldl -lm $(LIB_SYS)
	@echo "------------------------------"

coxline: coxline.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  coxline.cc"
	$(CC) $(CFLAGS) $(NR_INCLUDE)  coxline.cc -o coxline $(OBJ)  $(LIB_SYS) 
	@echo "------------------------------"

simgauss:: simgauss.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  simgauss.cc"
	$(CC) $(CFLAGS) $(NR_INCLUDE)  simgauss.cc -o simgauss $(OBJ)  $(LIB_SYS) 
	@echo "------------------------------"

cf_ana: cf_ana.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  cf_ana.cc"
	$(CC) $(CFLAGS) $(NR_INCLUDE)  cf_ana.cc -o cf_ana $(OBJ)  $(LIB_SYS) 
	@echo "------------------------------"

mr3d_mwt: mr3d_mwt.cc $(OBJ) 
	@echo "------------------------------"
	@echo "Compilation of  mr3d_mwt.cc"
	$(CC) $(CFLAGS) $(NR_INCLUDE)  mr3d_mwt.cc -o mr3d_mwt $(OBJ)  $(LIB_SYS)
	@echo "------------------------------"

clean:
	@echo "Remove: $(FILES) $(OBJ) "
	@rm -f $(FILES) $(OBJ) *% *~ ,*
