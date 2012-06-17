# Name of the file(s) we're compiling

TARGET = rossler

# flags
# for 64 bit
CPPFLAGS = -fPIC
OPT = -O3
#DEBUG = -g


# PYTHON = /usr/include/python2.6
# EXTERN = /home/jberwald/Dropbox/Projects/diffusion

# INSTALL_LIB = /usr/local/lib/python2.6/dist-packages/jb/diffusion 

# location of CAPD libs and include directories
CAPD = /Users/jberwald/src/capd/capdDynSys/include

# $(OBJECT).so: $(TARGET2).o
#         g++ -shared -Wall,--export-dynamic \
#         $(TARGET2).o -L$(BOOST_LIB) -L/usr/lib/python2.6/config \
#         -lpython2.6 -lboost_python -llapack -lann \
#         -o $(OBJECT).so

$(TARGET).o: $(TARGET2).cpp $(TARGET1).o
	g++ $(OPT) $(CPPFLAGS) $(DEBUG) -I$(CAPD) -I$(EXTERN) -I$(BOOST_INC) -c $(TARGET).cpp

clean:
        rm -f *.o *.dll *.a *~ core
