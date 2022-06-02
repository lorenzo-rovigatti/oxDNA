.SUFFIXES: .cu

# change these two lines to match your CUDA installation
SOURCEDIR = ../../
# 64 or 32 bit?
CPU_ARCH = 64
SOURCES = $(SOURCEDIR)Interaction.cpp  $(SOURCEDIR)Particle.cpp  $(SOURCEDIR)SimBackend.cpp  $(SOURCEDIR)Utils.cpp  $(SOURCEDIR)IOManager.cpp PrintDataManager.cpp  $(SOURCEDIR)MDBackend.cpp  $(SOURCEDIR)MD_CPUBackend.cpp  $(SOURCEDIR)SimManager.cpp ProcessData_Backend.cpp  $(SOURCEDIR)MCBackend.cpp  $(SOURCEDIR)ExternalForce.cpp  $(SOURCEDIR)MC_CPUBackend.cpp $(SOURCEDIR)CMMC_CPUBackend.cpp $(SOURCEDIR)OrderParameters.cpp $(SOURCEDIR)Histogram.cpp $(SOURCEDIR)Weights.cpp $(SOURCEDIR)FFS_MD_CPUBackend.cpp $(SOURCEDIR)FL_CMMC_CPUBackend.cpp $(SOURCEDIR)VMMC_CPUBackend.cpp 
OBJS =  $(SOURCEDIR)/parse_input/getline.o $(SOURCEDIR)/parse_input/parse_input.o $(SOURCEDIR)/timing/timing.o $(SOURCEDIR)/time_scales/time_scales.o
EXEC = output_bonds
CC = g++  
DEFINES = -DNOCUDA 
INCLUDES = -I$(SOURCEDIR)
FLAGS = -O3 -ffast-math
localSOURCES = PrintDataManager.cpp ProcessData_Backend.cpp
$(EXEC): $(SOURCES) 
	 $(CC) -o $(EXEC) $(INCLUDES) $(OBJS) $(SOURCES) $(DEFINES) $(FLAGS)
clean: 
	rm $(EXEC)
