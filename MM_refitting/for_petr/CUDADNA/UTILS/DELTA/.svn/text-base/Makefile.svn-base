CC = cc
LINKER = g++
CFLAGS = -Wall -Wshadow
LFLAGS = 
LIBS = -lm
DEFINES = -DBUILD_TIME="\"`date +'%D %R'`\"" -DSVN_VERSION="\"`svnversion`\"" -DNOCUDA
INCLUDES = "-I../../common"

ifeq ($(dbg), 1)
	CFLAGS += -O0 -g3 -pg
else
	CFLAGS += -O3 -ffast-math

	ifeq ($(g), 1)
		CFLAGS += -g3 -pg
		LFLAGS += -g3 -pg
	else
		DEFINES += -DNDEBUG
	endif
endif

ifeq ($(arch), 20)
	CUFLAGS += -arch=sm_20
else
	CUFLAGS += -arch=sm_13
	DEFINES += -DOLD_ARCH
endif

C_OBJS = IOManager.o MD_CPUBackend.o \
SimBackend.o SimManager.o Utils.o Particle.o \
MDBackend.o MCBackend.o MC_CPUBackend.o \
Interaction.o ExternalForce.o Strand.o

C_OBJS_DELTA = delta.o DeltaManager.o DeltaBackend.o $(C_OBJS)
C_OBJS_CSD = csd.o CSDManager.o CSDBackend.o $(C_OBJS)

EXT_OBJS = timing.o parse_input.o time_scales.o

#all: delta

delta: externals delta.o $(C_OBJS_DELTA)
	$(LINKER) -o delta $(C_OBJS_DELTA) $(EXT_OBJS) $(LFLAGS) $(LIBS) $(INCLUDES) $(DEFINES)
	
csd: externals csd.o $(C_OBJS_CSD)
	$(LINKER) -o csd $(C_OBJS_CSD) $(EXT_OBJS) $(LFLAGS) $(LIBS) $(INCLUDES) $(DEFINES)

%.o: ../../%.cpp ../../%.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
csd.o: csd.cpp
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
CSDManager.o: CSDManager.cpp CSDManager.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
CSDBackend.o: CSDBackend.cpp CSDBackend.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
delta.o: delta.cpp
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
DeltaManager.o: DeltaManager.cpp DeltaManager.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
DeltaBackend.o: DeltaBackend.cpp DeltaBackend.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)
	
Strand.o: Strand.cpp Strand.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)

Interaction.o: ../../Interaction.cpp ../../Interaction.h ../../model.h
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDES) $(DEFINES)

externals:
	$(CC) -o timing.o -c ../../timing/timing.c $(CFLAGS)
	$(CC) -o parse_input.o -c ../../parse_input/parse_input.c $(CFLAGS)
	$(CC) -o time_scales.o -c ../../time_scales/time_scales.c $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf delta csd
