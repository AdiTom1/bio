CC = g++

CCFLAGS = -Wall -O2
LDFLAGS = -lstdc++ -lm

# User defined classes and modules. (no file suffixes)

CLASSES = Vector3 Atom Rotation3 Matrix3 RigidTrans3 Triangle PDB Match numerics alignRand

# Prepare object and source file list using pattern substitution func.
ALL  = $(CLASSES)
OBJS = $(patsubst %, %.o,  $(ALL))
SRCS = $(patsubst %, %.cc, $(ALL))

TARGET = alignRand.Linux

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)  -o $(TARGET)

%.o: %.cc
	$(CC) $(CCFLAGS) -c $*.cc

clean:
	/bin/rm -f *.o *~ \#* core

depend:
	makedepend -- $(CCFLAGS) -- $(SRCS)
# DO NOT DELETE THIS LINE -- make depend depends on it.
