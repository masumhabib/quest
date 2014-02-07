#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=mpicc
CCC=mpic++
CXX=mpic++
FC=mpif90
AS=as

# Macros
CND_PLATFORM=MPI-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/atoms/Atoms.o \
	${OBJECTDIR}/src/atoms/Lattice.o \
	${OBJECTDIR}/src/negf/NEGF.o \
	${OBJECTDIR}/src/negf/NegfEloop.o \
	${OBJECTDIR}/src/negf/computegs.o \
	${OBJECTDIR}/src/parallel/parloop.o \
	${OBJECTDIR}/src/potential/linearPot.o \
	${OBJECTDIR}/src/potential/potential.o \
	${OBJECTDIR}/src/potential/terminal.o \
	${OBJECTDIR}/src/python/pyqmicad.o \
	${OBJECTDIR}/src/qm/kp/genKpAtoms.o \
	${OBJECTDIR}/src/qm/kp/graphenekp.o \
	${OBJECTDIR}/src/qm/kp/tikp.o \
	${OBJECTDIR}/src/simulations/Device.o \
	${OBJECTDIR}/src/utils/ConsoleProgressBar.o \
	${OBJECTDIR}/src/utils/stringutils.o


# C Compiler Flags
CFLAGS=-m64

# CC Compiler Flags
CCFLAGS=-m64
CXXFLAGS=-m64

# Fortran Compiler Flags
FFLAGS=-m64

# Assembler Flags
ASFLAGS=--64

# Link Libraries and Options
LDLIBSOPTIONS=-Wl,--start-group /home/kh8hk/usr/lib/libboost_mpi_python.a /home/kh8hk/usr/lib/libboost_mpi.a -Wl,--end-group -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -Wl,-Bdynamic -L/home/kh8hk/usr/lib -lboost_serialization -L/home/kh8hk/usr/lib -lboost_python -L/opt/intel/lib/intel64 -liomp5 -lpthread -lm -ldl -Wl,--no-as-needed

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk lib/qmicad.${CND_DLIB_EXT}

lib/qmicad.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p lib
	${LINK.cc} -o lib/qmicad.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/atoms/Atoms.o: src/atoms/Atoms.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/atoms
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atoms/Atoms.o src/atoms/Atoms.cpp

${OBJECTDIR}/src/atoms/Lattice.o: src/atoms/Lattice.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/atoms
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atoms/Lattice.o src/atoms/Lattice.cpp

${OBJECTDIR}/src/negf/NEGF.o: src/negf/NEGF.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/negf
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/negf/NEGF.o src/negf/NEGF.cpp

${OBJECTDIR}/src/negf/NegfEloop.o: src/negf/NegfEloop.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/negf
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/negf/NegfEloop.o src/negf/NegfEloop.cpp

${OBJECTDIR}/src/negf/computegs.o: src/negf/computegs.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/negf
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/negf/computegs.o src/negf/computegs.cpp

${OBJECTDIR}/src/parallel/parloop.o: src/parallel/parloop.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/parallel
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/parallel/parloop.o src/parallel/parloop.cpp

${OBJECTDIR}/src/potential/linearPot.o: src/potential/linearPot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/potential
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/potential/linearPot.o src/potential/linearPot.cpp

${OBJECTDIR}/src/potential/potential.o: src/potential/potential.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/potential
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/potential/potential.o src/potential/potential.cpp

${OBJECTDIR}/src/potential/terminal.o: src/potential/terminal.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/potential
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/potential/terminal.o src/potential/terminal.cpp

${OBJECTDIR}/src/python/pyqmicad.o: src/python/pyqmicad.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/python
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/python/pyqmicad.o src/python/pyqmicad.cpp

${OBJECTDIR}/src/qm/kp/genKpAtoms.o: src/qm/kp/genKpAtoms.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/qm/kp
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qm/kp/genKpAtoms.o src/qm/kp/genKpAtoms.cpp

${OBJECTDIR}/src/qm/kp/graphenekp.o: src/qm/kp/graphenekp.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/qm/kp
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qm/kp/graphenekp.o src/qm/kp/graphenekp.cpp

${OBJECTDIR}/src/qm/kp/tikp.o: src/qm/kp/tikp.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/qm/kp
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qm/kp/tikp.o src/qm/kp/tikp.cpp

${OBJECTDIR}/src/simulations/Device.o: src/simulations/Device.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/simulations
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/simulations/Device.o src/simulations/Device.cpp

${OBJECTDIR}/src/utils/ConsoleProgressBar.o: src/utils/ConsoleProgressBar.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/utils
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/utils/ConsoleProgressBar.o src/utils/ConsoleProgressBar.cpp

${OBJECTDIR}/src/utils/stringutils.o: src/utils/stringutils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/utils
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/usr/include/python2.6 -I/home/kh8hk/usr/include -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/utils/stringutils.o src/utils/stringutils.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} lib/qmicad.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
