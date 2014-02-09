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
CND_CONF=Release
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
	${OBJECTDIR}/src/python/PyNegfEloop.o \
	${OBJECTDIR}/src/python/PyNegfParams.o \
	${OBJECTDIR}/src/python/pyqmicad.o \
	${OBJECTDIR}/src/qm/kp/graphenekp.o \
	${OBJECTDIR}/src/qm/kp/tikp.o \
	${OBJECTDIR}/src/simulations/Device.o \
	${OBJECTDIR}/src/string/stringutils.o \
	${OBJECTDIR}/src/utils/ConsoleProgressBar.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libqmicad.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libqmicad.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libqmicad.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/atoms/Atoms.o: src/atoms/Atoms.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/atoms
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atoms/Atoms.o src/atoms/Atoms.cpp

${OBJECTDIR}/src/atoms/Lattice.o: src/atoms/Lattice.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/atoms
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/atoms/Lattice.o src/atoms/Lattice.cpp

${OBJECTDIR}/src/negf/NEGF.o: src/negf/NEGF.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/negf
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/negf/NEGF.o src/negf/NEGF.cpp

${OBJECTDIR}/src/negf/NegfEloop.o: src/negf/NegfEloop.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/negf
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/negf/NegfEloop.o src/negf/NegfEloop.cpp

${OBJECTDIR}/src/negf/computegs.o: src/negf/computegs.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/negf
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/negf/computegs.o src/negf/computegs.cpp

${OBJECTDIR}/src/parallel/parloop.o: src/parallel/parloop.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/parallel
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/parallel/parloop.o src/parallel/parloop.cpp

${OBJECTDIR}/src/potential/linearPot.o: src/potential/linearPot.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/potential
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/potential/linearPot.o src/potential/linearPot.cpp

${OBJECTDIR}/src/potential/potential.o: src/potential/potential.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/potential
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/potential/potential.o src/potential/potential.cpp

${OBJECTDIR}/src/potential/terminal.o: src/potential/terminal.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/potential
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/potential/terminal.o src/potential/terminal.cpp

${OBJECTDIR}/src/python/PyNegfEloop.o: src/python/PyNegfEloop.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/python
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/python/PyNegfEloop.o src/python/PyNegfEloop.cpp

${OBJECTDIR}/src/python/PyNegfParams.o: src/python/PyNegfParams.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/python
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/python/PyNegfParams.o src/python/PyNegfParams.cpp

${OBJECTDIR}/src/python/pyqmicad.o: src/python/pyqmicad.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/python
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/python/pyqmicad.o src/python/pyqmicad.cpp

${OBJECTDIR}/src/qm/kp/graphenekp.o: src/qm/kp/graphenekp.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/qm/kp
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qm/kp/graphenekp.o src/qm/kp/graphenekp.cpp

${OBJECTDIR}/src/qm/kp/tikp.o: src/qm/kp/tikp.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/qm/kp
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qm/kp/tikp.o src/qm/kp/tikp.cpp

${OBJECTDIR}/src/simulations/Device.o: src/simulations/Device.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/simulations
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/simulations/Device.o src/simulations/Device.cpp

${OBJECTDIR}/src/string/stringutils.o: src/string/stringutils.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/string
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/string/stringutils.o src/string/stringutils.cpp

${OBJECTDIR}/src/utils/ConsoleProgressBar.o: src/utils/ConsoleProgressBar.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/utils
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/utils/ConsoleProgressBar.o src/utils/ConsoleProgressBar.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libqmicad.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
