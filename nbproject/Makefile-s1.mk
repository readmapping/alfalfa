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
CC=gcc.exe
CCC=g++.exe
CXX=g++.exe
FC=gfortran
AS=as.exe

# Macros
CND_PLATFORM=Cygwin-Windows
CND_CONF=s1
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/_ext/1942341467/sparseSA.o \
	${OBJECTDIR}/_ext/1942341467/utils.o \
	${OBJECTDIR}/_ext/1942341467/fasta.o \
	${OBJECTDIR}/_ext/1942341467/qsufsort.o \
	${OBJECTDIR}/_ext/1942341467/mapper.o \
	${OBJECTDIR}/_ext/1942341467/dp.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-O3 -DNDEBUG
CXXFLAGS=-O3 -DNDEBUG

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/alfalfa.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/alfalfa.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/alfalfa ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/_ext/1942341467/sparseSA.o: /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/sparseSA.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1942341467
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1942341467/sparseSA.o /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/sparseSA.cpp

${OBJECTDIR}/_ext/1942341467/utils.o: /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/utils.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1942341467
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1942341467/utils.o /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/utils.cpp

${OBJECTDIR}/_ext/1942341467/fasta.o: /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/fasta.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1942341467
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1942341467/fasta.o /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/fasta.cpp

${OBJECTDIR}/_ext/1942341467/qsufsort.o: /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/qsufsort.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/1942341467
	${RM} $@.d
	$(COMPILE.c) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1942341467/qsufsort.o /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/qsufsort.c

${OBJECTDIR}/_ext/1942341467/mapper.o: /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/mapper.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1942341467
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1942341467/mapper.o /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/mapper.cpp

${OBJECTDIR}/_ext/1942341467/dp.o: /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/dp.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1942341467
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1942341467/dp.o /cygdrive/C/Users/mvyvermn/Documents/NetBeansProjects/alfalfa/dp.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/alfalfa.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
