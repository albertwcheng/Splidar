################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Commonf.cpp \
../ConfigFile.cpp \
../ConservedBlock.cpp \
../Debugger.cpp \
../GffEntry.cpp \
../RandomAccessFile.cpp \
../SpliceMMGraph.cpp \
../SplidarGraph.cpp \
../SuperConfigFile.cpp \
../snip.cpp \
../snipMain.cpp 

OBJS += \
./Commonf.o \
./ConfigFile.o \
./ConservedBlock.o \
./Debugger.o \
./GffEntry.o \
./RandomAccessFile.o \
./SpliceMMGraph.o \
./SplidarGraph.o \
./SuperConfigFile.o \
./snip.o \
./snipMain.o 

CPP_DEPS += \
./Commonf.d \
./ConfigFile.d \
./ConservedBlock.d \
./Debugger.d \
./GffEntry.d \
./RandomAccessFile.d \
./SpliceMMGraph.d \
./SplidarGraph.d \
./SuperConfigFile.d \
./snip.d \
./snipMain.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include -I/usr/include/python2.5/ -I../samtools-0.1.8/ -O3 -Wall -c -fmessage-length=0 -Wno-reorder -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<" #-Wno-builtin-macro-redefined 
	@echo 'Finished building: $<'
	@echo ' '


