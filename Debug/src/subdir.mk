################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/vt.cpp

C_SRCS += \
../src/functions.c \
../src/finitevol.c \
../src/closure.c \
../src/bc.c \
../src/init.c \
../src/interface.c \
../src/stream.c \
../src/collide.c \
../src/macro.c \
../src/spongezones.c \
../src/main.c 

OBJS += \
./src/vtk.o \
./src/functions.o \
./src/finitevol.o \
./src/closure.o \
./src/bc.o \
./src/init.o \
./src/interface.o \
./src/stream.o \
./src/collide.o \
./src/macro.o \
./src/spongezones.o \
./src/main.o \

C_DEPS += \
./src/functions.d \
./src/finitevol.d \
./src/closure.d \
./src/bc.d \
./src/init.d \
./src/interface.d \
./src/stream.d \
./src/collide.d \
./src/macro.d \
./src/spongezones.d \
./src/main.d 

CPP_DEPS += \
./src/vtk.d


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DD3Q19 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


