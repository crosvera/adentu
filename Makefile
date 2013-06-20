# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# These flags will override any settings
ifeq ($(i386),1)
    OS_SIZE = 32
    OS_ARCH = i686
endif

ifeq ($(x86_64),1)
    OS_SIZE = 64
    OS_ARCH = x86_64
endif

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

# Location of the CUDA Toolkit binaries and libraries
CUDA_PATH       ?= /Developer/NVIDIA/CUDA-5.0
CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
ifneq ($(DARWIN),)
    CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib
else
    ifeq ($(OS_SIZE),32)
        CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib
    else
        CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib64
    endif
endif

# Common binaries
NVCC            ?= $(CUDA_BIN_PATH)/nvcc
GCC             ?= gcc

# Extra user flags
EXTRA_NVCCFLAGS ?=
EXTRA_LDFLAGS   ?=
EXTRA_CCFLAGS   ?= -std=c99 


# CUDA code generation flags
GENCODE_SM10    := -gencode arch=compute_10,code=sm_10
GENCODE_SM13    := -gencode arch=compute_13,code=sm_13
GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35
GENCODE_FLAGS   := $(GENCODE_SM20) $(GENCODE_SM30)


# OS-specific build flags
ifneq ($(DARWIN),) 
      LDFLAGS   := -Xlinker -rpath $(CUDA_LIB_PATH) -L$(CUDA_LIB_PATH) -lcudart -lcufft -lcurand -framework OpenGL -framework GLUT
      CCFLAGS   := -arch $(OS_ARCH) 
else
  ifeq ($(OS_SIZE),32)
      LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart -lcufft -lcurand -lglut -lGLU
      CCFLAGS   := -m32
  else
      LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart -lcufft -lcurand -lglut -lGLU
      CCFLAGS   := -m64
  endif
endif


# OS-architecture specific flags
ifeq ($(OS_SIZE),32)
      NVCCFLAGS := -m32
else
      NVCCFLAGS := -m64
endif

dbg = 1

# Debug build flags
ifeq ($(dbg),1)
      CCFLAGS   += -g
      NVCCFLAGS += -g -G
      TARGET    := debug
else
      TARGET    := release
endif


#GLib
CCFLAGS += `pkg-config --cflags glib-2.0`
NVCCFLAGS += `pkg-config --cflags glib-2.0`
LDFLAGS += `pkg-config --libs glib-2.0` -lm


# Common includes and paths for CUDA and Adentu libs.
INCLUDES      := -I$(CUDA_INC_PATH) -I./src -I./src/event  



# Target rules

all: build

build: adentu


adentu-utils: vec3.h adentu-utils.h

vec3.o: src/vec3.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

vec3-cuda.o: src/vec3-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-atom.o: src/adentu-atom.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-atom-cuda.o: src/adentu-atom-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-grid.o: src/adentu-grid.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-grid-cuda.o: src/adentu-grid-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-model.o: src/adentu-model.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-neighbourhood.o: src/adentu-neighbourhood.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-neighbourhood-cuda.o: src/adentu-neighbourhood-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

# ADENTU EVENTS
adentu-event.o: src/adentu-event.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-bc.o: src/event/adentu-event-bc.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-bc-cuda.o: src/event/adentu-event-bc-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-mpc.o: src/event/adentu-event-mpc.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-mpc-cuda.o: src/event/adentu-event-mpc-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-gfc.o: src/event/adentu-event-gfc.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-gfc-cuda.o: src/event/adentu-event-gfc-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-ggc.o: src/event/adentu-event-ggc.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-ggc-cuda.o: src/event/adentu-event-ggc-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-event-usr.o: src/event/adentu-event-usr.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<


adentu-runnable.o: src/adentu-runnable.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<

# ADENTU GRAPHICS
adentu-graphic.o: src/adentu-graphic.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<



adentu-cuda-utils.o: src/adentu-cuda-utils.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu.o: src/adentu.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<



# USR MODULES
adentu-usr-atoms-pos-cuda.o: src/usr/atoms-pos-cuda.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ -c $<

adentu-usr-print-event-info.o: src/usr/print-event-info.c
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $(INCLUDES) -o $@ -c $<





adentu: vec3.o vec3-cuda.o adentu.o adentu-grid.o adentu-atom.o  adentu-model.o adentu-atom-cuda.o adentu-grid-cuda.o adentu-neighbourhood.o adentu-neighbourhood-cuda.o adentu-event.o adentu-event-bc.o adentu-event-bc-cuda.o adentu-event-mpc.o adentu-event-mpc-cuda.o adentu-runnable.o adentu-graphic.o adentu-cuda-utils.o adentu-event-gfc-cuda.o adentu-event-gfc.o adentu-event-ggc-cuda.o adentu-event-ggc.o adentu-event-usr.o    adentu-usr-atoms-pos-cuda.o adentu-usr-print-event-info.o 
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) -o $@ $+ $(LDFLAGS)  $(EXTRA_LDFLAGS)


run: build
	./adentu

clean:
	rm -f adentu *.o
