# use variables in submakes

OPTIMIZATION_ON ?= 1
ASAN ?= 0
DEPRECATION_ON ?= 1
DEBUG ?= 0
COPYCHECK_ARGS ?=
LLD ?= 0
WERROR ?= 0

CXX ?= cxx


CXXFLAGS := -std=c++17 -Wall -Wextra -O2
OPTFLAGS :=

SABER_SRC_DIRS := $(shell find src -type d)
SRC_DIRS = $(SABER_SRC_DIRS)

SABER_CPP_FILES := $(foreach dir,$(SABER_SRC_DIRS),$(wildcard $(dir)/*.cpp))

CPP_FILES += $(SABER_CPP_FILES) 
O_FILES   := $(foreach f,$(CPP_FILES:.cpp=.o),build/$f)

# create build directories
$(shell mkdir -p $(foreach dir,$(SRC_DIRS),build/$(dir)))


# Main targets
all: Saber.elf 
clean:
	rm -r build\

distclean: clean libclean

rebuild: clean all

format:
	clang-format-11 -i $(SABER_CPP_FILES)

.PHONY: all  clean rebuild format

build/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $(OUTPUT_OPTION) $<

# Linking
Saber.elf: $(O_FILES)
	$(CXX) $(CXXFLAGS) $(O_FILES) $(OUTPUT_OPTION) -L/usr/local/lib -lcryptopp