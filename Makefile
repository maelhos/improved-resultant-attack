
.PHONY: default
default: lib

# library files used for all attacks
cpp_lib_source_files := $(shell find lib/src -name *.cpp)
cpp_lib_object_files := $(patsubst lib/src/%.cpp, lib/build/%.o, $(cpp_lib_source_files))
LIB_INCLUDE_PATH := -I lib/include

# rescue prime attack
cpp_rescue_source_files := $(shell find rescue/src -name *.cpp)
cpp_rescue_object_files := $(patsubst rescue/src/%.cpp, rescue/build/%.o, $(cpp_rescue_source_files))
RESCUE_INCLUDE_PATH := -I rescue/include

# anemoi attack
cpp_anemoi_source_files := $(shell find anemoi/src -name *.cpp)
cpp_anemoi_object_files := $(patsubst anemoi/src/%.cpp, anemoi/build/%.o, $(cpp_anemoi_source_files))
ANEMOI_INCLUDE_PATH := -I anemoi/include

# griffin attack
cpp_griffin_source_files := $(shell find griffin/src -name *.cpp)
cpp_griffin_object_files := $(patsubst griffin/src/%.cpp, griffin/build/%.o, $(cpp_griffin_source_files))
GRIFFIN_INCLUDE_PATH := -I griffin/include

CPP_INCLUDE_PATHS := -I ntl/include

OUT := attack

LIBS := -lgmp -lm -lssl -lcrypto
OPTS := -march=native -O3 -ffast-math -std=c++20 -L ntl/src/ntl.a
CC := g++

###########################################

# lib files
$(cpp_lib_object_files): lib/build/%.o : lib/src/%.cpp
	@mkdir -p $(dir $@) && \
	$(CC) $(OPTS) -c -Wall $(LIB_INCLUDE_PATH) $(CPP_INCLUDE_PATHS) $(patsubst lib/build/%.o, lib/src/%.cpp, $@) -o $@ $(LIBS)

# anemoi files
$(cpp_anemoi_object_files): anemoi/build/%.o : anemoi/src/%.cpp
	@mkdir -p $(dir $@) && \
	$(CC) $(OPTS) -c -Wall $(LIB_INCLUDE_PATH) $(ANEMOI_INCLUDE_PATH) $(CPP_INCLUDE_PATHS) $(patsubst anemoi/build/%.o, anemoi/src/%.cpp, $@) -o $@ $(LIBS)

# rescue files
$(cpp_rescue_object_files): rescue/build/%.o : rescue/src/%.cpp
	@mkdir -p $(dir $@) && \
	$(CC) $(OPTS) -c -Wall $(LIB_INCLUDE_PATH) $(RESCUE_INCLUDE_PATH) $(CPP_INCLUDE_PATHS) $(patsubst rescue/build/%.o, rescue/src/%.cpp, $@) -o $@ $(LIBS)

# griffin files
$(cpp_griffin_object_files): griffin/build/%.o : griffin/src/%.cpp
	@mkdir -p $(dir $@) && \
	$(CC) $(OPTS) -c -Wall $(LIB_INCLUDE_PATH) $(GRIFFIN_INCLUDE_PATH) $(CPP_INCLUDE_PATHS) $(patsubst griffin/build/%.o, griffin/src/%.cpp, $@) -o $@ $(LIBS)

###########################################

.PHONY: anemoi
anemoi: $(cpp_lib_object_files) $(cpp_anemoi_object_files)
	@$(CC) $(OPTS) -Wall $(CPP_INCLUDE_PATHS) $(LIB_INCLUDE_PATH) -o $(OUT) $(cpp_lib_object_files) $(cpp_anemoi_object_files) ntl/src/ntl.a $(LIBS)

.PHONY: rescue
rescue: $(cpp_lib_object_files) $(cpp_rescue_object_files)
	@$(CC) $(OPTS) -Wall $(CPP_INCLUDE_PATHS) $(LIB_INCLUDE_PATH) -o $(OUT) $(cpp_lib_object_files) $(cpp_rescue_object_files) ntl/src/ntl.a $(LIBS)

.PHONY: griffin
griffin: $(cpp_lib_object_files) $(cpp_griffin_object_files)
	@$(CC) $(OPTS) -Wall $(CPP_INCLUDE_PATHS) $(LIB_INCLUDE_PATH) -o $(OUT) $(cpp_lib_object_files) $(cpp_griffin_object_files) ntl/src/ntl.a $(LIBS)

clean:
	@rm -rfv lib/build
	@mkdir lib/build

	@rm -rfv anemoi/build
	@mkdir anemoi/build

	@rm -rfv rescue/build
	@mkdir rescue/build

	@rm -rfv griffin/build
	@mkdir griffin/build

	@echo "CLEAN"