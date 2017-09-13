all: clean compile meshmonk

# Object files which need to be linked together
TARGETS = build/meshmonk.o \
build/BaseCorrespondenceFilter.o \
build/CorrespondenceFilter.o \
build/Downsampler.o \
build/helper_functions.o \
build/InlierDetector.o \
build/NeighbourFinder.o \
build/NonrigidRegistration.o \
build/PyramidNonrigidRegistration.o \
build/RigidRegistration.o \
build/RigidTransformer.o \
build/ScaleShifter.o \
build/SymmetricCorrespondenceFilter.o \
build/ViscoElasticTransformer.o

# Compile flags
## Flags to compile
M_FLAGS = --verbose -Wall -fexceptions -O2 -Wall -std=c++14 -g -Wl,-V -fPIC -I /usr/local/include/ -I vendor -c

## Flags to build the example
M_FLAGS2 = --verbose -Wall -fexceptions -O2 -Wall -std=c++14 -g -fPIC -I /usr/local/include/ -I vendor

# Build the meshmonk library for OSX. The output will be a .dynlib file
# Copy/paste the lib: cp libmeshmonk.dylib /usr/local/libe
meshmonk: $(TARGETS)
	g++ -shared $(TARGETS) -dynamiclib -o libmeshmonk.dylib -lOpenMeshCore -lOpenMeshTools -L/usr/local/lib

# Compile all the files explicitly
compile:
	mkdir -p build
	g++ $(M_FLAGS) meshmonk.cpp -o build/meshmonk.o
	g++ $(M_FLAGS) src/BaseCorrespondenceFilter.cpp -o build/BaseCorrespondenceFilter.o
	g++ $(M_FLAGS) src/CorrespondenceFilter.cpp -o build/CorrespondenceFilter.o
	g++ $(M_FLAGS) src/Downsampler.cpp -o build/Downsampler.o
	g++ $(M_FLAGS) src/helper_functions.cpp -o build/helper_functions.o
	g++ $(M_FLAGS) src/InlierDetector.cpp -o build/InlierDetector.o
	g++ $(M_FLAGS) src/NeighbourFinder.cpp -o build/NeighbourFinder.o
	g++ $(M_FLAGS) src/NonrigidRegistration.cpp -o build/NonrigidRegistration.o
	g++ $(M_FLAGS) src/PyramidNonrigidRegistration.cpp -o build/PyramidNonrigidRegistration.o
	g++ $(M_FLAGS) src/RigidRegistration.cpp -o build/RigidRegistration.o
	g++ $(M_FLAGS) src/RigidTransformer.cpp -o build/RigidTransformer.o
	g++ $(M_FLAGS) src/ScaleShifter.cpp -o build/ScaleShifter.o
	g++ $(M_FLAGS) src/SymmetricCorrespondenceFilter.cpp -o build/SymmetricCorrespondenceFilter.o
	g++ $(M_FLAGS) src/ViscoElasticTransformer.cpp -o build/ViscoElasticTransformer.o

# Build the example
# Run: ./example
example:
	g++ $(M_FLAGS2) -lOpenMeshCore -lOpenMeshTools -lmeshmonk example.cpp -o example

# Clean all the .o and .dynlib files
clean:
	rm -f build/*.o libmeshmonk.dylib
