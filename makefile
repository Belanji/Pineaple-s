########### Gnu:
COMPILER = g++
FLAGS=  -Ofast -fipa-pta -fopenmp
#FLAGS=  -O0 -g #-fipa-pta -fopenmp
LIB = -lgsl -lgslcblas -lm -lgomp -pthread
############ Intel:

#COMPILER = mpiicpc
#FLAGS=   -fPIC -wd1572  -O3  -no-prec-div -xHost -simd -fp-model fast=2  
#LIB = -mkl 

############ Intel 2 :

#COMPILER = icpc
#FLAGS= -ipo -O3  -no-prec-div -xAVX -simd -qopenmp -fp-model fast=2 -static
#LIB = -mkl -lgsl 


############ Intel Debug :

#COMPILER = icpc
#FLAGS= -O0 -g -static
#LIB = -mkl -lgsl 


######### For use with petsc:
include ${PETSC_DIR}/lib/petsc/conf/variables

print-%  : ; @echo $* = $($*)


CPPS := $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJS  := $(patsubst src/%.cpp,build/%.o,${CPPS})

Pineaple: ${OBJS}
	@${COMPILER}  ${FLAGS}  ${OBJS} ${LIB} ${PETSC_LIB} -o pineaple
${OBJS}: build/%.o: src/%.cpp | build
	@${COMPILER} -I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include  ${FLAGS} -c $< -o $@
${OBJS}: ${HEADER}
build:
	@mkdir build

clean:
	@rm -f pineaple
	@rm -fr build
	@rm -f *.csv
	@rm -f *.oo
