all: tdma_model tdma_inverse tdma_single csma_model preprocessor

CFLAGS	         = -O3 -Isrc/common/ -MMD
CPPFLAGS         = $(CFLAGS) -std=c++11 -Wfatal-errors -Wno-literal-suffix -g0
LIBS		 = -lboost_program_options -lboost_system -lboost_filesystem -lboost_graph

SRC_common = src/common/Route.cpp src/common/Experiment.cpp src/common/Relations.cpp src/common/Topology.cpp src/common/TDMASchedule.cpp
SRC_csma_model = $(SRC_common) src/csma_model/csma_model.cpp src/csma_model/Calculation.cpp src/csma_model/RelationsGenerator.cpp src/csma_model/ResultWriter.cpp
SRC_tdma_model = $(SRC_common) src/tdma_model/tdma_model.cpp src/tdma_model/Calculation.cpp src/tdma_model/ResultWriter.cpp src/tdma_model/Queue.cpp
SRC_tdma_inverse = $(SRC_common) src/tdma_model/tdma_inverse.cpp src/tdma_model/Queue.cpp
SRC_tdma_single = $(SRC_common) src/tdma_model/tdma_single.cpp src/tdma_model/Queue.cpp
SRC_preprocessor = $(SRC_common) src/preprocessor/preprocessor.cpp src/preprocessor/TopologyGenerator.cpp src/preprocessor/ConnectionsGenerator.cpp src/preprocessor/RouteGenerator.cpp src/preprocessor/TDMAGenerator.cpp src/common/Connections.cpp

OBJECTS_csma_model = $(SRC_csma_model:.cpp=.o)
OBJECTS_tdma_model = $(SRC_tdma_model:.cpp=.o)
OBJECTS_tdma_inverse = $(SRC_tdma_inverse:.cpp=.o)
OBJECTS_tdma_single = $(SRC_tdma_single:.cpp=.o)
OBJECTS_preprocessor = $(addprefix obj/preprocessor/,$(SRC_preprocessor:.cpp=.o))

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

csma_model: $(OBJECTS_csma_model) | chkopts
	-${CLINKER} -o csma_model $(OBJECTS_csma_model) ${PETSC_DM_LIB} $(LIBS)

tdma_model: $(OBJECTS_tdma_model) | chkopts
	-${CLINKER} -o tdma_model $(OBJECTS_tdma_model) ${PETSC_DM_LIB} $(LIBS)

tdma_inverse: $(OBJECTS_tdma_inverse) | chkopts
	-${CLINKER} -o tdma_inverse $(OBJECTS_tdma_inverse) ${PETSC_DM_LIB} $(LIBS)

tdma_single: $(OBJECTS_tdma_single) | chkopts
	-${CLINKER} -o tdma_single $(OBJECTS_tdma_single) ${PETSC_DM_LIB} $(LIBS)

preprocessor: $(OBJECTS_preprocessor) | chkopts
	-${CLINKER} -o preprocessor $(OBJECTS_preprocessor) ${PETSC_DM_LIB} $(LIBS)

obj/preprocessor/%.o: %.cpp
	mkdir -p $(dir $@)
	g++ $(CPPFLAGS) -c -o $@ $< $(LIBS)

cl:
	rm -f rm -f tdma_single tdma_inverse tdma_model preprocessor csma_model $(OBJECTS_preprocessor) $(OBJECTS_csma_model) $(OBJECTS_tdma_model) $(OBJECTS_tdma_inverse) $(OBJECTS_tdma_single) $(SRC_csma_model:.cpp=.d) $(SRC_tdma_inverse:.cpp=.d) $(SRC_tdma_single:.cpp=.d) $(SRC_tdma_model:.cpp=.d) $(SRC_preprocessor:.cpp=.d)

-include $(SRC_csma_model:.cpp=.d)
-include $(SRC_tdma_model:.cpp=.d)
-include $(SRC_tdma_single:.cpp=.d)
-include $(SRC_tdma_inverse:.cpp=.d)
-include $(addprefix obj/preprocessor/,$(SRC_preprocessor:.cpp=.d))

