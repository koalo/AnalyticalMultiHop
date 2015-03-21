all: model preprocessor

CFLAGS	         = -O3 -Isrc/common/ -MMD
CPPFLAGS         = $(CFLAGS) -std=c++11 -Wfatal-errors -g
LIBS		 = -lboost_program_options -lboost_system -lboost_filesystem -lboost_graph

SRC_common = src/common/Experiment.cpp src/common/Route.cpp  src/common/Relations.cpp src/common/Topology.cpp
SRC_model = $(SRC_common) src/model/model.cpp src/model/Calculation.cpp src/model/RelationsGenerator.cpp src/model/ResultWriter.cpp
SRC_preprocessor = $(SRC_common) src/preprocessor/preprocessor.cpp src/preprocessor/TopologyGenerator.cpp src/preprocessor/ConnectionsGenerator.cpp src/preprocessor/RouteGenerator.cpp src/common/Connections.cpp

OBJECTS_model = $(SRC_model:.cpp=.o)
OBJECTS_preprocessor = $(addprefix obj/preprocessor/,$(SRC_preprocessor:.cpp=.o))

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test

model: $(OBJECTS_model) chkopts
	-${CLINKER} -o model $(OBJECTS_model) ${PETSC_DM_LIB} $(LIBS)

preprocessor: $(OBJECTS_preprocessor) chkopts
	-${CLINKER} -o preprocessor $(OBJECTS_preprocessor) ${PETSC_DM_LIB} $(LIBS)

obj/preprocessor/%.o: %.cpp
	mkdir -p $(dir $@)
	g++ $(CPPFLAGS) -c -o $@ $< $(LIBS)

cl:
	rm -f $(OBJECTS_preprocessor) $(OBJECTS_model)

-include $(SRC_model:.cpp=.d)
-include $(addprefix obj/preprocessor/,$(SRC_preprocessor:.cpp=.d))

