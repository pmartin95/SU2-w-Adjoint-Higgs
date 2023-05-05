CXX=g++ -std=c++17
SRCDIR=src
OBJDIR=obj
BINDIR=bin
DEPDIR=includes
SRCS=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
DEP=$(wildcard $(DEPDIR)/*.hpp)
BOOST=-lboost_system -lboost_filesystem
INC=-I$(DEPDIR)/eigen -I$(DEPDIR)
CXXFLAGS=-lm -fopenmp -O3 # -ggdb3 -pg
TARGET=$(BINDIR)/main

all: $(TARGET)

$(TARGET): $(OBJS) $(DEP)
	@echo "Linking objects into main.exe..."
	@$(CXX) $(INC) $(CXXFLAGS) -o $@ $(OBJS) $(BOOST)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(DEP)
	@echo "Creating $(<:$(SRCDIR)/%.cpp=%) object file..."
	@$(CXX) $(INC) $(CXXFLAGS) -c -o $@ $< $(BOOST)

.PHONY: clean
clean:
	@rm -rf ./obj/*.o ./bin/*.exe .bin/*.out $(TARGET)

.PHONY: run
run: $(TARGET).out
	@echo "Executing binary..."
	@./$(TARGET).out

.PHONY: timerun
timerun: $(TARGET).out
	@echo "Executing a timed binary..."
	@time ./$(TARGET).out

.PHONY: profile
profile: $(TARGET).out gmon.out
	gprof $(TARGET) gmon.out > ./dat/profdata.txt

.PHONY: graph
graph:
	@echo "Graphing results..."
	@gnuplot dat/graph_plaq.p

.PHONY: cleanconf
cleanconf:
	@rm -rf ./configurations/*.bin
