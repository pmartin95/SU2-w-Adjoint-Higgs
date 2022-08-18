CC=g++ -std=c++17
SRCDIR=src
OBJDIR=obj
BINDIR=bin
DEPDIR=includes
SRCS=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
DEP=$(wildcard $(DEPDIR)/*.hpp)

INC=-I$(DEPDIR)/eigen -I$(DEPDIR)
FLAGS=-lm -fopenmp -O4# -ggdb3 -pg

all: bin/main.out

bin/main: $(OBJS) $(DEP)
	@echo "Linking objects into main.exe..."
	@$(CC) $(INC) $(FLAGS) -o $(BINDIR)/main $(OBJS)

obj/%.o : src/%.cpp $(DEP)
	@echo "Creating $(<:src/%.cpp=%) object file..."
	@$(CC) $(INC) $(FLAGS) -c -o $@ $<

.PHONY:clean 

clean:
	@rm -rf ./obj/*.o ./bin/*.exe .bin/*.out

.PHONY:run 

run : bin/main.out
	@echo "Executing binary..."
	@./bin/main.out	

.PHONY:timerun

timerun : ./bin/main.out	
	@echo "Executing a timed binary..."
	@time ./bin/main.out	

.PHONY:profile

profile : ./bin/main.out gmon.out
	gprof ./main.exe gmon.out > ./dat/profdata.txt 

graph : 
	@echo "Graphing results..."
	@gnuplot dat/graph_plaq.p

.PHONY:cleanconf
cleanconf :
	@rm -rf ./configurations/*.bin 
	
