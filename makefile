default: all done
parallel: par done

# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -Wall -O3
CC_PAR = g++-12
CFLAGS_PAR = -std=c++17 -Wall -O3 -fopenmp
LIBFLAG = -larmadillo

# Directories
SRC_DIR = src
OBJ_DIR = objectFiles
OBJ_DIR_PAR = objectFilesPar

# Files
MAIN_FILE = $(SRC_DIR)/main.cpp
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
OBJ_FILES_PAR = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR_PAR)/%.o,$(SRC_FILES))

# Target default (serial)
EXECUTABLE = $(SRC_DIR)/main
# Target parallel
EXECUTABLE_PAR = $(SRC_DIR)/main_par

# Rules for default (serial)
$(EXECUTABLE): $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBFLAG)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

# Rules for parallel
$(EXECUTABLE_PAR): $(OBJ_FILES_PAR)
	$(CC_PAR) $(CFLAGS_PAR) -o $@ $^ -I/usr/local/include -L/usr/local/lib -larmadillo
#$(LIBFLAG)

$(OBJ_DIR_PAR)/%.o: $(SRC_DIR)/%.cpp
	$(CC_PAR) $(CFLAGS_PAR) -c -o $@ $< -I/usr/local/include -L/usr/local/lib -larmadillo

# Phony targets
.PHONY: clean clean_par

# Clean up serial (default)
clean:
	rm -f $(OBJ_DIR)/*.o $(EXECUTABLE)

# Clean up parallel
clean_par:
	rm -f $(OBJ_DIR_PAR)/*.o $(EXECUTABLE_PAR)

# Run the program in serial (default)
all: $(EXECUTABLE)
	./$(EXECUTABLE)

# Run the program in parallel
par: $(EXECUTABLE_PAR)
	./$(EXECUTABLE_PAR)

# Default target for serial
build: $(EXECUTABLE)

# Default target for parallel
build_par: $(EXECUTABLE_PAR)


###############OLD MAKEFILE################


# .DEFAULT_GOAL := default
# CC = g++
# CFLAG = --std=c++17 -Wall -O3 -larmadillo
# LIBFLAG = -larmadillo

# default: compile move link run done

	
# all: compile move link run clean done

# build_main: compile_main move_main link run

# compile_main: 
# 	g++ $(CFLAG) -c src/main.cpp

# compile:
# 	g++ $(CFLAG) -c $(wildcard src/*.cpp) $(LIBFLAG)

# link:
# 	g++ $(CFLAG) -o src/main $(wildcard objectFiles/*.o) $(LIBFLAG)

# run:
# 	./src/main

# clean:
# 	rm -f *.o 

# move_main:
# 	mv main.o objectFiles

# move:
# 	mv *.o objectFiles

done:
	@echo "Done! \(@^0^@)/"


# Create a new project with all the necessary folders and files
create_project: make_dir make_files

make_dir:
	mkdir -p data include objectFiles objectFilesPar scripts src tex

make_files:
	touch README.md src/main.cpp data/.keep include/.keep objectFiles/.keep objectFilesPar/.keep scripts/.keep tex/.keep