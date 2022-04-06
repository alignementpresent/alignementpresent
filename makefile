CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++11 -fopenmp #-Wno-deprecated-copy #-g -fsanitize=address -fno-omit-frame-pointer

# IFLAGS = -I 
# LFLAGS = -L 

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@ 

main :  main.o PresentData.o PresentData_2R.o PresentData_3R.o PresentData_4R.o permGen.o aux_function.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o main main.o PresentData.o PresentData_2R.o PresentData_3R.o PresentData_4R.o permGen.o aux_function.o $(LFLAGS)

clean :
	rm -rf *.o