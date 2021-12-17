# path_planning
This project implements an ant colony optimization algorithm for finding the shortest path in a grid with surmountable (but highly weighted) barriers. The serial implementation is ant.cpp, and the parallelized version is ant_par.cpp. The latter uses pthreads but is otherwise unchanged from the former.

to run:
make
./ant
./ant_par

This should run on *nix systems and probably on windows.