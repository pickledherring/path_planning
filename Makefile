todo: ant ant_par
ant: ant.cpp
	g++ -o ant ant.cpp -std=c++2a
ant_par: ant_par.cpp
	g++ -pthread -o ant_par ant_par.cpp -std=c++2a
clean:
	rm ant ant_par