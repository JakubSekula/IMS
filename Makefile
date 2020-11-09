proj = sir

make:
	@g++ -Wall -Wextra -std=c++17 $(proj).cpp -o $(proj)
	pip3 install pandas
	pip3 install matplotlib
run:
	@./$(proj) -r 1.0 > sir.out
	@python3 plotit.py
clean:
	rm sir
	rm sir.out