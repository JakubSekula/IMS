proj = ims2

make:
	@g++ -Wall -std=c++17 $(proj).cpp -o $(proj)

run:
	@make
	@./$(proj) -r 1.0 > sir_out
	@python3 plotit.py

clean:
	rm $(proj)
	rm $(proj)_out