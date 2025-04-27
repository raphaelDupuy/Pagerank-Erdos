run: Main
	./Main

Main: main.c
	gcc main.c -Wall -o Main

clean:
	rm -f Main
	ls -l
