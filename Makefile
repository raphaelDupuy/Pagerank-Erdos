run: Td3
	./Td3

Td3: td3.c
	gcc td3.c -Wall -o Td3

clean:
	rm -f Td3
	ls -l
