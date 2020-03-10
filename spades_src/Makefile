  
#This is for building spades (Used for my MULKSG project)

all: targz spades 

targz:
	./make-targz.sh
	tar -xzvf SPAdes-3.13.0.tar.gz 

spades:
	cd SPAdes-3.13.0/ ; \
	./spades_compile.sh ; \
	export PATH=$(PATH):$(PWD)/bin ;\
	cd .. 

clean:
	rm -rf SPAdes-3.13.0 SPAdes-3.13.0.tar.gz 
