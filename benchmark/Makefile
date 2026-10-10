SHELL := /bin/bash
.ONESHELL:

.PHONY: all extract dbs search m8extract sensitivity querylib querylib-extract querylib-search querylib-analyze querylib-plot analyze clean

all: extract dbs search sensitivity m8extract analyze

extract:
	bash benchmark.sh extract

dbs:
	bash benchmark.sh dbs

search:
	bash benchmark.sh search

m8extract:
	bash benchmark.sh m8extract

sensitivity:
	bash benchmark.sh sensitivity

querylib:
	bash benchmark.sh querylib

querylib-extract:
	bash benchmark.sh querylib-extract

querylib-search:
	bash benchmark.sh querylib-search

querylib-analyze:
	bash benchmark.sh querylib-analyze

querylib-plot:
	bash benchmark.sh querylib-plot

analyze:
	bash benchmark.sh analyze

clean:
	bash benchmark.sh clean
