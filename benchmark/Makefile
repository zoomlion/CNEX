SHELL := /bin/bash
.ONESHELL:

.PHONY: all extract dbs search m8extract analyze clean

all: extract dbs search m8extract analyze

extract:
	bash benchmark.sh extract

dbs:
	bash benchmark.sh dbs

search:
	bash benchmark.sh search

m8extract:
	bash benchmark.sh m8extract

analyze:
	bash benchmark.sh analyze

clean:
	bash benchmark.sh clean
