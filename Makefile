PROJECT = sample

.PHONY: all dist test version install clean

all: version

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > transition_path_sampling/_commit.py
	@echo __date__ = \'$(DATE)\' >> transition_path_sampling/_commit.py

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

install: version
	python setup.py install --user

develop:
	python setup.py develop --user

clean:
	rm -f ${PROJECT}/*pyc ${PROJECT}/*/*pyc tests/*pyc
