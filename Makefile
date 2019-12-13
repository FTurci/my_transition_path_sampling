PROJECT = transition_path_sampling
COMMIT = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
COMMIT_CLEAN = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT_CLEAN) | cut -d " " -f 1)

.PHONY: all dist test version install clean

all: version

version:
	@echo __commit__ = \'$(COMMIT)\' > atooms/${PROJECT}/_commit.py
	@echo __date__ = \'$(DATE)\' >> atooms/${PROJECT}/_commit.py

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

install: version
	python setup.py install

user: version
	python setup.py install --user

develop:
	python setup.py develop --user

clean:
	rm -f atooms/${PROJECT}/*pyc atooms/${PROJECT}/*/*pyc tests/*pyc
