# Makefile for Syndirella

# Configuration
VENV_NAME?=venv
PYTHON=${VENV_NAME}/bin/python

.PHONY: all install clean test

all: install

$(VENV_NAME): requirements.txt
	test -d $(VENV_NAME) || python3 -m venv $(VENV_NAME)
	${PYTHON} -m pip install -U pip
	${PYTHON} -m pip install -r requirements.txt
	touch $(VENV_NAME)

install: $(VENV_NAME)

test: install
	${PYTHON} -m pytest

clean:
	rm -rf $(VENV_NAME)
	find . -type f -name '*.pyc' -delete
	find . -type d -name '__pycache__' -delete
