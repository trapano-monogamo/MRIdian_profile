activation:
	Scripts\activate

run: activation
	python src\main.py

.PHONY: clear

clear:
	rmdir /s ./out/*
	rmdir /s ./src/__pycache__
