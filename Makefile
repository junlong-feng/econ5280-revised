.PHONY: pdf live clean

pdf:
	mkdir -p output/pdf
	cd latex && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	cp latex/main.pdf output/pdf/ECON5280_Revised_Lecture_Notes.pdf

live:
	cd live && quarto render

clean:
	cd latex && latexmk -C
	rm -rf live/_site
