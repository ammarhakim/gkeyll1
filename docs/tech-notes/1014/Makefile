NAME=1014-multifluid-equil
$(NAME).pdf: $(NAME).tex
	pdflatex $(NAME).tex
	bibtex $(NAME)
	pdflatex $(NAME).tex
	pdflatex $(NAME).tex

clean:
	rm -rf $(NAME).pdf *~ $(NAME).aux $(NAME).bbl $(NAME).blg $(NAME).log