all: html pdf

html : lucee-ref.texi
	makeinfo --html lucee-ref.texi -o html

pdf : lucee-ref.texi
	pdftex lucee-ref.texi

clean:
	rm -rf lucee-ref.pdf html *.dvi *.tp *.ky *.toc *.vt *.aux *.cp *.fn *.pg *.vr