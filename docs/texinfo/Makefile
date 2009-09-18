all: html pdf

html : lucee-ref.texi array-classes.texi implementation.texi
	makeinfo --html lucee-ref.texi -o html

pdf : lucee-ref.texi array-classes.texi implementation.texi
	pdftex lucee-ref.texi
	texindex lucee-ref.??
	pdftex lucee-ref.texi
	texindex lucee-ref.??
	pdftex lucee-ref.texi

clean:
	rm -rf lucee-ref.pdf html *.dvi *.tp *.ky *.toc *.vt *.aux *.cp *.fn *.pg *.vr *.log *~