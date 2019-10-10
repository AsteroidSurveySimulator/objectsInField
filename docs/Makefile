filename=manual
pdfreader=open

pdf: ps
	ps2pdf ${filename}.ps

pdf-print: ps
	ps2pdf -dColorConversionStrategy=/LeaveColorUnchanged -dPDFSETTINGS=/printer ${filename}.ps

ps:	dvi
	dvips -t letter ${filename}.dvi

dvi:
	latex ${filename}
	bibtex ${filename}
	latex ${filename}
	latex ${filename}

read:
	${pdfreader} ${filename}.pdf &

clean:
	rm -f ${filename}.{ps,pdf,log,aux,out,dvi,bbl,blg}
