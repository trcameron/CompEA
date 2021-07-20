#####################################################
#   CompHorner Make File   							#
#####################################################
include make.inc

install_horner_test:
	@$(CCC) $(CCOPT) $(CCLNFLAGS) -o horner_test C/horner_test.c
	
install_illcond_test:
	@$(CCC) $(CCOPT) $(CCLNFLAGS) -o illcond_test C/illcond_test.c
	
install_specpoly_test:
	@$(CCC) $(CCOPT) $(CCLNFLAGS) -o specpoly_test C/specpoly_test.c
	
install_wellcond_test:
	@$(CCC) $(CCOPT) $(CCLNFLAGS) -o wellcond_test C/wellcond_test.c
	
pdf_horner_test:
	@pdflatex --interaction=batchmode TeX/horner_test.tex
	@mv horner_test.pdf figures/horner_test.pdf
	@rm horner_test.log
	@rm horner_test.aux
	@open figures/horner_test.pdf
	
pdf_illcond_test:
	@pdflatex --interaction=batchmode TeX/illcond_err.tex
	@mv illcond_err.pdf figures/illcond_err.pdf
	@rm illcond_err.log
	@rm illcond_err.aux
	@open figures/illcond_err.pdf
	
	@pdflatex --interaction=batchmode TeX/illcond_time.tex
	@mv illcond_time.pdf figures/illcond_time.pdf
	@rm illcond_time.log
	@rm illcond_time.aux
	@open figures/illcond_time.pdf
	
pdf_wellcond_test:
	@pdflatex --interaction=batchmode TeX/wellcond_err.tex
	@mv wellcond_err.pdf figures/wellcond_err.pdf
	@rm wellcond_err.log
	@rm wellcond_err.aux
	@open figures/wellcond_err.pdf
	
	@pdflatex --interaction=batchmode TeX/wellcond_time.tex
	@mv wellcond_time.pdf figures/wellcond_time.pdf
	@rm wellcond_time.log
	@rm wellcond_time.aux
	@open figures/wellcond_time.pdf
	
pdf_specpoly_test:
	@pdflatex TeX/specpoly_test.tex
	@mv specpoly_test.pdf figures/specpoly_test.pdf
	@rm specpoly_test.log
	@rm specpoly_test.aux
	@open figures/specpoly_test.pdf
	
uninstall:
	@rm -f horner_test
	@rm -f illcond_test
	@rm -f specpoly_test
	@rm -f wellcond_test