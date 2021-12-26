#####################################################
#   CompHorner Make File   							#
#####################################################
include make.inc

install_horner_test:
	@$(CCC) $(CCOPT) -o horner_test C/horner_test.c $(CCLNFLAGS)
	
install_illcond_test:
	@$(CCC) $(CCOPT) -o illcond_test C/illcond_test.c $(CCLNFLAGS)
	
install_specpoly_test:
	@$(CCC) $(CCOPT) -o specpoly_test C/specpoly_test.c $(CCLNFLAGS)
	
install_wellcond_test:
	@$(CCC) $(CCOPT) -o wellcond_test C/wellcond_test.c $(CCLNFLAGS)
	
pdf_horner_test:
	@pdflatex --interaction=batchmode TeX/horner_errbound.tex
	@mv horner_errbound.pdf figures/horner_errbound.pdf
	@rm horner_errbound.log
	@rm horner_errbound.aux
	@open figures/horner_errbound.pdf
	
	@pdflatex TeX/horner_time.tex
	@mv horner_time.pdf figures/horner_time.pdf
	@rm horner_time.log
	@rm horner_time.aux
	@open figures/horner_time.pdf
	
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
	
pdf_specpoly_test:
	@pdflatex --interaction=batchmode TeX/specpoly_test.tex
	@mv specpoly_test.pdf figures/specpoly_test.pdf
	@rm specpoly_test.log
	@rm specpoly_test.aux
	@open figures/specpoly_test.pdf

	@pdflatex --interaction=batchmode TeX/limit_acc.tex
	@mv limit_acc.pdf figures/limit_acc.pdf
	@rm limit_acc.log
	@rm limit_acc.aux
	@open figures/limit_acc.pdf
	
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
	
run_all: install_horner_test install_illcond_test install_specpoly_test install_wellcond_test
	./horner_test 5
	./illcond_test
	./specpoly_test
	./wellcond_test
	
tex_all: pdf_horner_test pdf_illcond_test pdf_specpoly_test pdf_wellcond_test
	
uninstall:
	@rm -f horner_test
	@rm -f illcond_test
	@rm -f specpoly_test
	@rm -f wellcond_test