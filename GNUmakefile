DATADIR= /scratch/sylju
#HOMEDIR= /mn/felt/u1/sylju
HOMEDIR= $(HOME)
ABELDIR= $(HOME)
SAGADIR= $(HOME)
OSXDIR = $(HOME)

%.x64 : force
	cp Makefile.x64 Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

%.osx : force
	cp Makefile.osx Makefile.local
	@echo `git describe`
	make -f makefile $*.exec


%.abel : force
	cp Makefile.abel Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

%.saga : force
	cp Makefile.saga Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

# The envirmoment variable MYMAKEFILE contains the name of the appropriate Makefile.local
%    :  force
	cp $(MYMAKEFILE) Makefile.local 
	make -f makefile $*




force : ;

triangular_x64: triangular.x64 triangular_corr.x64 triangular_obc_corr.x64 triangular_fixed_corr.x64 triangular_corr_noadia.x64  triangular_obc_corr_noadia.x64 triangular_corr_reduced.x64

square_x64: square.x64 square_corr.x64

triangular_saga: triangular.saga triangular_corr.saga triangular_obc_corr.saga triangular_fixed_corr.saga triangular_fixed.saga triangular_corr_noadia.saga  triangular_obc_corr_noadia.saga triangular_corr_reduced.saga 

square_saga: square.saga square_corr.saga






