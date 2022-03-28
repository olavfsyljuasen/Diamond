include Makefile.local


DIAMONDSOURCES =  dia.C dia.h globalheader.h RunParameter.h


%.run  : %.exec
	./$<

%.exec : %.o 
	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)
	@echo $*.`git describe`$(EXTENSION)
	mv $*.exec $(BINDIR)/$*.`git describe`$(EXTENSION)


dia.o : $(DIAMONDSOURCES)
	$(CCC) $(CCFLAGS) -c -o $@ $<	

dia_shell.o : $(DIAMONDSOURCES) 
	$(CCC) -DBZSHELL $(CCFLAGS) -c -o $@ $<	



spotless:	
	make clean
	rm -f *.ps

clean	:
	rm -f *.dvi
	rm -f *.aux
	rm -f *.log
	rm -f *~
	rm -f core
	rm -f *.o
	rm -f *.exec
	rm -f *.d








