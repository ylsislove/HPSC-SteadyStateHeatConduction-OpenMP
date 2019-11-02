
FC = gfortran
FFLAGS = -fbounds-check -fopenmp
LFLAGS = $(FFLAGS)

.PHONY: plots, clean, clobber

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

GS_WF_SOR2d.exe:  GS_WF_SOR2d_main_omp.o GS_WF_SOR2d_sub_omp.o
	$(FC) $(LFLAGS) GS_WF_SOR2d_main_omp.o GS_WF_SOR2d_sub_omp.o -o $@

solution.txt: GS_WF_SOR2d.exe
	@echo 
	@echo Running code...
	./GS_WF_SOR2d.exe

plots: solution.txt
	@echo 
	@echo Plotting results...
	python plot_solution.py

clean:
	rm -f *.o *.exe

clobber: clean
	rm -f solution.txt *.png

