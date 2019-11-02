
FC = gfortran
FFLAGS = -fbounds-check -fopenmp
LFLAGS = $(FFLAGS)

.PHONY: plots, clean, clobber

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

GS_RB2d.exe:  GS_RB2d_main_omp.o GS_RB2d_sub_omp.o
	$(FC) $(LFLAGS) GS_RB2d_main_omp.o GS_RB2d_sub_omp.o -o $@

solution.txt: GS_RB2d.exe
	@echo 
	@echo Running code...
	./GS_RB2d.exe

plots: solution.txt
	@echo 
	@echo Plotting results...
	python plot_solution.py

clean:
	rm -f *.o *.exe

clobber: clean
	rm -f solution.txt *.png

