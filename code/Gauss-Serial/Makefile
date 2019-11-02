
FC = gfortran
FFLAGS = -fbounds-check -fopenmp
LFLAGS = $(FFLAGS)

.PHONY: plots, clean, clobber

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

GS_Serial2d.exe:  GS_Serial2d_main_omp.o GS_Serial2d_sub_omp.o
	$(FC) $(LFLAGS) GS_Serial2d_main_omp.o GS_Serial2d_sub_omp.o -o $@

solution.txt: GS_Serial2d.exe
	@echo 
	@echo Running code...
	./GS_Serial2d.exe

plots: solution.txt
	@echo 
	@echo Plotting results...
	python plot_solution.py

clean:
	rm -f *.o *.exe

clobber: clean
	rm -f solution.txt *.png

