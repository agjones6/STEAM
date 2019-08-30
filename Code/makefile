# To check dependencies, run [gfortran -M -cpp *.f95] untill it runs clean
# Note that each command line must be tabbed
# '-c' -- do not discard comments



scode: main.o globals.o state.o files.o controllers.o conditional.o
	gfortran -o scode -llapack -lblas main.o globals.o state.o heat.o \
	files.o controllers.o conditional.o safety.o steamdome.o \
	pressuredrop.o matrices.o startup.o; rm -f *.o *.mod

state.o state.mod: state.f95
	gfortran -c state.f95

files.o files.mod: files.f95
	gfortran -c files.f95

steamdome.o steamdome.mod: steamdome.f95 state.o
	gfortran -c steamdome.f95

globals.o globals.mod: globals.f95 state.o files.o
	gfortran -c globals.f95

conditional.o conditional.mod: conditional.f95 globals.o state.o steamdome.o files.o
	gfortran -c conditional.f95

controllers.o controllers.mod: controllers.f95 globals.o state.o files.o
	gfortran -c controllers.f95

pressuredrop.o pressuredrop.mod: pressuredrop.f95 state.o globals.o \
    conditional.o
	gfortran -c pressuredrop.f95

matrices.o matrices.mod: matrices.f95 state.o globals.o conditional.o pressuredrop.o
	gfortran -c matrices.f95

safety.o safety.mod: safety.f95 state.o globals.o pressuredrop.o
	gfortran -c safety.f95

startup.o startup.mod: startup.f95
	gfortran -c startup.f95

heat.o heat.mod: heat.f95 globals.o state.o
	gfortran -c heat.f95

main.o: main.f95 globals.o state.o files.o controllers.o conditional.o \
	 	steamdome.o pressuredrop.o matrices.o safety.o startup.o heat.o
	gfortran -c main.f95

.PHONY: clean

clean:
	rm -f *.o *.mod
