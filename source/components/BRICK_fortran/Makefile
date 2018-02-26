# Define "all"
all: obj dais_fastdyn.so simple.so brick_te.so brick_tee.so gsic_magicc.so libbrick.so

# Fortran compiler:
F90 = $(shell which gfortran) # g95

# Directories: objpath, srcpath
OP = ./obj
SP = ./src

# Compiler flags:
Optimization = -O3
SOFlags = -fPIC
#Modpath = -module $(OP)
#Flags = ${SOFlags} ${Optimization} -fno-range-check #$(Modpath)
Flags = ${SOFlags} -fno-range-check #$(Modpath)

# Preprocessor options:
#CPP = cpp
#CPPFLAGS = -DMOC_NOISE

# ==============================================================
## The main programs/subroutines ##
# ==============================================================

## create ./obj folder
obj:
	mkdir -p $@

dais_fastdyn.so: $(OP)/global.o $(OP)/dais_fastdyn.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

simple.so: $(OP)/global.o $(OP)/simple.o $(OP)/run_simple.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

brick_te.so: $(OP)/global.o $(OP)/brick_te.o $(OP)/run_brick_te.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

brick_tee.so: $(OP)/global.o $(OP)/brick_tee.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

gsic_magicc.so: $(OP)/global.o $(OP)/gsic_magicc.o $(OP)/run_gsic_magicc.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

libbrick.so: $(OP)/global.o $(OP)/dais_fastdyn.o $(OP)/brick_te.o $(OP)/brick_tee.o $(OP)/simple.o $(OP)/gsic_magicc.o $(OP)/brick.o $(OP)/run_brick.o
	$(F90) -o $@ $(Flags) -shared $(Incl) $^

# ==============================================================
## Main Object files ##
# ==============================================================

## BRICK ##
$(OP)/run_brick.o: $(OP)/global.o $(OP)/brick.o $(SP)/run_brick_slr.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_brick_slr.f90 -o $@

$(OP)/brick.o: $(OP)/global.o $(OP)/dais_fastdyn.o $(OP)/brick_te.o $(OP)/brick_tee.o $(OP)/simple.o $(OP)/gsic_magicc.o $(SP)/brick_slr.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/brick_slr.f90 -o $@

## DAIS ##
$(OP)/dais_fastdyn.o: $(OP)/global.o $(SP)/dais_fastdyn.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/dais_fastdyn.f90 -o $@

## SIMPLE ##
$(OP)/run_simple.o: $(OP)/global.o $(OP)/simple.o $(SP)/run_simple.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_simple.f90 -o $@

$(OP)/simple.o: $(OP)/global.o $(SP)/simple.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/simple.f90 -o $@

## brick_te ##
$(OP)/run_brick_te.o: $(OP)/global.o $(OP)/brick_te.o $(SP)/run_brick_te.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_brick_te.f90 -o $@

$(OP)/brick_te.o: $(OP)/global.o $(SP)/brick_te.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/brick_te.f90 -o $@

## brick_tee ##
$(OP)/brick_tee.o: $(OP)/global.o $(SP)/brick_tee.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/brick_tee.f90 -o $@

## gsic_magicc ##
$(OP)/run_gsic_magicc.o: $(OP)/global.o $(OP)/gsic_magicc.o $(SP)/run_gsic_magicc.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/run_gsic_magicc.f90 -o $@

$(OP)/gsic_magicc.o: $(OP)/global.o $(SP)/gsic_magicc.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/gsic_magicc.f90 -o $@

## Some global variables used in all/most modules ##
$(OP)/global.o: $(SP)/global.f90
	$(F90) -c $(Flags) $(Incl) $(SP)/global.f90 -o $@

.PHONY: clean

#Cleaning up
clean:
	rm -f ./*.so ./*.mod $(OP)/*.o $(OP)/*.mod $(OP)/*.F90
