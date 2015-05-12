FC = gfortran
FFLAGS =  -march=native -O3 -Wall -Wextra -Wtabs -fcheck=all  -mno-avx -ffast-math -fopenmp -I /usr/include

LDFLAGS = -fopenmp
LIBS = -llapack -lblas -g  -lfftw3

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)




OBJS = 
OBJS += modules.o
OBJS += tunneling.o




all: tunneling



tunneling: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<


.PHONY: clean
clean:
	$(RM) tunneling  $(OBJS) *.mod





