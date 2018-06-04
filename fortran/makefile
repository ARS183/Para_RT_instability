f77=mpif90
#f77 = ifort
opt= -O3   
srs= openNS2d.f90  openNS2d_mpi.f90  VanLeer.f90  \
	openNS2d_init_RT.f90  openNS2d_write_data.f90 \
	openNS3d_DF_AECDS.f90 openNS3d_DF_BOUND.f90 \
	 openNS2d_solver.f90 Du_UCC.f90

OBJS=$(srs:.f90=.o)

%.o:%.f90
	$(f77) $(opt) -c $<

default: $(OBJS)
	$(f77) -O3  -o openNS2d.out $(OBJS)
	rm -f *.o 
clean:
	rm -f *.out *.o 
