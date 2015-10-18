CFLAGS= -O2 -arch=sm_13
MPI_LIBS= -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl  -lnsl -lutil -lm -ldl
MPI_INCLUDES= -I/usr/include/openmpi/
all: fourier_transform_vec

fourier_transform_vec : fourier_transform_vec.o host_fourier_transform.o device_fourier_transform.o #device_fourier_transform.cubin
	nvcc ${CFLAGS} -o fourier_transform_vec fourier_transform_vec.o  device_fourier_transform.o host_fourier_transform.o ${MPI_LIBS} 

fourier_transform_vec.o : fourier_transform_vec.cu 
	nvcc ${CFLAGS} -c -o fourier_transform_vec.o fourier_transform_vec.cu ${MPI_INCLUDES}

host_fourier_transform.o: host_fourier_transform.cu host_fourier_transform.h
	nvcc ${CFLAGS} -c -o host_fourier_transform.o host_fourier_transform.cu ${MPI_INCLUDES}

device_fourier_transform.o: device_fourier_transform.cu device_fourier_transform.h
	nvcc ${CFLAGS} -c -o device_fourier_transform.o device_fourier_transform.cu ${MPI_INCLUDES} -Xptxas -v

#device_fourier_transform.cubin: device_fourier_transform.cu device_fourier_transform.h
#	nvcc ${CFLAGS} -keep -cubin device_fourier_transform.cu ${MPI_INCLUDES} 

