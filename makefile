CC         =  mpicc
CCFLAGS    =  -O3 -march=native
CCGFLAGS   =  -g
LIBS       =  -lmpi -lm

BINS= a3 

a3: a3_v2.c
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS)


clean:
	$(RM) $(BINS)
