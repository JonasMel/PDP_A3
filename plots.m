% Array size = 100000000, Number of processes = 1, Lowest time (of 1 runs) to sort = 10.587318 seconds. 
% Array size = 100000000, Number of processes = 2, Lowest time (of 1 runs) to sort = 6.151277 seconds. 
% Array size = 100000000, Number of processes = 4, Lowest time (of 1 runs) to sort = 3.963363 seconds. 
% Array size = 100000000, Number of processes = 8, Lowest time (of 1 runs) to sort = 3.066966 seconds. 
% Array size = 100000000, Number of processes = 16, Lowest time (of 1 runs) to sort = 2.255670 seconds. 
% Array size = 100000000, Number of processes = 32, Lowest time (of 1 runs) to sort = 2.083552 seconds. 
% Array size = 100000000, Number of processes = 64, Lowest time (of 1 runs) to sort = 2.238509 seconds.
% Array size = 1000000000, Number of processes = 1, Lowest time (of 1 runs) to sort = 108.506589 seconds. 
% Array size = 1000000000, Number of processes = 2, Lowest time (of 1 runs) to sort = 63.433258 seconds. 
% Array size = 1000000000, Number of processes = 4, Lowest time (of 1 runs) to sort = 42.103183 seconds. 
% Array size = 1000000000, Number of processes = 8, Lowest time (of 1 runs) to sort = 29.133287 seconds. 
% Array size = 1000000000, Number of processes = 16, Lowest time (of 1 runs) to sort = 23.660104 seconds. 
% Array size = 1000000000, Number of processes = 32, Lowest time (of 1 runs) to sort = 18.765329 seconds. 
% Array size = 1000000000, Number of processes = 64, Lowest time (of 1 runs) to sort = 15.405089 seconds. 

times_100mil = [10.587318 6.151277 3.963363 3.066966 2.255670 2.083552 2.238509];
times_1bil = [108.506589 63.433258 42.103183 29.133287 23.660104 18.765329 15.405089];
processes = [1 2 4 8 16 32 64];

plot(processes, processes, processes, times_100mil(1)./times_100mil, processes, times_1bil(1)./times_1bil);
xlabel('Processes');
ylabel('speed-up');