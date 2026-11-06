%% Report temporal resolution of the time-frequency wavelets
% Computes the temporal full width at half maximum (FWHM), in seconds, for
% the logarithmically spaced frequencies and cycle counts used in the EEG
% time-frequency analyses. The result is returned in vector `h`.
% Author: Davide Gheza


min_freq = 1.5;
max_freq = 40;
num_frex = 20;   
wavelet_cycles_lims = [2,10];  % n = number of cycles 

frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
n_cycles = logspace(log10(wavelet_cycles_lims(1)),log10(wavelet_cycles_lims(2)),num_frex);

h = n_cycles*(sqrt(2*log(2))) ./ (pi*frequencies); 
