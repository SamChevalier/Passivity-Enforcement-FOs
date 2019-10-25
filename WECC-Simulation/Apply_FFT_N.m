function [f_vec,y_fft] = Apply_FFT_N(signal,tstep,N)
% INPUT a signal and its associated parameters, and OUTPUT the fft data.

% Ensure that N is odd
if mod(N,2) == 0
    N = N-1;
end

% Now, take the N point FFT
y     = fft(signal/N,N);
fs    = 1/tstep;
f_vec = fs*(0:((N-1)/2))/N;

% Analyze Signal
y_fft = y(1:((N+1)/2));

% Double all but DC
y_fft(2:end) = 2*y_fft(2:end);
end