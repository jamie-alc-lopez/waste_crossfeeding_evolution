function [C,tlags,fs,tlin] = spline_autocorr(t,f)

N =5*numel(t);
tlin = linspace(t(1), t(end), N);
fs = interp1(t,f,tlin,'spline');
C = xcorr(fs,fs,'unbiased');
tlags = [fliplr(-tlin(2:end)+t(1)),tlin-t(1)];

end