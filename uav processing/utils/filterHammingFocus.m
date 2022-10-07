function [A] = filterHammingFocus(F,L)
%HAMMINGFILTERFOCUS filter the focused image with hamming filter
%   [A] = filterHammingFocus(F,L)
filt = hamming(L);
filt = filt ./L;
F = abs(F);
A = conv2(F,filt,'same');
A = conv2(A,filt.','same');
end

