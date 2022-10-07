function [out] = gaussActivFunc(x,sigma)
%GAUSSACTIVFUNC gaussian activation function
%   [out] = gaussActivFunc(x,sigma)
out = exp(-0.5*((x)./sigma).^2);
end

