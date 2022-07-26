function [out] = gaussActivFunc(x,sigma)
%GAUSSACTIVFUNC gaussian activation function
%   [out] = gaussActivFunc(x,sigma)
    out = 1/(sigma*sqrt(2*pi)) * exp(-0.5*((x)./sigma).^2); 
    out = out./max(out);
end

