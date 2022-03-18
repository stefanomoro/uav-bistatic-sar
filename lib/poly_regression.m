function [fit,coeff]= poly_regression(x,y,order)
    X = ones(length(x),order + 1);
    for i = 0:order
        X(:,i+1) = x.^(i);
    end
    %compute coefficients
    coeff =  inv(X'*X)*X'* y;

    %computes actual fit
    fit = zeros(size(y));
    for i = 0 : order
        fit = fit + (coeff(i+1).*x.^(i));
    end
end