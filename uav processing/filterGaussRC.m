function [RC_filtr] = filterGaussRC(RC,dt,OSF,sigma)
%FILTERGAUSSRC filter with gaussian filter the range compressed matrix

L = 101;
x_fil = ((0:L-1)-floor(L/2))*dt/OSF;

g = 1/(sigma*sqrt(2*pi)) * exp(-0.5*(x_fil/sigma).^2); 
g = g/ max(g);

%% Filtering
RC_filtr =zeros(size(RC));
for n = 1:size(RC,2)
    RC_filtr(:,n) = conv(RC(:,n),g,'same');
end
%% Plot
figure
subplot(2,1,1)
imagesc([],[],abs(RC))
title('Range compressed original' )

subplot(2,1,2)
imagesc([],[],abs(RC_filtr))

title('Range compressed matrix side-lobe filtering - absolute value' )

end


