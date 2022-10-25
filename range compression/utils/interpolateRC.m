function [RC_OS] = interpolateRC(RC,OSF,norm_B)
%INTERPOLATERC interpolate RC wit OSF 
 

t = 0:size(RC,1) -1;
tf = 0:1/OSF:(size(RC,1) -1);

tf_men_t =tf(:) - t;
W = sinc(norm_B*tf_men_t);

% interpolate whole RC signal
RC_OS = W*RC;


end

