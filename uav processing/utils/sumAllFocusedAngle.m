function [A] = sumAllFocusedAngle(focus_vec,idxs)
%SUMALLFOCUSEDANGLE sum together not coherently all the squint image
%   [A] = sumAllFocusedAngle(focus_vec,idxs)
A = zeros(size(focus_vec(:,:,1)));


figure
for i = idxs
    temp = abs(focus_vec(:,:,i));
%     imagesc(20*log10(temp.')),colorbar
%     caxis([140 200])
%     colormap jet;
%     pause(1)
    A = A + temp;
end
end

