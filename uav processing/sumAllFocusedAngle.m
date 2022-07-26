function [A] = sumAllFocusedAngle(focus_vec,idxs)
%SUMALLFOCUSEDANGLE sum together not coherently all the squint image
%   [A] = sumAllFocusedAngle(focus_vec,idxs)
A = zeros(size(focus_vec{1}));

% for i = 1:length(focus_vec)
for i = idxs
    temp = abs(focus_vec{i});

    A = A + temp;
end
end

