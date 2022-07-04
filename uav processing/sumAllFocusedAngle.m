function [A] = sumAllFocusedAngle(focus_vec)
%SUMALLFOCUSEDANGLE sum together not coherently all the focused image
%   Detailed explanation goes here
A = zeros(size(focus_vec{1}));
for i = 1:size(focus_vec,1)
    A = A + abs(focus_vec{i});
end
end

