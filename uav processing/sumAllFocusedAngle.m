function [A] = sumAllFocusedAngle(focus_vec)
%SUMALLFOCUSEDANGLE sum together not coherently all the focused image
%   [A] = sumAllFocusedAngle(focus_vec)
A = zeros(size(focus_vec{1}));

% for i = 1:length(focus_vec)
for i = 6:10
    temp = abs(focus_vec{i}).^2;

    A = A + temp;
end
end

