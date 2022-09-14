function [focus] = convertCellToMatrix(focus)
%CONVERTCELLTOMATRIX convert cell to matrix
%   Detailed explanation goes here
F = zeros(size(focus.Focused_vec{1},1),size(focus.Focused_vec{1},2),length(focus.Focused_vec));
Q = zeros(size(F));

for i = 1:length(focus.Focused_vec)
    F(:,:,i) = focus.Focused_vec{i};
    if(iscell(focus.Focus_eq))
        Q(:,:,i) = focus.Focus_eq{i};
    end
end
focus.Focused_vec = F;
if(iscell(focus.Focus_eq))
    focus.Focus_eq = Q;
end

end

