function [stab] = allzerostabstate(n)
% This function creates a n-qubit all-zero stabilizer state

% INPUT:
% n: an integer for the size of the system
% OUTPUT:
% stab: representa a pure stabilizer state, stabilized by +Z_i, for i = 1,...,n

% Version: v2.0, Date: 04/2024

stab.Tableau = sparse([zeros(n),eye(n)]);
stab.SignVector = sparse(zeros(n,1));

end