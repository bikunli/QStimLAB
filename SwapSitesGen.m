function gen = SwapSitesGen(gen, perm, inv)
% This function performs SWAP operation on a set of Pauli strings

% INPUT:
% gen is a structure variable, with following properties:
% gen.Tableau is the tableau representation of generators
% gen.SignVector is a binary vector records the sign: (-1)^SignVector
% perm: a permutation vector that permutes subsystems
% inv: 0 by default, if inv = 1, performs the inverted permutation

if nargin == 2
    inv = false;
end
if mod(size(gen.Tableau, 2),2) == 1
    error(['The input ',inputname(1),'.Tableau must have even columns!']);
end
n = size(gen.Tableau, 2)/2;
permmat = speye(n);

switch inv
    case true
        Q = blkdiag(permmat(:,perm),permmat(:,perm));
    case false
        Q = blkdiag(permmat(perm,:),permmat(perm,:));
    otherwise
        error(['Unexpected input of ',inputname(3),'!']);
end
gen.Tableau = mod(gen.Tableau * Q,2);
end