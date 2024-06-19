function x = isPauliStrCommutative(ps1, ps2)
% this function determines whether two sets of Pauli strings are commutative
% ps1 and ps2 are tableau matrices of same number of columns

% INPUT:
% ps1 and ps2: tableau representations of pauli strings
% % Example of input:
% % (1) ps1 = [0,1,1,0, 0,0,1,1] represents the pauli string 'IXYZ'
% % (2) ps1 = [1,1, 0,0; 0,0, 1,1] represents two pauli strings 'XX' and 'ZZ'
% % (3) ps2 = 'XXZZII'

% OUTPUT:
% x is either true or false

% v1.0. Date: 04/2024

if nargin == 1
    ps2 = ps1;
end

% if ps1 or ps2 are char variables, then convert ps1 and ps2 to tableau representations
str_list = ['I','Z';'X','Y'];
if all(ismember(ps1,str_list))
    ps1 = PauliStrtoTableau(ps1);
end
if all(ismember(ps2,str_list))
    ps2 = PauliStrtoTableau(ps2);
end

% rule out the case that either ps1 or ps2 is non-binary
if ~isequal(ps1, mod(ps1,2))
    error('ps1 is not a binary matrix!')
elseif ~isequal(ps2, mod(ps2,2))
    error('ps2 is not a binary matrix!')
end

if mod(size(ps1,2),2) == 1
    error(['The input ',inputname(1),' does NOT have even number of columns!'])
elseif mod(size(ps2,2),2) == 1
    error(['The input ',inputname(2),' does NOT have even number of columns!'])
end

if size(ps1,2) ~= size(ps2,2)
    error(['The inputs ',inputname(1),' and ',inputname(2),' do NOT have the same number of columns!']);
else
    n = size(ps1,2)/2;
    P = sparse([zeros(n),eye(n); eye(n), zeros(n)]);
    x = all(mod(ps1 * P * ps2.',2) == 0,'all');
end
end