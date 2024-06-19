function [x] = isGenStabGroup(gen)
% This funciton determines whether gen generates a stabilizer group
% % % A stabilizer group is an Abelian Pauli subgroup, that does not contain -I.

% INPUT
% gen: a struct variable that has properties Tableau and SignVector
% % stabilizers in gen is not necessarily independent
% % Tableau is a binary matrix with 2*n columns and m rows
% % SignVector is an mx1 binary vector

% OUTPUT: 
% x is either true or false

% Version: v2.0, Date: 04/2024

x = false;
if mod(size(gen.Tableau,2),2)==1
    error('gen.Tableau cannot have odd columns!!')
end
if ~isequal(gen.Tableau, mod(gen.Tableau,2))
    error('gen.Tableau is not a binary matrix!')
end
if ~isequal(gen.SignVector, mod(gen.SignVector,2))
    error('gen.SignVector is not a binary vector!')
end

if ~isPauliStrCommutative(gen.Tableau, gen.Tableau)
    warning('gen is not a commutative set!')
elseif exist_minus_I(gen)
    warning('the generator generates operator -I !')
else
    x = true;
end
end

function [x] = exist_minus_I(gen)
% this function determines whether the generator set gen generates -I
% output: x is either true or false

x = false;
[~,R] = rrefTableau(gen.Tableau);
gen = GaugeTransformation(gen, R);
ind_I = ~any(gen.Tableau,2); % get the indices of identity operators
if ~isempty(ind_I)
    x = any(gen.SignVector(ind_I)); % false if all I are +I
end
end