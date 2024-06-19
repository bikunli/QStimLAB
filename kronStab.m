function [stab] = kronStab(stab1, stab2)
% Thus function finds the kronecker product (tensor product) of two
% stabilizer states.

% INPUT:
% stab1, stab2: two stabilizer generator set that represent two stabilizer states.

% stab: a stabilizer group, which has following properties:
% stab.Tableau is the tableau representation of generators
% stab.SignVector is a binary vector records the sign: (-1)^SignVector

% OUTPUT:
% stab: the product state 'stab1 x stab2'

if ~isGenStabGroup(stab1)
    error(['the input ',inputname(1),' does not represent a stabilizer group!'])
end
if ~isGenStabGroup(stab2)
    error(['the input ',inputname(2),' does not represent a stabilizer group!'])
end

tab1 = stab1.Tableau;
n_1 = size(tab1,2)/2;
tab1_X = tab1(:,1:n_1);
tab1_Z = tab1(:, n_1 + (1:n_1));

tab2 = stab2.Tableau;
n_2 = size(tab2,2)/2;
tab2_X = tab2(:,1:n_2);
tab2_Z = tab2(:, n_2 + (1:n_2));

stab.Tableau = [blkdiag(tab1_X,tab2_X),blkdiag(tab1_Z,tab2_Z)];
stab.SignVector = [stab1.SignVector; stab2.SignVector];

end