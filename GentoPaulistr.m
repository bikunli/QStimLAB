function paulistr = GentoPaulistr(gen)
% This function converts a set of generators to readable Pauli strings
%
% INPUT:
% gen: a set of Hermitian Pauli operators, which has following properties:
% gen.Tableau: the tableau matrix representation of pauli strings, which is
% not necessary commutative
% gen.SignVector is a binary vector records the sign: (-1)^SignVector

% OUTPUT:
% paulistr: a char array, each row explicitly represents an Hermitian Pauli operator

% EXAMPLE:
% paulistr = ['XZIIZ';'ZXZII';'IZXZI';'IIZXZ';'ZIIZX'];
% gen = GentoPaulistr(paulistr);
% GentoPaulistr(gen)

% Version: v2.0, Date: 04/2024

[r,c] = size(gen.Tableau);
if mod(c,2) == 1
    error([inputname(1),'.Tableau cannot have odd number of columns!']);
end
n = c/2;
paulistr_sign = repmat('?',[r,1]);
paulistr_letter = repmat('I',[r,n]);

sign_str = '+-';
str_ls = {'I','Z';'X','Y'};

for i_r = 1:r
    paulistr_sign(i_r) = sign_str(gen.SignVector(i_r)+1);
    for i_c = 1:n
        paulistr_letter(i_r,i_c) = str_ls{gen.Tableau(i_r,i_c)+1,gen.Tableau(i_r,n+i_c)+1};
    end
end

paulistr = [paulistr_sign, paulistr_letter];

end