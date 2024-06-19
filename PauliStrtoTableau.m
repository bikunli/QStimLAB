function [tableau] = PauliStrtoTableau(pauli_str)

% this function converts the input Pauli string pauli_str into a tableau
% tableau is a binary matrix of 2*n columns

str_list = ['I','Z';'X','Y'];
if ~all(ismember(pauli_str,str_list))
    error('Unexpected symbol: All elements of pauli_str should be a member of {I,X,Y,Z}!')
end
[rownum,colnum] = size(pauli_str);
tableau = zeros(rownum,2*colnum);
for i_r = 1:rownum
    for i_c = 1:colnum
        [r,c] = find(str_list == pauli_str(i_r,i_c), colnum);
        tableau(i_r,[i_c,colnum+i_c]) = [r-1,c-1];
    end
end
end