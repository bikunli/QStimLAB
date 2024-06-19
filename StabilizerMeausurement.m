function [stab] = StabilizerMeausurement(ps, s, stab)
% this function performs general stabilizer measurment on a stabilizer group stab
% if ps commutes with the stabilizer group, but the sign s is incompatible,
% then the output stab will be the same as stab

% INPUT:
% ps: a pauli string, e.g. [1,1,0,0] or 'XX', it must have one row
% s: the sign of the pauli string

% stab: a stabilizer group, which has following properties:
% stab.Tableau is the tableau representation of generators
% stab.SignVector is a binary vector records the sign: (-1)^SignVector

% OUTPUT:
% stab: the updated stabilizer group

% Version: v2.0, Date: 04/2024

if ~isGenStabGroup(stab)
    error(['The input ',inputname(1),' is NOT a legit stabilizer group!']);
end
if size(ps,1) ~= 1
    error(['The input ',inputname(1),' can only have one row.']);
end
if all(ismember('+-',ps))
    error(['Please remove +/- in',inputname(1),'.'])
end
if all(ismember(ps,'IXYZ')==1)
    ps = paulistr2tableau(ps);
end
switch s
    case '+'
        s = 0;
    case '-'
        s = 1;
    case {0,1}
    otherwise
        error(['Unexpected input: ',inputname(2),' !'])
end

n = size(stab.Tableau,2)/2;

if isPauliStrCommutative(ps, stab.Tableau)
    if gfrank([stab.Tableau; ps],2) > gfrank([stab.Tableau],2)
        % in this case, stab is commutative with paulistr, and stab does not
        % contain either +paulistr or -paulistr.
        stab.Tableau = [stab.Tableau; ps];
        stab.SignVector = [stab.SignVector; s];
    else
        % in this case, stab already contains a stabilizer with opposite
        % sign_ps, so stab is not updated
    end
else
    % in this case, paulistr is not commutative with stab
    P = [zeros(n), eye(n); eye(n), zeros(n)];
    R = eye(size(stab.Tableau,1));
    anticomm_list = mod(stab.Tableau * P * ps',2);
    ind = find(anticomm_list,1,'first');
    R(:,ind) = anticomm_list;
    stab = GaugeTransformation(stab, R);
    stab.Tableau(ind,:) = ps;
    stab.SignVector(ind) = s;
end
end


function [tableau] = paulistr2tableau(pauli_str)

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