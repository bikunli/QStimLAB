function [stab] = SinglePauliMeasurement(p,s,i,stab)
% This function is a handy version of StabilizerMeasurement,
% which performs single qubit pauli measurement on a stabilizer state (pure or mixed)

% INPUT:
% paulistr: a pauli string, e.g. [1,0] or 'X' or 'x', it must have one row.
% s: the sign of the pauli string
% i: a scalar, which is the coordinate of the qubit to be measured

% stab: a stabilizer group, which has following properties:
% stab.Tableau is the tableau representation of generators
% stab.SignVector is a binary vector records the sign: (-1)^SignVector

% OUTPUT:
% stab: the updated stabilizer group

% Version: v2.0, Date: 04/2024

if ~isGenStabGroup(stab)
    error(['The input ',inputname(4),' is NOT a legit stabilizer group!']);
end

n = size(stab.Tableau,2)/2;
if any(~ismember(i,1:n))
    error(['The input ',inputname(3),' is out of the range 1:%d of the system!'],n);
end

switch p
    case {'x','X',[1,0]}
        p = 'X';
    case {'y','Y',[1,1]}
        p = 'Y';
    case {'z','Z',[0,1]}
        p = 'Z';
    otherwise
        error(['Unexpected input: ',inputname(1),'.']);
end

paulistr = repmat('I',[1,n]);
paulistr(i) = p;
stab = StabilizerMeausurement(paulistr, s, stab);
end