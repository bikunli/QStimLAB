function [stab] = BellMeausurement(sites, basis, stab)
% this function performs Bell state measurment on a pair of qubits of a stabilizer group stab

% INPUT:
% sites: [s1, s2], two sites to be measured
% basis: [b1, b2], which are the signs of 'XX' and 'ZZ'

% stab: a stabilizer group, which has following properties:
% stab.Tableau is the tableau representation of generators
% stab.SignVector is a binary vector records the sign: (-1)^SignVector

% OUTPUT:
% stab: the updated stabilizer group

% Version: v2.0, Date: 04/2024

if ~isGenStabGroup(stab)
    error(['The input ',inputname(4),' is NOT a legit stabilizer group!']);
end
if numel(sites) ~= 2
    error('sites has to be an array with two elements.')
end

n = size(stab.Tableau,2)/2;
if any(~ismember(sites,1:n))
    error(['The input ',inputname(3),' is out of the range 1:%d of the system!'],n);
end

gen_Bell_str = repmat('I',[2,n]);
gen_Bell_str(1,sites) = 'XX';
gen_Bell_str(2,sites) = 'ZZ';

for i = 1:2
    stab = StabilizerMeausurement(gen_Bell_str(i,:),basis(i),stab);
end
end