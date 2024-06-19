function x = isEqualStabGroup(stab1, stab2)
% this function determines whether stabilizer group stab1 is equivalent to stab2

% INPUT: stab1, stab2: two sets of Pauli strings, which satisfy:
%    true = isGenStabGroup(stab1)
%    true = isGenStabGroup(stab2)
% OUTPUT: x = true/false

if isStabSubgroup(stab1, stab2) && isStabSubgroup(stab2, stab1)
    x = true;
else
    x = false;
end
end