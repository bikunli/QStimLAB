function [h] = SubRegionEntropy(Tableau, SubRegionIndex)
% This function computes the Renyi entropy of a subsystem (SubRegion)
% the total system can be a pure or mixed state

% INPUT:
% Tableau: the tableau matrix representation of a set of Pauli operators
% SubRegionIndex is a subset of set {1,2,3,...,n}, where n is the system size.
% Tableau is a binary matrix of even columns, which may not be full rank.
% example: [1,1, 0,0; 0,0, 1,1] represents two pauli strings 'XX' and 'ZZ'
% this function returns error if Tableau does not represents a set of
% commutative pauli operators.

% OUTPUT:
% h: the entropy of subsytem

% EXAMPLE:
% paulistr = ['XZIIZ';'ZXZII';'IZXZI';'IIZXZ';'ZIIZX'];
% gen = paulistr2stab(paulistr);
% SubRegionEntropy(gen1.Tableau,[3,4]) % --> 2

% Version: v2.0, Date: 04/2024

if isstruct(Tableau)
    Tableau = Tableau.Tableau;
end
if ~isPauliStrCommutative(Tableau,Tableau)
    error(['The input ',inputname(1),' does NOT represent a stabilizer group!'])
end


if isempty(Tableau) || all(Tableau == 0,'all') % if Tableau represents a completely mixed state
    h = numel(SubRegionIndex);
else
    n = size(Tableau,2)/2;
    P = [zeros(n),eye(n); eye(n),zeros(n)];

    if any(~ismember(SubRegionIndex,1:n))
        error(['The input ',inputname(2),' is out of the range of the system!']);
    end

    if ~all(mod(Tableau * P * Tableau.',2)==0,"all")
        error(['The input ',inputname(1),' does NOT represent a set of commutative pauli operators!'])
    end

    CompRegionIndex = setdiff(1:n,SubRegionIndex); % The complementary Region
    CompRegionSize = numel(CompRegionIndex); % The length of the complementary region
    perm_ind = [CompRegionIndex,SubRegionIndex]; % new arrangement of sites
    Tableau_perm = Tableau(:,[perm_ind,n+perm_ind]);

    Tableau_perm_rref = int8(rrefTableau(Tableau_perm));
    Tableau_perm_rref( ~any(Tableau_perm_rref,2), : ) = []; % remove zeros rows

    if isempty(Tableau_perm_rref)
        numstab = 0;
        numstab_c = 0;
    else
        numstab = size(Tableau_perm_rref,1);
        Tableau_perm_rref_gf4 = Tableau_perm_rref(:,1:n) + 2*Tableau_perm_rref(:, n + (1:n));
        numstab_c = sum(any(Tableau_perm_rref_gf4(:,1:CompRegionSize),2));
    end

    h = (n - CompRegionSize) - (numstab - numstab_c);
end
end