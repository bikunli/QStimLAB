function [stab_pt] = PartialTraceStab(stab, SubRegionIndex)
% This function performs partial trace on a stabilizer state defined by
% stab, the input stabilizer state can be either a pure or mixed state

% INPUT:
% stab: a stabilizer group, which has following properties:
% stab.Tableau is the tableau representation of generators
% stab.SignVector is a binary vector records the sign: (-1)^SignVector
% SubRegionIndex is the sub-region to be partial-traced, which is a subset of set {1,2,3,...,n},
% where n is the system size.

% OUTPUT
% stab_pt: The partial-traced stabilizer group

% % EXAMPLE:
% % paulistr = ['+XZIIZ';'-ZXZII';'+IZXZI';'+IIZXZ';'+ZIIZX'];
% % gen = PauliStrtoStab(paulistr);
% % [gen1_pt] = PartialTraceStab(gen, [1,2]);
% % GentoPaulistr(gen1_pt) % --> '+IIZXZ'

% Version: v2.0, Date: 04/2024

if ~isGenStabGroup(stab)
    error(['The input ',inputname(1),' is NOT a legit stabilizer group!']);
else
    stab = GetIndepStab(stab);
end
if isempty(SubRegionIndex)
    stab_pt = stab;
else
    n = size(stab.Tableau,2)/2;
    n_sub = numel(unique(SubRegionIndex));

    CompRegionIndex = setdiff(1:n,SubRegionIndex); % The complementary Region
    n_comp = numel(CompRegionIndex); % The size of the complementary region
    perm_ind = [SubRegionIndex, CompRegionIndex]; % new arrangement of sites

    remsystem_ind = n_sub + (1:n_comp); % the index of the remaining system of interest, after permutation
    stab.Tableau = stab.Tableau(:,[perm_ind,n+perm_ind]);
 
    [~,R] = rrefTableau(stab.Tableau);
    stab_rref = GaugeTransformation(stab,R);
    Tableau_perm_rref = int8(stab_rref.Tableau);
    Tableau_perm_rref_gf4 = Tableau_perm_rref(:,1:n) + 2*Tableau_perm_rref(:, n + (1:n));
    ind_remain = find(~any(Tableau_perm_rref_gf4(:,1:n_sub),2));

    stab_pt.Tableau = stab_rref.Tableau(ind_remain, [remsystem_ind, n+remsystem_ind]);
    stab_pt.SignVector = stab_rref.SignVector(ind_remain);
end
end