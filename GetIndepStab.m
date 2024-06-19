function [gen_trimmed] = GetIndepStab(gen)
% it is possible that the Pauli operators in gen are not independent
% this function trims the generator set gen as an independent set

% INPUT
% gen: a struct variable which has properties Tableau and SignVector, it must be a legit stabilizer group

% OUTPUT:
% gen_trimmed: equivalent to gen, but without redundancy

% Version: v2.0, Date: 04/2024

if ~isGenStabGroup(gen)
    error(['The input ',inputname(1),' does NOT generate a stabilizer group!']);
end

ind = find(~all(gen.Tableau==0, 2)); % find all the identity stabilizers
gen.Tableau = gen.Tableau(ind,:); % keep those non-identity stabilizers
gen.SignVector = gen.SignVector(ind);


% trivially remove the repetition of stabilizers
[~, ind_u] = unique(gen.Tableau,"rows"); 
if ~all(ismember(1:size(gen.Tableau,1),ind_u))
    gen.Tableau = gen.Tableau(ind_u,:);
    gen.SignVector = gen.SignVector(ind_u);
end

% if the remaining tableau is still not full (row) rank
if ~((gfrank(gen.Tableau,2) == size(gen.Tableau,1)))
    [r,c] = size(gen.Tableau);
    tab_ex = [gen.Tableau, eye(r)]; % the extended tableau
    tab_ex_rref = rrefgf(tab_ex,2);
    R = tab_ex_rref(:,c+(1:r));
    gen = GaugeTransformation(gen, R); % the redundant stabilizers becomes +I
    gen = GetIndepStab(gen); % self-iteration
end
gen_trimmed = gen;
end