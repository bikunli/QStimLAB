function [x] = isStabSubgroup(gen1,gen2)
% this function determines whether G1 is a subgroup of G2
% where G1 and G2 are stabilizer groups determined by gen1 and gen2
% INPUT:
% gen1 and gen2: struct variables with properties Tableau and SignVector

% OUTPUT:
% x: true/false
% this function returns error is either gen1 or gen2 are illegal input

% Version: v2.0, Date: 04/2024

x = false;
gen1 = GetIndepStab(gen1);
gen2 = GetIndepStab(gen2);
if isPauliStrCommutative(gen1.Tableau, gen2.Tableau) && size(gen1.Tableau,1) <= size(gen2.Tableau,1) % it returns error if gen1 and gen2 act on systems of different sizes
    gen12.Tableau = [gen1.Tableau; gen2.Tableau]; % the union set of gen1 and gen2
    gen12.SignVector = [gen1.SignVector; gen2.SignVector];
    if isGenStabGroup(gen12)
        gen12_indep = GetIndepStab(gen12);
        if size(gen12_indep.Tableau,1) ==  size(gen2.Tableau,1)
        % in this case, gen12 and gen2 generates the same stabilizer group
        % which implies gen1 generates a subgroup of what gen2 generates
            x = true;
        end
    end
else
    
end

end


