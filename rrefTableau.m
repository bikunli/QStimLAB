function [Tableau_echelon,R] = rrefTableau(Tableau)
% The echelon form Algorithm introduced in paper: New J. Phys. 7, 170 (2005)

% INPUT:
% Tableau: tableau matrix representations of pauli strings

% OUTPUT:
% Tableau_echelon: the transformed tableau matrix under echelon gauge
% R: a binary matrix that transforms rows of tableau: Tabluea_echelon = R * Tableau

% Version: v2.0, Date: 04/2024

[r,c] = size(Tableau);

if mod(c,2) == 1
    error(['The input ',inputname(1),' does NOT have even number of columns']);
end

tab_x = Tableau(:,1:(c/2));
tab_z = Tableau(:,(c/2)+(1:(c/2)));

% interleaving the tab_x and tab_z to get a alternative column arrangement
Tableau_alt = zeros(size(Tableau));
Tableau_alt(:,1:2:c) = tab_x;
Tableau_alt(:,2:2:c) = tab_z;

Tableau_alt_rref = rrefgf([Tableau_alt, eye(r)],2);

tab_x = Tableau_alt_rref(:,1:2:c);
tab_z = Tableau_alt_rref(:,2:2:c);
R = logical(Tableau_alt_rref(:,c + (1:r)));

% Return to conventional Tableau colomn order
Tableau_echelon = logical([tab_x,tab_z]);

end