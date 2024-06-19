function [Tableau_gauged,Tableau_lclean, Hadamard_v, PhaseGate_v] = GraphGauge(Tableau)
% This function convert any full rank tableau into an equivalent graph
% state Tableau: Tableau_gauged
% The second output Tableau_lclean is a intermediate Tableau for debugging.
% The last two output: Hadamard_v and PhaseGate_v are logical vector with
% length L, which indicates Hadamard and PhaseGate needed to convert
% original Tableau to Tableau_gauged:
% (HP)^\dagger.Tableau.(HP) ~ Tableau_gauged
L = size(Tableau,1);
if gfrank(Tableau) ~= L
    warning('The input Tableau is not full ranked!');
end
if size(Tableau,1)*2 ~= size(Tableau,2)
    warning('The inpute Tableau has wrong size!');
end

Hadamard_v = zeros(L,1);
PhaseGate_v = zeros(L,1);
for i_r = 1:L
    rpiv = find(Tableau(i_r:L,i_r) == 1,1,'first');
    rpiv = rpiv + i_r - 1;
    if isempty(rpiv) == 1
        Tableau(:,[i_r,i_r+L]) = Tableau(:,[i_r+L,i_r]); % the pre-cleaning step
        Hadamard_v(i_r) = 1;
        rpiv = find(Tableau(i_r:L,i_r) == 1,1);
        rpiv = rpiv + i_r - 1;
    end
    if i_r < L
        if rpiv  ~= i_r
            Tableau([i_r,rpiv],:) = Tableau([rpiv,i_r],:);
        end
        
        ind2belim = find(Tableau((i_r+1):L,i_r) == 1);
        ind2belim = ind2belim + i_r;
        Tableau(ind2belim,:) = mod(Tableau(ind2belim,:)+Tableau(i_r,:),2);
    end
end
Tableau_lclean = Tableau;

for i_r = L:-1:2
    ind2belim = find(Tableau(1:(i_r-1),i_r) == 1);
    Tableau(ind2belim,:) = mod(Tableau(ind2belim,:)+Tableau(i_r,:),2);
end

Tableau_gauged = Tableau;
ind2belim = find(diag(Tableau(:,L+(1:L))) == 1 );
PhaseGate_v(ind2belim) = 1;
Tableau_gauged(ind2belim,L+(1:L)) = mod(Tableau_gauged(ind2belim,1:L) + Tableau_gauged(ind2belim,L+(1:L)),2);

Hadamard_v = logical(Hadamard_v);
PhaseGate_v = logical(PhaseGate_v);
end