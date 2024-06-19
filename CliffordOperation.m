function [gen] = CliffordOperation(gen,gate_str,site)
% This function operates Generators by an elementary clifford gate

% INPUT:
% gen is a structure variable, with following properties:
% gen.Tableau is the tableau representation of generators
% gen.SignVector is a binary vector records the sign: (-1)^SignVector
% % the length of gen.SignVector must be the same as the row number of gen.Tableau
% % gen needs NOT to be a commutative set of Pauli operators
% gate_str can be 'I','H','P','PH','HP','HPH','X','Y','Z','CNOT',.......
% * but 'SWAP' is not included
% (more gate_str see the dictionary: OtherLC)
% Other operation should be realized repeatedly using this function
% site is a scalar for single qubit gate
% site is a vector [i,j] for controlled-gate CNOT_{i,j}

% OUTPUT:
% gen: the transformed gen

% Version v2.0,  Date: 04/2024

[n_r,n_c] = size(gen.Tableau);
if mod(n_c,2) == 1
    error('n_c cannot be an odd number!')
end
if n_r>(n_c/2)
    error('n_r>(n_c/2)!')
end
n = n_c/2;

OtherLC = {'HP','PH','HPH', ...
    'HX','PX','HPX','PHX','HPHX', ...
    'HY','PY','HPY','PHY','HPHY', ...
    'HZ','PZ','HPZ','PHZ','HPHZ'}; % other single-qubit Clifford gates
switch gate_str
    case {'I',''}
    case 'H'
        gen.Tableau = mod(gen.Tableau * Qmatrix(n, gate_str ,site),2);
        gen.SignVector = mod(gen.SignVector + gen.Tableau(:,site).*gen.Tableau(:,n+site),2);
    case 'P'
        gen.Tableau = mod(gen.Tableau * Qmatrix(n, gate_str ,site),2);
        gen.SignVector = mod(gen.SignVector + gen.Tableau(:,site).*gen.Tableau(:,n+site),2);
    case 'X'
        gen.SignVector = mod(gen.SignVector + gen.Tableau(:,n+site),2);
    case 'Y'
        gen.SignVector = mod(gen.SignVector + gen.Tableau(:,site) + gen.Tableau(:,n+site),2);
    case 'Z'
        gen.SignVector = mod(gen.SignVector + gen.Tableau(:,site),2);
    case OtherLC
        for i_str = 1:length(gate_str)
            gen = CliffordOperation(gen,gate_str(i_str),site);
        end
    case {'CNOT','CX'}
        gen.Tableau = mod(gen.Tableau * Qmatrix(n,'CNOT',site(1), site(2)),2);
        gen.SignVector = mod(gen.SignVector ...
            + gen.Tableau(:,site(1)).*gen.Tableau(:,n+site(2)) ...
            .*(gen.Tableau(:,site(2)) + gen.Tableau(:,n+site(1)) + 1),2);
    case {'CZ'}
        gen = CliffordOperation(gen,'H',site(2));
        gen = CliffordOperation(gen,'CNOT',site);
        gen = CliffordOperation(gen,'H',site(2));
    otherwise
        fprintf(1,['gate_str = %s .\n'],gate_str);
        error('Unexpected string of variable ''gate_str''  ');
end
end

function [Q]=Qmatrix(N,type,i,j)
% Q matrix is a matrix that right-multiply on the tableau matrix
if strcmp(type,'CZ') == 1 || strcmp(type,'CNOT') == 1 || strcmp(type,'CX') == 1
    if nargin < 4
        error('Insufficient input arguments!');
    end
end
Q = speye(2*N);

switch type
    case 'I'
        Q([i,N+i],[i,N+i]) = eye(2);
    case 'H'
        Q([i,N+i],[i,N+i]) = [0,1;1,0];
    case 'P'
        Q([i,N+i],[i,N+i]) = [1,1;0,1];
    case 'HP'
        Q([i,N+i],[i,N+i]) = [0,1;1,1];
    case 'PH'
        Q([i,N+i],[i,N+i]) = [1,1;1,0];
    case 'HPH'
        Q([i,N+i],[i,N+i]) = [1,0;1,1];
    case {'CNOT','CX'}
        Q(i,j) = 1;
        Q(N+j,N+i) = 1;
    case 'CZ'
        Q(j,[N+i]) = 1;
        Q(i,[N+j]) = 1;
    otherwise
        error('Unexpected string of variable ''type''! ');
end

end
