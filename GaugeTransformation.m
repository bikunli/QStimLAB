function [gen] = GaugeTransformation(gen,R)
% This function performs gauge transformation (row transformation).

% INPUT:
% gen is a structure variable, with following properties:
% gen.Tableau is the tableau representation of generators
% gen.SignVector is a binary vector records the sign: (-1)^SignVector
% % the number of elements must be the same as the row number of gen.Tableau

% R is a binary r x r matrix that completes the row transformation, where r
% is the row number of gen.Tableau

% OUTPUT:
% gen: the transformed gen

% Version: v2.0, Date: 04/2024

% the transformation of the property SignVector is decomposed as two parts:
% the linear transformed part SignVector_temp and the non-linear part ExtraSign
SignVector_temp = mod(R*gen.SignVector, 2);

if mod(size(gen.Tableau,2),2) == 1
   error('gen.Tableau cannot have odd columns!!')
end
n = size(gen.Tableau,2)/2;
if size(R,1) ~= size(gen.Tableau,1) || size(R,1) ~= size(R,2)
    error('size(R,1) ~= size(gen.Tableau,1) || size(R,1) ~= size(R,2)');
end
if gfrank(R,2) ~= size(R,1)
    error('R does not have full row-rank')
end
if gfrank(R,2) ~= size(R,2)
    error('R does not have full column-rank')
end

for k = 1:size(R,1)
    % collect the rows to be operated, which is found in the k-th row of R.
    g_mat = gen.Tableau(logical(R(k,:)),:); 
    ExtraSign = 0; % phase factor: (-1)^ExtraSign
    m = size(g_mat,1); 
    % The follow steps mutiply the i-th row with the reference row,
    rr = g_mat(1,:); % Reference row: the first row
    for i = 2:m
        r4ins = g_mat(i,:); % row for inspection: will be multipied with the rr row.
        for j = 1:n
            % update the phase factor: (-1)^ph (the so called row-sum in PRA 70, 052328 (2004))
            ExtraSign = ExtraSign + ((rr(j)==1)*(rr(j+n)==1)*(r4ins(j+n) - r4ins(j)) ...
                + (rr(j)==1)*(rr(j+n)==0)*(r4ins(j+n)*(2*r4ins(j)-1)) ...
                + (rr(j)==0)*(rr(j+n)==1)*(r4ins(j)*(1-2*r4ins(j+n))))/2;
        end
        rr = mod(rr + r4ins,2); % row multiplication in tableau
    end

    if mod(ExtraSign,1) ~= 0 % ExtraSign must be an interger !
        error('ExtraSign = %2.2f is not an interger!',ExtraSign);    % which shouldn't happen if g_mat is legal
    else
        ExtraSign = round(ExtraSign);
    end
    gen.SignVector(k) = mod(SignVector_temp(k) + ExtraSign,2);
end
gen.Tableau = mod( R*gen.Tableau,2 );
end
