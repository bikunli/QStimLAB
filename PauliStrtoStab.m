function gen = PauliStrtoStab(paulistr_signed)
%  This function converts a set of signed Pauli string to a struct
%  variable gen, which has properties Tableau and SignVector

% INPUT:
%  The legit input of paulistr can be:
%  paulistr_signed = '+IXYZ' or ['+XX';'-ZZ'] or ['XZZ';'ZXZ';'ZZX'];
%  paulistr_signed must be compatible to generate a stabilizer group, otherwise
%  this function returns error

% OUTPUT:
% gen:  a struct variable, which has properties Tableau and SignVector

% EXAMPLE:
% paulistr_signed = ['+XX';'-ZZ'];
% gen3 = PauliStrtoStab(paulistr_signed)

% Version: v2.0, Date: 04/2024

signset = '+-';
letterset = 'IXYZ';
if all(ismember(paulistr_signed,letterset),'all')
    paulistr_sign = repmat('+',[size(paulistr_signed,1),1]);
    paulistr_letter = paulistr_signed;
elseif all(ismember(paulistr_signed(:,2:end),letterset),'all') ...
        && all(ismember(paulistr_signed(:,1),signset),'all')
    paulistr_sign = paulistr_signed(:,1);
    paulistr_letter = paulistr_signed(:,2:end);
else
    error(['Unexpected input is found in ',inputname(1),'!'])
end

[~,ind] = ismember(paulistr_sign, signset);
gen.SignVector = ind(:)-1;
gen.Tableau = PauliStrtoTableau(paulistr_letter);
if ~isGenStabGroup(gen)
    warning('paulistr does NOT generate a stabilizer group!');
end
end


