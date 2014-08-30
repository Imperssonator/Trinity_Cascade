function out = Gij(species,VF)

out = zeros(3);

out(1,2) = get_X(species{1},species{2});
out(1,3) = get_X(species{1},species{3});
out(2,3) = get_X(species{2},species{3});

end

function out = get_X(s1,s2)
%% Get X
% s1 and s2 are strings that correspond to species. Get X will search
% a cell array to find which position each species occupies. Then it will
% use those indices to grab the binary parameter from a matrix located at
% i,j where i and j are the cell array positions of the species strings...

SPLIST = struct('P3HT',1, ...
                'PS',2, ...
                'CHCl3',3, ...
                'hexane',4, ...
                'DCB',5, ...
                'acetone',6 ...
                );

i1 = SPLIST.(s1);
i2 = SPLIST.(s2);
BINPARAMS = sparse(zeros(6));
BINPARAMS(1,2) = 0.48;
BINPARAMS(1,3) = 0.99;
BINPARAMS(1,4) = 1.4;
BINPARAMS(1,5) = 0.5; %not right
BINPARAMS(2,3) = 0.39;
BINPARAMS(3,4) = 0.55;

BINPARAMS = BINPARAMS+BINPARAMS'; % make it so i,j gives the same value as j,i

out = BINPARAMS(i1,i2);

end