function out = PFDIFF(PF)
%% sum squared phase fraction difference

out = sum((PF(:,1)-PF(:,2)).^2);

end