function M = Mincheol()

M = zeros(5,6);
M(:,1) = [100 50 30 10 5];
M(:,2) = M(:,1)./10;
M(:,3) = M(:,2)./(1.1*1000);
M(:,4) = (10-M(:,2))./(1000);
M(:,5) = M(:,3)./(M(:,3)+M(:,4));
M(:,6) = 1-M(:,5);

end