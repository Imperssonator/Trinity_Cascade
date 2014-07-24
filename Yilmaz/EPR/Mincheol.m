function M = Mincheol()

M = zeros(5,6);
M(:,1) = [100 50 30 10 5]; %starting weight percent P3HT out of total polymer weight
M(:,2) = M(:,1)./10; %starting mg/mL P3HT
M(:,3) = M(:,2)./(1.1*1000); % starting vol frac P3HT
M(:,4) = (10-M(:,2))./(1000); % starting vol frac PS
M(:,5) = M(:,3)./(M(:,3)+M(:,4))*0.5; % formula for final vol frac P3HT.. ending at 0.5 vol frac solvent
M(:,6) = 0.5-M(:,5);

end