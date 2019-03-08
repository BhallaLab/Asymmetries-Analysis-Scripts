clearvars
clc

grid29 = zeros(841,4);
grid29(:,1)=linspace(1,841,841);
grid29(:,2)=linspace(29,29,841);
grid29(:,3)=linspace(29,29,841);
grid29(:,4)= randperm(841,841);

dlmwrite('Random_grid_29x_coordinates.txt',grid29,'delimiter',' ','newline','pc')