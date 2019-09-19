clearvars
clc

grid59 = zeros(3841,4);
grid59(:,1)=linspace(1,3841,3841);
grid59(:,2)=linspace(59,59,3841);
grid59(:,3)=linspace(59,59,3841);
grid59(:,4)= randperm(3841,3841);

dlmwrite('Random_grid_59x_coordinates.txt',grid59,'delimiter',' ','newline','pc')

