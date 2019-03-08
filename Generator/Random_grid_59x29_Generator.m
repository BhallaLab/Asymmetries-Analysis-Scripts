clearvars
clc

grid59 = zeros(841,4);
grid59(:,1)=linspace(1,841,841);
grid59(:,2)=linspace(59,59,841);
grid59(:,3)=linspace(59,59,841);

coord = zeros(841,1);
randomCoord = randperm(841,841);
randomCoord =randomCoord';
k = 1;

for i =1:29
    for j=1:29
        coord(k)=(15+(i-1))*59+15+j;
        k=k+1;
    end
end

grid59(:,4)= coord(randperm(841));

dlmwrite('Random_grid_59x29_coordinates.txt',grid59,'delimiter',' ','newline','pc')

clear i j k