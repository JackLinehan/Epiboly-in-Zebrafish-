%% EVL Analysis 

%% Load Data 

% find the cell outline data saved as .mat file. 
my_data = dir ('**/*.mat'); 

% n is the number of cells that will be rendered in the animation - number
% of cells that were tracked. 
n = max(size(my_data)); 


% this is the data that was collected by the ROI project
area_all = cell(n,1); 
centroids_all = cell(n,1); 
intensity_all = cell(n,1); 
outline_all = cell(n,1); 
perimeter_all = cell(n,1); 


% go through and collect the data for each cell into a cell array so that
% it can be used in the animation. Generate a cell array containing the
% area, center, perimeter and outline for each cell. It is easier to work
% with this data in the cell format. 
for k = 1:n 
    
    load(my_data(k).name); 
    
    area_all{k,1} = my_area; 
    centroids_all{k,1} = my_centroids; 
    intensity_all{k,1} = my_intensity; 
    perimeter_all{k,1} = my_perimeter; 
    outline_all{k,1} = my_outline; 
    
end 

%% Cell Position and Velocity 
VX = cell(1,n); 
VY = cell(1,n); 
V1 = cell(1,n); 
VMEAN = zeros(1,n); 
for k = 1:n 
x = centroids_all{k,1}(:,1); 
y = centroids_all{k,1}(:,2); 

vx = zeros(1,39); 
vy = zeros(1,39); 
V = zeros(1,39); 

for j = 2:length(x)
    
    vx(1,j-1) = x(j) - x(j-1); 
    
    vy(1,j-1) = y(j) - y(j-1); 
    
    V(1,j-1) = sqrt(vx(1,j-1)^2 + vy(1,j-1)^2); 
    
end 

    VX{1,k} = vx; 
    VY{1,k} = vy; 
    V1{1,k}= V ; 
    VMEAN(1,k) = mean(V); 

end 

% cell centroid (x,y) coordinates 
plot(x,y,'*'); 
title('Centroids'); 
pause(); 
% change in location 
plot(diff(x)); 
title('speed in x'); 
pause(); 
plot(diff(y),'*-'); 
title('speed in y'); 
pause();
plot(V); 
title('average speed'); 

%% Cell Surface Area 

A = area_all{1,1}; 

plot3(x,y,A); 
grid on 
xlabel('x'); ylabel('y'); zlabel('Surface Area'); 

plot(x,A); 
grid on 
xlabel('x-position'); ylabel('Surface Area'); 
pause(); 
plot(y,A); 
grid on 
xlabel('y'); 
ylabel('Surface Area'); 









