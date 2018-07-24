%% Animate all data from file: load it in and display on the same graph
% Load in data from file 
% Last edited 7.24.18
%% Edits 
% - creates new folder 
% - generates text files for each cell of stats to that folder 
% - 7.24.18 

my_data = dir ('**/*.mat'); 

n = max(size(my_data)); 

area_all = cell(n,1); 
centroids_all = cell(n,1); 
intensity_all = cell(n,1); 
outline_all = cell(n,1); 
perimeter_all = cell(n,1); 

for k = 1:n 
    
    load(my_data(k).name); 
    
    area_all{k,1} = my_area; 
    centroids_all{k,1} = my_centroids; 
    intensity_all{k,1} = my_intensity; 
    perimeter_all{k,1} = my_perimeter; 
    outline_all{k,1} = my_outline; 
    
end 
clear vars my_area my_intensity my_centroids my_perimeter my_outline k my_data
%% Animate the Cells 
f = figure(1); 
axis([0 500 0 475]); 
set(gca,'color','k'); 
N = abs(length(area_all{1,1})); 
x = zeros(n,N); 
y = zeros(n,N); 
for j = 1:N-1
    
    centroids = cell2mat(cellfun(@(x) x(j,:),centroids_all,'UniformOutput',false)); 
    x(:,j) = centroids(:,2); 
    y(:,j) = centroids(:,1); 
    
    outlines = cellfun(@(x) x{:,j},outline_all,'UniformOutput',false); 
    
    for m = 1:n
        
        hold on 
        
        x_1 = outlines{m,1}(:,2); 
        y_1 = outlines{m,1}(:,1);
        
        plot(y_1,x_1,'g'); 
        ax = f.Children; 
        b = ax.Children; 
        hold off 
    end 
        
    
    hold on 
    plot(x(:,1:j),y(:,1:j),'o');
    hold off; 
    
    disp(j); 
    pause(0.13); 
    
    if (j < N-1)
    delete(b); 
    end 
end 

%% Let's look at some measures, change in area n'position 
 
my_diff_area = cellfun(@(x) diff(x(1,1:end-1)),area_all,'UniformOutput',false); 
my_diff_luminosity = cellfun(@(x) diff(x(1,1:end-1)),intensity_all,'UniformOutput',false);

my_diff_x = diff(x(:,1:end),1,2); 
my_diff_y = diff(y(:,1:end),1,2); 

displacement = sqrt(my_diff_x.^2 + my_diff_y.^2); 


%% Relative Luminosity 
% Compare the luminosity of each cell to the others at any one time by
% through normalization: luminosity / mean(all_luminosity)

relative_luminosity = zeros(n,N-1); 

for j = 1:N-1 
    
    my_lum = cell2mat(cellfun(@(x) x(:,j),intensity_all,'UniformOutput',false)); 
    
    norm_fac = mean(my_lum); 
    
    relative_luminosity(:,j) = my_lum./norm_fac; 
    
end 
clear vars ax b j m x_1 y_1 x y 


%% Write text sheet to export stats 
filename = 'Trace Data'; 
mkdir(filename); 
cd(filename); 
for j = 1:abs(length(area_all))
    
      filename = ['cell' num2str(j) '.txt']; 
    
      x = [relative_luminosity(j,:),0]; 
      y = [displacement(j,:),0]; 
      to_write = [area_all{j,1};intensity_all{j,1}; perimeter_all{j,1}; x; y]; 
      

    
        fid = fopen(filename, 'w'); 
        
        fprintf(fid,[' Area    ' '     Intensity ' '  perimeter ' '  Norm Intensity ' 'Displacement \n']); 
        
        fprintf(fid,'%f    %f    %f      %f     %f \n',to_write); 




end 








    