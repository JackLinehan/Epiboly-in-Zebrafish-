%% Animate all data from file: load it in and display on the same graph
% Load in data from file 
% Last edited 7.24.18
%% Edits 
% - creates new folder 
% - generates text files for each cell of stats to that folder 
% - 7.24.18 

%% Document

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
clear vars my_area my_intensity my_centroids my_perimeter my_outline k my_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Animate the Cells 
f = figure(1); % generate the figure space 
axis([-500 300 -300 300 -700 0]); %  200 450]);
view([154 19]);  %view angle 
ax = gca; 
%camroll(ax,-60); 
MP = get(0, 'MonitorPositions'); % How much monitor space is their? (i.e. one monitor or two?) 
%grid on 
set(gcf, 'Position',[MP(1,4)+1000 1 MP(1,4) MP(1,3)]); % When you have two screens
% sets up the figure location on screen 
set(gca,'color',[.5 .5 .5]); 

axis square
% N is the number of cells we will render
N = abs(length(area_all{1,1})); 
% x coordinate space matrix 
x = zeros(n,N); 
% y coordinate space matrix 
y = zeros(n,N); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spherical coordinates to generate sphere
% this section generates the "yolk sac" that our epidermal cells will be
% migrating over through epiboly. Including a sphere helps add a bit of
% depth perception. 

% Try to generate this sphere using as accurate a number as possible. 
r_1 =pi+1; % this is the number were taking to be the radius of the zfish embryo 

% basic steps to generating a sphere: set the limits 
phi_1 = linspace(-pi,pi,40);
theta_1 = linspace(0,pi,40);
[phi_1,theta_1]=meshgrid(phi_1,theta_1);

%[x_1,y_1,z_1] = sph2cart(theta_1,phi_1,.000250); 
x_1 = sin(phi_1).*cos(theta_1) ; 
y_1 = r_1*sin(phi_1).*sin(theta_1); 
z_1 = r_1*cos(phi_1); 
x_sphere = x_1.*70 - 25; 
y_sphere = y_1.*70 - 15; 
z_sphere = z_1.*70 - 385; 



%% rotate command 
direction = [0 1 0]; 
%% 
% j is the outer loop variable, and indictes the frame or point in time. 
for j = 1:N-1 
    hold on 
     s = surf(4.*x_sphere,y_sphere,z_sphere);
     s.EdgeColor = 'none'; % This removes the edge color from the 'yolk sac' sphere
     s.FaceColor = 'w';
    centroids = cell2mat(cellfun(@(x) x(j,:),centroids_all,'UniformOutput',false)); 
    x(:,j) = centroids(:,2)./1000000; 
    y(:,j) = centroids(:,1)./1000000; 
    
    outlines = cellfun(@(x) x{:,j},outline_all,'UniformOutput',false); 

    for m = 1:n
        
        hold on 
        
        x_1 = outlines{m,1}(:,2); 
        y_1 = outlines{m,1}(:,1);
        z = sqrt(pi+.01-(x_1.^2 + y_1.^2)); 
        %z = 1.*z.^2; 
        z=-1.*imag(z); 
        %z = z./1000000; 
        
       membranes = plot3(x_1,y_1,z,'r'); 
       set(gca,'Fontweight','bold'); 
       box on 
      %rotate(membranes,direction,15); 
        %plot3(xwrap,ywrap,zwrap); 
        
        xlabel('x');
        ylabel('y'); 
        zlabel('z'); 
        ax = f.Children; 
        sphere_axis = gca; 
        
        b = ax.Children; 
        hold off 
        
    end 
    hold on 
    y_2 =1000000.* x(:,1:j); % 1:j to include previous points
    x_2 = 1000000.*y(:,1:j); 
    
    z_2 = sqrt(pi^2-(x_2.^2 + y_2.^2)); 
    %z_2 = z_2; 
    z_2 = -1.*imag(z_2); 
    %[az_1, elev_1, r_1] = cart2sph(x_2,y_2,z_2);  
    
    plot3(x_2,y_2,z_2,'k.'); % centroids plots
    %rotate(cell_centroids,direction,-15); 
    hold off; 
    drawnow 
    disp(j); 
    pause(.1); 
    
    
    frame = getframe(f); % save the image in the figure 
    im = frame2im(frame); % convert that frame into a true image
    [imind, cm] = rgb2ind(im,256); % add some color to it 
    
    % Put the color image we just generated in the frame into a gif 
    if (j == 1)
        imwrite(imind, cm,'my_gif.gif','Loopcount',inf); 
    else 
        imwrite(imind, cm, 'my_gif.gif','WriteMode','append'); 
    end 
    
    if (j < N-1)
        delete(b); 
    end 
    
     
end 


%% Commented this section out 
% 
% %% Let's look at some measures, change in area n'position 
%  
% my_diff_area = cellfun(@(x) diff(x(1,1:end-1)),area_all,'UniformOutput',false); 
% my_diff_luminosity = cellfun(@(x) diff(x(1,1:end-1)),intensity_all,'UniformOutput',false);
% 
% my_diff_x = diff(x(:,1:end),1,2); 
% my_diff_y = diff(y(:,1:end),1,2); 
% 
% displacement = sqrt(my_diff_x.^2 + my_diff_y.^2); 
% 
% 
% %% Relative Luminosity 
% % Compare the luminosity of each cell to the others at any one time by
% % through normalization: luminosity / mean(all_luminosity)
% 
% relative_luminosity = zeros(n,N-1); 
% 
% for j = 1:N-1 
%     
%     my_lum = cell2mat(cellfun(@(x) x(:,j),intensity_all,'UniformOutput',false)); 
%     
%     norm_fac = mean(my_lum); 
%     
%     relative_luminosity(:,j) = my_lum./norm_fac; 
%     
% end 
% clear vars ax b j m x_1 y_1 x y 
% 
% 
% %% Write text sheet to export stats 
% filename = 'Trace Data'; 
% mkdir(filename); 
% cd(filename); 
% for j = 1:abs(length(area_all))
%     
%       filename = ['cell' num2str(j) '.txt']; 
%     
%       x = [relative_luminosity(j,:),0]; 
%       y = [displacement(j,:),0]; 
%       to_write = [area_all{j,1};intensity_all{j,1}; perimeter_all{j,1}; x; y]; 
%       
% 
%     
%         fid = fopen(filename, 'w'); 
%         
%         fprintf(fid,[' Area    ' '     Intensity ' '  perimeter ' '  Norm Intensity ' 'Displacement \n']); 
%         
%         fprintf(fid,'%f    %f    %f      %f     %f \n',to_write); 
% 
% 
% 
% 
% end 








    