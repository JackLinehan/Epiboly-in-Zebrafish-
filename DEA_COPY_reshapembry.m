%% Deformed Embryo Copy - Making adjustments to embyro shape. 
rot_my_object = 55; %y-component
rot_x = -47; 
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
axis([-.4 2 -1 1.5 0 2.3]); %  200 450]);
view([-37 21]);  
MP = get(0, 'MonitorPositions'); 

set(gcf, 'Position',[MP(1,4)+1000 1 MP(1,4)+200 MP(1,3)+200]); % When you have two screens
set(gcf,'color','k'); 

N = abs(length(area_all{1,1})); 
%x = zeros(n,N); 
%y = zeros(n,N); 
%% Embryo shape as written in evl.m - original
% 
% r = .5; 
% phi=linspace(0,pi,40); % Input to meshgrid - [0,pi] is 180 degrees total 
% theta=linspace(-pi,pi,40);
% [phi,theta]=meshgrid(phi,theta); % this is where our grid is created 
% x_n = r*sin(phi(:,1:end/2)).*cos(theta(:,1:end/2)); 
% x_m = r*sin(phi(:,end/2:end)).*cos(theta(:,end/2:end)); 
% x = [x_n x_m]; 
% y_n = r*sin(phi(:,1:end/2)).*sin(theta(:,1:end/2)); 
% y_m = r*sin(phi(:,end/2:end)).*sin(theta(:,end/2:end)); 
% y = [y_n y_m]; 
% 
% gradient =[ 1.21, 1.21, 1.205, 1.2005,1.2,flip(1.0:.05:1.2)]; 
% %gradient_flip = [.1 ones(1,length(phi) - 2*length(gradient)) flip(gradient)]; %original copy
% gradient_flip = [.1 ones(1,length(phi) - 2*length(gradient)) ones(1,length(gradient))];
% z_part_one = gradient.*(r*cos(phi(:,1:10)));
% z_part_two = r*cos(phi(:,10:end)).*gradient_flip; 
% z = [z_part_one z_part_two]; 
% 
% 
% x_sphere = x + 1.2; 
% y_sphere = y; 
% %z_sphere = rot90(z + 2,2); 
% z_sphere = -1.*z; 
% clear vars x y z 

%% Reshape Origina YolkSac
r = .5; 
phi=linspace(0,pi,40); % Input to meshgrid - [0,pi] is 180 degrees total 
theta=linspace(-pi,pi,40);
[phi,theta]=meshgrid(phi,theta); % this is where our grid is created 
x_n = r*sin(phi(:,1:end/2)).*cos(theta(:,1:end/2)); 
x_m = r*sin(phi(:,end/2:end)).*cos(theta(:,end/2:end)); 
x = [x_n x_m]; 
y_n = r*sin(phi(:,1:end/2)).*sin(theta(:,1:end/2)); 
y_m = r*sin(phi(:,end/2:end)).*sin(theta(:,end/2:end)); 
y = [y_n y_m]; 

gradient =[ 1.21, 1.21, 1.205, 1.2005,1.2,flip(1.0:.05:1.2)]; 
%gradient_flip = [.1 ones(1,length(phi) - 2*length(gradient)) flip(gradient)]; %original copy
gradient_flip = [.1 ones(1,length(phi) - 2*length(gradient)) ones(1,length(gradient))];
z_part_one = gradient.*(r*cos(phi(:,1:10)));
z_part_two = r*cos(phi(:,10:end)).*gradient_flip; 
z = [z_part_one z_part_two]; 


x_sphere = x + 1.2; 
y_sphere = y; 
z_sphere = z; 
%z_sphere = rot90(z + 2,2); 
%z_sphere = -1.*z; % For embryo laying along the x axis with animal pole to
%the left. 
clear vars x y z 
%% Begin for loop for animation 

for j = 1:N-10
    hold on  
    s = surf(x_sphere - .7,y_sphere - .2,z_sphere + 1.1);
    rotate(s,[0 0 1],-38); 
    
    s.EdgeColor = 'none'; % This removes the edge color from the 'yolk sac' sphere
    s.FaceColor = 'w';
    centroids = cell2mat(cellfun(@(x) x(j,:),centroids_all,'UniformOutput',false)); 
    x(:,j) = centroids(:,2); % removed ./1000000
    y(:,j) = centroids(:,1); 
    
    outlines = cellfun(@(x) x{:,j},outline_all,'UniformOutput',false); 

    for m = 1:n
        
        hold on 
        
        x_1 = outlines{m,1}(:,2); 
        y_1 = outlines{m,1}(:,1);
        z = sqrt(2.23*1795-(x_1.^2 + y_1.^2)); % was 2.3 
        z = z.^2; 
        z = z./1000; 
        
        [az,elev,r] = cart2sph(x_1,y_1,z); 
        r = r./900; 
        
        membranes = plot3(az, elev,r,'Color',[0.8 0.17 .17]); 
        rotate(membranes, [0 1 0], rot_my_object);  
        rotate(membranes, [1 0 0], rot_x);

        set(gca,'Box','off','Projection','perspective','visible','off'); 
        ax = f.Children; 
        sphere_axis = gca; 
        b = ax.Children; 
        hold off 

    end 
    hold on 
    y_2 = x(:,1:j); % 1:j to include previous points
    x_2 = y(:,1:j); 
    
    z_2 = sqrt(2.23*1795-(x_2.^2 + y_2.^2)); % was 2.2
    z_2 = (z_2).^2; 
    z_2 = z_2./1000; % removed 1000 
    [az_1, elev_1, r_1] = cart2sph(x_2,y_2,z_2);
    r_1 = r_1./900; 
    
    %az_1 = rot90(wrapTo2Pi(az_1)); 
    %elev_1 = rot90(wrapToPi(elev_1)); 
    
    centroids = plot3(az_1,elev_1,r_1,'k.');
    rotate(centroids,[0 1 0], rot_my_object); 
    rotate(centroids, [1 0 0], rot_x);
    hold off; 
    drawnow 
    disp(j); 
    %pause(.1); 
    
    
    frame = getframe(f); % save the image in the figure 
    im = frame2im(frame); % convert that frame into a true image
    [imind, cm] = rgb2ind(im,256); % add some color to it 
    
    % Put the color image we just generated in the frame into a gif 
    if (j == 1)
        imwrite(imind, cm,'my_gif.gif','Loopcount',inf); 
    else 
        imwrite(imind, cm, 'my_gif.gif','WriteMode','append'); 
    end 
    
    if (j < N-10)
    delete(b); 
    end 
end 

