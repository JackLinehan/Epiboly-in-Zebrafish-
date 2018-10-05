% This code animates a sphere forming from the north pole to its
% equator over another sphere (or all the way to the 'vegetal pole') 


%% Setting up our "space" - this is a meshgrid, matrix/grid with y rows and x columns 

% Setting up the figure window 
figure
%ax1 = axes('Position',[0.1 0.1 0.7 0.7]); 
%ax2 = axes('Position', [0.65 0.65 0.28 0.28]); 

axis([-1,1,-1,1,-1.3,1.3])% set axis for all figures 
r = 1; 
phi=linspace(0,pi,40); % Input to meshgrid - [0,pi] is 180 degrees total 
theta=linspace(-pi,pi,40); % input to meshgrid - pi to -pi is 90 to 270 degrees 
% These coordinates give us a sphere 


[phi,theta]=meshgrid(phi,theta); % this is where our grid is created 

% These are to mess around with the grid pattern, hopefully not destroying
% the shape too much (you have to deform the shape to change the grid) 

%sinh_operator = (exp(phi) - exp(-phi)) / 2; 
%tanh_operator = (exp(theta) - exp(theta))/ (exp(theta) + exp(-theta)); 

% This section allows you to play around with any operators you may have
%written 
%phi = sinh_operator; 
%theta = tanh_operator; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructiong our x,y,z coordinates 

x_n = r*sin(phi(:,1:end/2)).*cos(theta(:,1:end/2)); 
x_m = r*sin(phi(:,end/2:end)).*cos(theta(:,end/2:end)); 
x = [x_n x_m]; 
y_n = r*sin(phi(:,1:end/2)).*sin(theta(:,1:end/2)); 
y_m = r*sin(phi(:,end/2:end)).*sin(theta(:,end/2:end)); 
y = [y_n y_m]; 



gradient =[ 1.21, 1.21, 1.205, 1.2005,1.2,flip(1.0:.05:1.2)]; 
gradient_flip = [ones(1,length(phi) - 2*length(gradient))    flip(gradient)]; 

z_part_one = gradient.*(r*cos(phi(:,1:10)));
z_part_two = r*cos(phi(:,11:end)).*gradient_flip; 
z = [z_part_one z_part_two]; 

% phi_1, theta_1 
r_1 = .99; 
phi_1 = linspace(-pi,pi,40);
theta_1 = linspace(-pi,pi,40);
[phi_1,theta_1]=meshgrid(phi_1,theta_1);


x_1 = r_1*sin(phi_1).*cos(theta_1); 
y_1 = r_1*sin(phi_1).*sin(theta_1); 
z_1 = r_1*cos(phi_1); 
% phi=linspace(0,pi,30);
% theta=linspace(0,2*pi,10);
% [phi,theta]=meshgrid(phi,theta);
% 
% 
% r_1 = 1; 
% x_1 = r_1*sin(phi).*cos(theta); 
% y_1 = r_1*sin(phi).*sin(theta); 
% z_1 = r_1*cos(phi); 
% Using Get Frame to record "movie" 

ax = gca; 
ax.NextPlot = 'replaceChildren';

F(40) = struct('cdata',[],'colormap',[]); 

%hold on % Uncommenting this reveals another sphere that is eveloped by the
%'EVL'
hold on
s = surf(x_1,y_1,z_1);
s.EdgeColor = 'none'; % This removes the edge color from the 'yolk sac' sphere
s.FaceColor = 'w';

count = 1; 
 
for j = 2:40
     
    x_in = x(:,1:j); 
    y_in = y(:,1:j); 
    z_in = z(:,1:j); 
    
    h = surf(x_in,y_in,z_in,'FaceAlpha',0.5); 
    h.FaceAlpha = 1; 
    h.FaceColor = 'k'; 
    h.EdgeColor = 'g'; 
    h.LineStyle = ':'; 
    set(gca,'Color','k'); 
    drawnow
    F(j) = getframe(gcf,[0 0 560 420]); 
   
    count = count + 1; 
    
    
     
    %x_in1 = x(:,j); 
    %y_in1 = y(:,j); 
    %z_in1 = z(:,j); 
    %input = [x_in1, y_in1]; 
    %input1 = [x_in,y_in,z_in]; 
    %plot3(ax1,x_in,y_in,z_in); 
    pause(0.01);  % Adjust timing 
    
 
end 
clear vars h j gradient theta theta_1 sinh_operator z_part_one z_part_two
clear vars ax count F gradient input input1 phi phi_1 r r_1 s x_in1 y_in1 z_in1
clear vars gradient_flip x_n x_m y_n y_m 
% fig = figure; 
% movie(fig,F(:,3:end),1)


%% Construct 4d data set by time frame 
% to note: 
% x,y,z - coordinates for 'EVL'
% x_1, y_1, z_1 - coordinates for sphere 

% How long is our 4D stack? It is length of z-1, because we want to start
% with some evl, not just have it magically appear. 

x_cum = x_1 + x; 
y_cum = y_1 +y; 





