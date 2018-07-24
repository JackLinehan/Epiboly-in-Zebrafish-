%% BUILD: Trace boundary through movie 
% This program allows you to trace the outline of a cell, and provides that
% outline's centroid, area, perimeter. You proceed
% through one full run of the movie tracking one cell, get to the end and
% are given the option to do it all again for another cell. 
% Jack Linehan (6.8.18)
% Follow That Cell: Motility Analysis 
%% Edits log 
% - user decides where to start and end tracing 
% - if the figure window is closed data is still saved to generated file
% - 7.24.18 

%% Frame By Frame movie analysis 
% this script loads in a movie saved as a .tif 

fname = uigetfile(); 

[info] = imfinfo(fname); % Returns graphics file information (i.e. format,filesize,colortype)

%% Let the user decide where to start tracing and stop tracing 

first_frame = input('At what frame do you want to begin tracing cells?:'); 
num_images = input('At what frame do you want to finish tracing cells?:'); 
enhance = 1.5; % brighten the image up a bit
flag = 1; % counts the number of cells you track 

%% Create file for storing shape data
% We create a file whose name is the date and time of its creation use the
% following to avoid issues with colons, which cannnot be used in strings

filename = datestr(now,'mm-dd-yyyy HH-MM'); 
mkdir(filename); 

cell_number = 0; 
while (cell_number == 0) 
%% Preallocate the data were intersted in looking at
total_space = (num_images - first_frame)+1; 

my_movie = cell(1,total_space);      % Preallocate space to store images 
my_outline = {1,total_space};        % store outline 
my_centroids = zeros(total_space,2); % store centroid 
my_area = zeros(1,total_space);      % preallocate area 
my_perimeter = zeros(1,total_space); % preallocate perimeter
my_intensity = zeros(1,total_space); % Preallocates intensity 
%% Running through the movie 

count = 1;

for k = first_frame:num_images
   
    [A,map] = imread(fname, k, 'Info', info); %This reads in the movie frame by frame 
    
    my_movie{1,count} = A;      % save the frame as a matrix element of a preallocated cell array 
    
    count = count +1;
    
end 

count = 1;

handle = figure(1); 
 
while (count <= total_space) 
    
        if (count >1)
        direction = waitforbuttonpress; 
        % if your satisified with the shape you have drawn then left click
        % the mouse, if not press any key and you erase your work and get
        % to do it again. 
            if (direction ==1 )
                count = count -1; 
            end
        end 
        
fontSize = 16; 
grayImage = my_movie{1,count}; %imread(fullFileName);

imshow(my_movie{1,count}*enhance,map);
set(handle,'Position',get(0,'Screensize')); 
% set(handle,'pointer','crosshair'); changes the pointer style of the
% mouse: works on the first frame but does not on the followers 

disp(count+first_frame); 
axis on;

if (count == first_frame)
    message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish. \n If your upset with your work hit space bar and you will get to try again');
    uiwait(msgbox(message)); 
end 

hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
    if (ishghandle(handle) == 0) % check to make sure the figure still exists, if they delete it the program saves whatever data they collected for however far they made it through the movie
        count = total_space; % finish the while loop and save whatever data is in the workspace 
    else 
     binaryImage = hFH.createMask();
     %xy = hFH.getPosition;
    end 


% Label the binary image and compute the centroid and center of mass.

measurements = regionprops(binaryImage, grayImage, ...
    'area', 'Centroid', 'WeightedCentroid', 'Perimeter');
my_area(1,count) = measurements.Area; 
my_centroids(count,:) = measurements.Centroid; 
centerOfMass = measurements.WeightedCentroid; 
my_perimeter(1,count) = measurements.Perimeter; 
my_intensity(1,count) = mean(my_movie{1,count}(binaryImage)); 

% Calculate the area, in pixels, that they drew.
numberOfPixels1 = sum(binaryImage(:)); 
% Another way to calculate it that takes fractional pixels into account.
numberOfPixels2 = bwarea(binaryImage); 

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries = bwboundaries(binaryImage);
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
my_outline{1,count} = xy; 
x = xy(:, 2); % Columns.
y = xy(:, 1); % Rows.
y = -1.*y; 

count = count + 1;  

end 

%% Now we can save our workspace variables to the folder
 % We really only want to save the following: 
 %
 % my_outline
 % my_centroids
 % my_area
 % my_perimeter 
 % my_intensity 

 
 h = pwd(); % get current directory 
 cd(filename); % go to new folder 
 
 save(['shape_data', num2str(flag)],'my_outline','my_area','my_centroids','my_perimeter','my_intensity'); % save the following data as shape_data(n)
 
 flag = flag + 1; % Move to the next cell 
 
 cd(h); % Go back to the directory fullfilmboundarytrace is saved in 
 
message = sprintf('Left click to track a new cell'); % prompt to be displayed to the user 

uiwait(msgbox(message)); % prompt user that a new run has begun

h = waitforbuttonpress; % Wait for them to move forward 

if (ishghandle(handle) == 0) % Do they want to continue tracing cells, or do they want to quit 
    break 
end 

if (h == 0)
    cell_number = h; 
end 
    
end 




