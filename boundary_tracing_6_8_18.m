%% BUILD: Trace boundary through movie 
% Jack Linehan (6.8.18)

%% Frame By Frame movie analysis 
% this script loads in a movie saved as a .tif and lets you to scroll
% through it frame by frame. 

% try to make so you can rewrite trace an object within the same frame, 
% save the data when your done to the directory so you don't overwrite it! 
% start saving the centorids

% Load in image data
fname = 'C:\Users\JLinehan\Desktop\EVL- Epiboly\MAX_Ex1HDRF1timelapse5_12_18.tif'; % name of .tif file 
%fname = 'C:\Users\JLinehan\Desktop\EVL- Epiboly\boundry tracing versions\MAX_003.lif - WT_6 timelapse.tif'; 
[info] = imfinfo(fname); % Returns graphics file information (i.e. format,filesize,colortype)
num_images = numel(info); % Returns the number of elements in info, i.e. the number of frames in the movie 
enhance = 2; % brighten the image up a bit
flag = 1; 

%% Create file for storing shape data
% We create a file whose name is the date and time of its creation use the
% following to avoid issues with colons, which cannnot be used in strings

filename = ['ShapeInfo_',datestr(now,'mm-dd-yyyy HH-MM')]; 
mkdir(filename); 
% I create a function to wrap around the bulk of the code. I do this
% because when your done tracing a cell through the entire movie, you
% probably intend to trace another. This way, we save all the tracing data
% to a new folder, and the data has its label (its packaged together). To
% go back and track another cell, I use the return command which brings you
% back to the start of the function that called it. 

% num_images = 4;  
cell_number = 0; 
while (cell_number == 0) 

%% Preallocate the data were intersted in looking at
my_movie = cell(1,num_images);      % Preallocate space to store images 
my_outline = {1,num_images};        % store outline 
my_centroids = zeros(num_images,2); % store centroid 
my_area = zeros(1,num_images);      % preallocate area 
my_perimeter = zeros(1,num_images); % preallocate perimeter
%% Running through the movie 
for k = 1:num_images
   
    [A,map] = imread(fname, k, 'Info', info); %This reads in the movie frame by frame 
    
    my_movie{1,k} = A;      % save the frame as a matrix element of a preallocated cell array 
end 
count = 1; 
figure(1); 
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

while (count <= num_images) 
    
        if (count >1)
        direction = waitforbuttonpress; 
        % if your satisified with the shape you have drawn than left click
        % the mouse, if not press any key and you erase your work and get
        % to do it again. 
            if (direction ==1 )
                count = count -1; 
                disp(direction)     
            end
        end 
fontSize = 16; 
grayImage = my_movie{1,count}; %imread(fullFileName);
imshow(my_movie{1,count}*enhance,map);
%set(figure(gcf), 'Position', get(0,'Screensize')); % Maximize figure.
title=(['Frame:',num2str(count)]); 
axis on;

%set(my_figure, 'Position', get(0,'Screensize')); % Maximize figure.
if (count == 1)
    message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish. \n If your upset with your work hit space bar and you will get to try again');
    uiwait(msgbox(message));
    %count = count + 1; 
end 
hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
binaryImage = hFH.createMask();
xy = hFH.getPosition;

% Label the binary image and computer the centroid and center of mass.
%labeledImage = bwlabel(binaryImage);
measurements = regionprops(binaryImage, grayImage, ...
    'area', 'Centroid', 'WeightedCentroid', 'Perimeter');
my_area(1,count) = measurements.Area; 
my_centroids(count,:) = measurements.Centroid; 
centerOfMass = measurements.WeightedCentroid; 
my_perimeter(1,count) = measurements.Perimeter; 

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
count = count + 1; 

end 

%% Now we can save our workspace variables to the folder
 % We really only want to save the following: 
 %
 % my_outline
 % my_centroids
 % my_area
 % my_perimeter 

 
 h = pwd(); 
 cd(filename); 
 
 save(['shape_data', num2str(flag)],'my_outline','my_area','my_centroids','my_perimeter'); 
 
 flag = flag + 1; 
 
 cd(h); 
 
message = sprintf('Left click to track a new cell'); 

uiwait(msgbox(message)); 

h = waitforbuttonpress; 

if (h == 0)
    cell_number = h; 
end 
end 




