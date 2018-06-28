%% This program finds the pixel values within the regions drawn on the image
% This is a test program to work out the method that will be used in
% boundary_intensity_measure.m

[A,map] = imread('MAX_Ex1HDRF1timelapse5_12_18.tif', 1); % load in a frame 

imshow(A,map); % show that frame

Xif= imfreehand();  % Call free hand 

binaryImage = Xif.createMask(); % Create a mask as a binary image 

structBoundaries = bwboundaries(binaryImage); % Get the boundaries of that mask 

xy = structBoundaries{1}; % Get the xy coordinates of the boundary structure 

% At this point we load in just one frame, add a boundary, extract that
% boundary, and begin to investigate how we access the actual pixels in the
% ROI. 

% select all elements of A that exist within the space enclosed by the
% points in xy. xy demarcates the subset within A. 

% Treat xy as if it is a map. We need a key to make sense of that map, and 
% we take y to be our key, which will essentially act as our index to find
% correpsonding x values. 

Y_bar = xy(:,1); % this indicates the set of all Y values in xy 

my_unique_y = unique(sort(Y_bar)); 

% the max and min values do not contain subsets of A that correspond to
% pixel values, so luckily we get to skip them 

my_space = {1,abs(length(my_unique_y))-2}; 

for j = 2:abs(length(my_unique_y))-1
    
    index = find(xy(:,1) == my_unique_y(j)); % this tells me what row in xy my corresponding Y_bar value exists in
    
    my_pits = find(abs(diff(xy(index,2))) >1); % this tells me if there is space between the points along x 
    
    x_set = xy(index,2); % this is the set of x values for some Y_bar
    
    my_subspace = {1,abs(length(my_pits))}; 
    
    for k = 2:abs(length(my_pits))
        
        my_subspace{k-1} = A(x_set(my_pits(k-1)):x_set(my_pits(k)),my_unique_y(j)); % Now we grab the subset of A that exists at Y betwen x_1 and x_2 
        
    end 
    
    my_space{j} = my_subspace; 

end 



% Now we have to actually check to see if we have, what are essentially
% chunks, in our boundary, something that would look like this
% |_____|------|________| we do this by looking for "pits" in our set. i.e
% continuos sequences of numbers in index. This is actually straight
% forward, all i need to do is check to make sure that index does not
% contain more than 2 elements, otherwise it is the max value or the min
% value. 

% We do not want to include the boundary in our intensity measurement, So
% we need find that there is a difference of at least one between x values
% in our index. I.E. if for some Y, x_2 - x_1 != 0, then we have the
% locations of pixels at Y. (this would be (Y, x_1+1:x_2-1); 












