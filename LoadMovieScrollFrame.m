%% Frame By Frame movie analysis 
% this script loads in a movie saved as a .tif and lets you to scroll
% through it frame by frame. 
% Scrolling forward is accomplished by clikcing the mouse 
% Scrolling backward is done by clicking any key 
% You will be warned if you approach the begging of the movie (by moving
% backwards) or the end of the movie 

fname = 'MAX_Ex1HDRF1timelapse5_12_18.tif'; % name of .tif file 
info = imfinfo(fname); % Returns graphics file information (i.e. format,filesize,colortype)
num_images = numel(info); % Returns the number of elements in info, i.e. the number of frames in the movie 
my_movie = cell(1,num_images); % Preallocate space to store images 



%% Playing around with image processing 
% This section was used to play around with some of the thresholding and
% filtering techniques built into matlab. 
% for k = 1:num_images
%     
%     [A,map] = imread(fname, k, 'Info', info); 
%     
%     %colormap map
%     %set(figure(1), 'Units', 'Normalized', 'OuterPosition',[0, 0.4,1,0.96]);  
%     KAverage = filter2(fspecial('average',3),A)/255; 
%     %histeq_image = histeq(A);
%     %histeq_image = medfilt2(histeq_image); 
%     %Kmedian = medfilt2(A); 
%     my_movie{1,k} = A; 
%     %fullfig
%     
%     imshow(KAverage, map); 
%     pause(0.1); 
% end 
%% Running through the movie 

for k = 1:num_images

    [A,map] = imread(fname, k, 'Info', info); %This reads in the movie frame by frame 
    
    %background = imopen(A,strel('disk',15)); 
       
    %A = A - background/2; 
  
    my_movie{1,k} = A; % save the frame as a matrix element of a preallocated cell array 
    
   
end 
count = 1; 
enhance = 1.5; %Increases brightness (intensity of image) 
imshow(my_movie{1,1}*enhance,map);
while (count>=0) 
 
direction = waitforbuttonpress; % How do you want to proceed 
    
   if(direction == 1 && count ~= 1)
        count = count -1; 
        imshow(my_movie{1,count}*enhance,map)
   elseif (direction == 0 && count ~= num_images) 
       count = count +1; 
       imshow(my_movie{1,count}*enhance,map)
   elseif (count == 1 || count == num_images) 
       if (count == 1)
            beep; 
            disp('warning you are at the begging of the movie and cannot move any further back, please click the mouse to move ahead');  
            direction = waitforbuttonpress; 
            count = count + direction; 
            imshow(my_movie{1,count}*enhance,map)  
       elseif (count == num_images)
           beep; 
           disp('warning, you are at the end of the movie'); 
           imshow(my_movie{1,count}*enhance,map); 
           direction = waitforbuttonpress; 
           count = count - direction; 
       end 
   end 
        
   
        
   
end     

