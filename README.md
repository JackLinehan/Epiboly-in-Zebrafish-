# Epiboly-in-Zebrafish-
MatLab codes for analysis of collective cell migration during Zfish epiboly

FullFilmBoundary Trace (Protocol) 
The movies you use in this program need to be max projections. (Do this in Fiji). 
To use the ROI Manager, fullfilmboundarytrace.m, you need to open MatLab. Double click the MatLab icon on the desktop. 
Once the MatLab home screen is open click on the folder with an arrow on top of it. This is directly below the New Live Script tab on the toolbar. Navigate to the folder “EVL-Epiboly”. Double click it. Now open the folder ROI Manager. 
In the command window, type the following:
">>  fullfilmboundarytrace"
A new window will appear for you to select your movie. Above the cancel button in the bottom right hand corner there is a drop down menu. Select “ALL Files”, navigate to your movie, and double click it. If you do not select “ALL Files”, you will not be able to find your movie. Double Click on it. 
Your movie will now appear in a separate window with instructions as to how to trace the cells in the image. Each time you run fullfilmlboundarytrace a new folder is generated. In this folder, you will find shape data stored as .m files. Each file has the date and time fullfilmboundarytrace.m was run as its’ name. Feel free to rename these files; this is where your tracking data is stored. 
Proceed through the movie tracking the same cell. Once you decide to move forward, you cannot move back. The command window will display the current frame. 
Stats: fullfilmboundarytrace.m collects the following statistics for the shape you have outlined. 
Area – computed based of the shape you have drawn. 
Centroid – Center of the shape in terms of the image. 
Outline – Gives you the locations of each pixel in the image that has been selected using the trace tool you use to draw the outline. 
Perimeter – The number of pixels highlighted by the trace tool. 
Intensity – The average value of each pixel in the selected region. 

Each of these is saved as a .m file into the folder generated by the program. Each cell is saved as its own .m file (shape_datan.mat).
Once you’ve finished tracing cells you can proceed to the analysis section. This is pretty easy. When you first ran fullfilmboundarytrace.m you generated a new folder, and all you need to do is to copy paste the function animate.m into it. So, move that cursor over to the Current Folder tab on the left hand side of the matlab interface. Look and find the function animate.m. Copy it using ctrl + c, or if for some reason your on a mac, hold ctrl and left click the mouse. You should see a drop down menu, select copy. Now find the new folder you just generated (it’s the date-time that you first ran fullfilmboundarytrace.m) and double click it. Now click in the white space within the Current Folder tab so that the title bar is in blue. Now, paste animate over using ctrl + v, on either mac or pc. 

Animate\\
Animate displays the data collected using fullfilmboundarytrace.m as an animation of the cell outlines and centroids through time. It does all the data handling for you, so you do not have to ever see any numbers (unless you want to), just graphs and animations. 2D data is projected into 3D with spherical coordinate transformation. 
To run animate, enter the following into the command line: 
">> animate"
You’ll see a figure pop up and an animation run of the outlines you drew, and there locations through time. 
This function also computes various statistics for each cell. 
•	Computes the change in area with time 
•	Computes displacement 
•	Relative luminosity – compares luminosity of each cell to one another at each point in time
•	Change in luminosity 


DEA_COPY_reshapembryo.m\\

This code projects the 2D coordinates collected from the max projection into 3D using spherical coordinate transformations. To help visualize the process, a digital yolk sac is created. The code generates an animation, and a gif, of the animation, of the EVL cells as they proceed through epiboly and into gastrulation. 

