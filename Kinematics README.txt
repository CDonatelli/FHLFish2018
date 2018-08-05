Kinematis Code README

Here are the rough instructions for how to use the kinematics
code in this repository. They will be updates as we go and
figure out what functions work better.

0) Setting up Matlab
 - Open Matlab
 - Set the path containining the GitRepo
	- Under the "Home" tab select "Set Path"
	- Select "Add with Subfolders" and navigate to the path
	  containing the functions
	- It will ask if you want to save it for later but you
	  can't on these computers since we don't have admin
  	  access -_-

1) Creating struct files
 - for this you will need the videos and the body lengths (in 
   mm)of all of the fish.
 - Make sure your working directory is the folder with all of
   yourvideos
 - Run Create2DFishStruct2018
 - SpeciesNumber = Create2DFishStruct2018('SpeciesNumber')
 - EX Aflav1 = Create2DFishStruct2018('Aflav1')
 - It will ask how many videos you have for this individual
 - Then you will enter a string with the names of the files
   for each trial. Be sure to use the 'apostrophe' to denote 
   a string
 - EX 'gp0321.avi'

2) Find Midlines
 - Your workspace should be clear for this to work
 - Again, make sure you are in the working directory containing
   your video files