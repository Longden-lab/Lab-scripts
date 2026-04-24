//Image that you want to analyze should be already open, and positioned at the forward
path=getTitle();
width = getWidth();
height = getHeight();

//This is for images recorded in 5th floor 2P
if (endsWith(path,".tif")){
run("Deinterleave", "how=2 keep");
selectImage(path+" #1"); 

// This will close the channel used for the picosprtizing. 
//Based on tyour image set up, you might need to close the other channel
close(); 

selectImage(path);
close(); // Close original file
run("Smooth", "stack");
run("Gaussian Blur...", "sigma=2 stack");
selectImage(path+" #2");
run("Grouped Z Project...", "projection=[Average Intensity] group=10");
selectImage(path+" #2");
}

//This is for images recorded in 6th floor 2P
else{
	run("Split Channels");
	selectImage("C1-"+path);
	close();
	run("Smooth", "stack");
	run("Gaussian Blur...", "sigma=2 stack");
}


run("Z Project...", "start=49 stop=50 projection=[Max Intensity]");

//After this, you may or may not run TurboReg after this script, to reduce XY movement
//TurboReg doesn't fit in the Macro, so you should run it manually.
