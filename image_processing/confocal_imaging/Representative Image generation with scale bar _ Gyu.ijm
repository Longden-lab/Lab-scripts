// Original Image should be opened and positioned at the forward. 
path=getTitle();
width = getWidth();
height = getHeight();

//Noise reduction and Channel seperation
if (endsWith(path,".nd2")){
run("Despeckle", "stack");
run("Deinterleave", "how=2 keep");


selectImage(path+" #1");
run("Measure");
min = getResult("Min", nResults -1);
max = getResult("Max", nResults - 1);
setMinAndMax(min+5, max); // Set the current min and max

selectImage(path+" #2");
run("Measure");
min = getResult("Min", nResults -1);
max = getResult("Max", nResults - 1);
setMinAndMax(min, max+5); // Set the current min and max


//automatic scale bar setting

scalebarsize = 0.1; // approximate size of the scale bar relative to image width

getPixelSize(unit,w,h);
if (unit == "pixels") exit("Image not spatially calibrated");

imagewidth = w*getWidth();  // image width in measurement units
scalebarlen = 1; // initial scale bar length in measurement units

// recursively calculate a 1-2-5 series until the length reaches scalebarsize, default to 1/10th of image width
// 1-2-5 series is calculated by repeated multiplication with 2.3, rounded to one significant digit
while (scalebarlen < imagewidth * scalebarsize) {
	scalebarlen = round((scalebarlen*2.3)/(Math.pow(10,(floor(Math.log10(abs(scalebarlen*2.3)))))))*(Math.pow(10,(floor(Math.log10(abs(scalebarlen*2.3))))));
}

scalebarsettings = "height=10 font=50 color=White background=None location=[Lower Right] bold label"; 

//run("Enhance Contrast", "saturated=0.0"); 
run("Merge Channels...", " c3=["+path+" #1] c1=["+path+" #2] ignore");
run("Scale Bar...", "width=&scalebarlen "+scalebarsettings);

}

saveAs("Tiff", "C:/Users/GB Park/OneDrive - University of Maryland School of Medicine/Desktop/Lab/data/confocal/251220 cFOS rb toy 1hr naive/Processed images/"+path+"with scale bar.tif");
