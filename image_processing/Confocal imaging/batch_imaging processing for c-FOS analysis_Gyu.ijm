// This script is based on the conocal images recorded with two channels
// Number of channels, color of channel, threshold can be adjusted
// After running this script in FIJI, Qupath is needed for actual analysis

inputDir = getDirectory("Choose Source Directory"); // select Directory where cFOS images are
outputDir = inputDir + "Processed images" + File.separator; // Select Directory for the output file
if(!File.exists(outputDir)){
	File.makeDirectory(outputDir);
}

list=getFileList(inputDir);

setBatchMode(true);

for(i=0;i<list.length;i++){
	filename=list[i];
	if (endsWith(filename, ".nd2")){
		open(inputDir+filename);	
		path=getTitle();
		width = getWidth();
		height = getHeight();

		run("Despeckle", "stack");
		run("Deinterleave", "how=2 keep");

		selectImage(path+" #1"); //Select Channel for DAPI
		run("Measure");
		min = getResult("Min", nResults -1);
		max = getResult("Max", nResults - 1);
		setMinAndMax(min+10, max); // Set the current min and max

		selectImage(path+" #2"); // Select Channel for cFOS
		run("Measure");
		min = getResult("Min", nResults -1);
		max = getResult("Max", nResults - 1);
		setMinAndMax(700, 1500); // Set the current min and max


		run("Merge Channels...", " c3=["+path+" #1] c1=["+path+" #2] ignore"); // Merge DAPI and cFOS - DAPI as blue, cFOS as red
		saveAs("Tiff",outputDir+filename+".tif");
		run("Close All");
		run("Clear Results");
	}
}
setBatchMode(false);
showMessage("Batch Processing Complete!");
