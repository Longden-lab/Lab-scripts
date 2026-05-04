// Before running this script, 
//You should draw a line on the PA in the direction of most perpendicular to the PA and less affected by Z movement

run("8-bit");


if (!(selectionType()==0 || selectionType==5 || selectionType==6))
       exit("Line or Rectangle Selection Required");
    Roi.setStrokeWidth(7);
	run("Measure");
	min = getResult("Min", nResults - 1);
	max = getResult("Max", nResults - 1);
	setMinAndMax(min, max); // Set the current min and max
	run("Apply LUT"); 
	run("Enhance Contrast", "saturated=0.0");       
    setBatchMode(true);
	setOption("ScaleConversions", true);
	 
     run("Plot Profile");
     Plot.getValues(x, y);
     run("Clear Results");
     for (i=0; i<x.length; i++)
         setResult("x", i, x[i]);
     close();

     n = nSlices;
     for (slice=1; slice<=n; slice++) {
         showProgress(slice, n);
         setSlice(slice);
         profile = getProfile();
         sliceLabel = toString(slice);
         sliceData = split(getMetadata("Label"),"\n");
         if (sliceData.length>0) {
             line0 = sliceData[0];
             if (lengthOf(sliceLabel) > 0)
                 sliceLabel = sliceLabel+ " ("+ line0 + ")";
         }
         for (i=0; i<profile.length; i++)
             setResult(sliceLabel, i, profile[i]);
     }
     setBatchMode(false);
     updateResults;


dir = getDirectory("Choose a Directory ");
name=getString("Enter the name","default");


saveAs("Results", dir+name+".csv");
close("*");

