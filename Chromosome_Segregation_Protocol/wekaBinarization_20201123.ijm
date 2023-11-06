//wekaBinarization.ijm

//Measure chromatin masses during anaphase from live cells stained for DNA
//20201111
//Megan E. Kelley kelley.e.megan@gmail.com

//Setup and cleanup
run("Close All");
run("Collect Garbage");
roiManager("reset");
run("Bio-Formats Macro Extensions");
run("Brightness/Contrast...");

//Analyzes a directory of files
Dir = getDirectory("Choose an Input Directory ");
list = getFileList(Dir);
iterations=list.length;
AnalysisDir=Dir+"Analysis"+File.separator;
if (!File.exists(AnalysisDir)) File.makeDirectory(AnalysisDir);
MatDir=AnalysisDir+"MatLab"+File.separator;
if (!File.exists(MatDir)) File.makeDirectory(MatDir);
DNADir=AnalysisDir+"DNAresults"+File.separator;
if (!File.exists(DNADir)) File.makeDirectory(DNADir);
CropDir=AnalysisDir+"Crops"+File.separator;
if (!File.exists(CropDir)) File.makeDirectory(CropDir);

for (i=0; i<iterations; i++) {
if (endsWith(list[i], ".nd2")||endsWith(list[i],".tif")){
		run("Bio-Formats Importer", "open="+Dir+list[i]+" autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		//Get image information
		image_window=getImageID();
		Stack.getDimensions(width, height, channels, slices, frames);
		Filename = list[i];
		FileTitle = split(list[i], ".");
		aImage=getImageID();
		selectImage(aImage);
		aName=AnalysisDir+FileTitle[0];
		mName=MatDir+FileTitle[0];
		dName=DNADir+FileTitle[0];
		cName=CropDir+FileTitle[0];
		//Isolate DNA signal
		run("Duplicate...", "duplicate channels=2");
		run("Z Project...", "projection=[Max Intensity] all");
		maxDNA=getImageID();
		run("Grays");
		//Select region to analyze to avoid irrelevant signal
		run("Z Project...", "projection=[Max Intensity]");
		run("Enhance Contrast", "saturated=0.35");
		setTool("freehand");
		waitForUser("Select boundary");
		selectImage(maxDNA);
		run("Restore Selection");
		setMinAndMax(100,250);
		waitForUser("Adjust boundary");
		//Select relevant timepoints
		waitForUser("Scroll to 1 frame prior to anaphase-onset");
		Stack.getPosition(channel, slice, frame);
		MinRange=frame;
		Stack.getDimensions(width, height, channels, slices, frames);
		MaxRange=frames;
		selectImage(aImage);
		run("Z Project...", "projection=[Max Intensity] all");
		run("Restore Selection");
		roiManager("reset");
		roiManager("Add");
		roiManager("Select", 0);
		roiManager("Save", aName+"_analysisROI.roi");
		run("Duplicate...", "duplicate frames="+MinRange+"-"+MaxRange+"");
		roiManager("Add");
		roiManager("Select", 1);
		roiManager("Save", aName+"_cropROI.roi");
		saveAs("Tiff", cName+"_crop");
		cropIm = getTitle();
		run("Duplicate...", "duplicate channels=1");
		saveAs("Tiff", aName+"_cropDIC");
		dicIm=getTitle();
		// //DIC segmentation
		// run("32-bit");
		// getStatistics(mean);
		// run("Subtract...", "value="+mean+" stack");
		// run("Abs", "stack");
		// run("Sobel Filter");
		// run("Subtract...", "value=1000 stack");
		// run("Mean...", "radius=25 stack");
		// run("8-bit");
		// run("Threshold...");
		// waitForUser("Apply Threshold");
		// // setAutoThreshold("Default dark stack");
		// // run("Convert to Mask", "method=Default background=Dark calculate black");
		// run("Fill Holes", "stack");
		// roiManager("Select", 1);
		// run("Make Inverse");
		// run("Set...", "value=0 stack");
		// run("Make Inverse");
		// saveAs("Tiff", mName+"_cropDICBinary");
		// cropDICIm=getTitle();

		selectImage(cropIm);
		roiManager("Select", 1);
		run("Duplicate...", "duplicate channels=2");
		saveAs("Tiff", aName+"_cropDNA");
		cropDNAIm=getTitle();
		//WeKa time
		selectImage(cropDNAIm);
		run("Trainable Weka Segmentation");
		setTool("wand");
		run("Wand Tool...", "tolerance=25 mode=Legacy");
		waitForUser("Classify DNA (class1) and background (class2)");
		call("trainableSegmentation.Weka_Segmentation.trainClassifier");
		waitForUser("Training ok?");
		call("trainableSegmentation.Weka_Segmentation.getResult");
		selectWindow("Classified image");
		run("Multiply...", "value=255 stack");
		run("Invert", "stack");
		run("Make Binary", "method=Default background=Dark calculate black");
		selectWindow("Classified image");
		run("Fill Holes", "stack");
		roiManager("Select", 1);
		run("Make Inverse");
		run("Set...", "value=0 stack");
		run("Make Inverse");
		waitForUser("Scroll to last usable frame");
		Stack.getPosition(NEWchannel, NEWslice, NEWframe);
		NewMaxRange=NEWslice;
		run("Duplicate...", "duplicate range=1-"+NewMaxRange+"");
		//run("Duplicate...","duplicate");
		saveAs("Tiff", mName+"_cropDNABinary");
		roiManager("Select", 1);
		run("Set Measurements...", "area centroid center perimeter bounding shape stack display redirect=None decimal=3");
		run("Analyze Particles...", "size=20.00-Infinity show=Outlines display exclude clear add stack");
		selectWindow("Results");
		saveAs("Text", dName+"_DNAresults.csv");
		close("Results");
		roiManager("Save", aName+"_DNAROIs.zip");
		// roiManager("reset");
		// selectImage(cropDICIm);
		// roiManager("Select", 1);
		// run("Set Measurements...", "area centroid center bounding stack display redirect=None decimal=3");
		// run("Analyze Particles...", "size=100.00-Infinity show=Outlines display exclude clear add stack");
		// selectWindow("Results");
		// saveAs("Text", aName+"_DICresults.csv");
		// close("Results");
		// roiManager("Save", aName+"_DICROIs.zip");
		roiManager("reset");
		run("Close All");
	}
}
