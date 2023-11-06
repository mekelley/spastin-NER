run("Close All");
run("Collect Garbage");
roiManager("reset");
run("Clear Results");
if (isOpen("Log")) {
		 selectWindow("Log");
		 run("Close" );
	 }
run("Bio-Formats Macro Extensions");
run("Brightness/Contrast...");

setBatchMode(true);

Dir = getDirectory("Choose an Input Directory ");
list = getFileList(Dir);
iterations=list.length;

SplitDir=Dir+"Split"+File.separator;
if (!File.exists(SplitDir)) File.makeDirectory(SplitDir);

for (i=0; i<iterations; i++) {
	if (endsWith(list[i], ".tif")){
		run("Bio-Formats Importer", "open="+Dir+list[i]+" autoscale color_mode=Composite open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
		Stack.setChannel(1);
		setMinAndMax(0, 1500);
		Stack.setChannel(2);
		setMinAndMax(25, 75);
		run("Split Channels");
		title = getTitle();
		saveAs("Tiff", SplitDir+title);
    close(title);
		title2 = getTitle();
		saveAs("Tiff", SplitDir+title2);
    close(title2);
    }
	}



beep();
