//Hyperstack.ijm

//Reslice hyperstacks that don't properly load (first export multipoints from NIS-Elements to TIFF)
//20201001
//Megan E. Kelley kelley.e.megan@gmail.com

//Setup and cleanup
run("Close All");
run("Collect Garbage");
roiManager("reset");
run("Bio-Formats Macro Extensions");
run("Brightness/Contrast...");

NewChannels=3;
NewSlices=2;

setBatchMode(true);

//Analyzes a directory of nd2 files
Dir = getDirectory("Choose an Input Directory ");
list = getFileList(Dir);
iterations=list.length;

for (i=0; i<iterations; i++) {
	if (endsWith(list[i], ".tif")){
		open(Dir+list[i]);
		run("Set Scale...", "distance=9.0909 known=1 unit=micron global");
		Stack.getDimensions(width, height, channels, slices, frames);
		TimePoints=slices/(NewChannels*NewSlices);
		run("Stack to Hyperstack...", "order=xyczt(default) channels="+NewChannels+" slices="+NewSlices+" frames="+TimePoints+" display=Grayscale");
		run("Save");
		close();
}
}

beep();
