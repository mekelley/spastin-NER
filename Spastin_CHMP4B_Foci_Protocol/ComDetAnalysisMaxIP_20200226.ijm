run("Close All");
run("Bio-Formats Macro Extensions");
run("Brightness/Contrast...");
roiManager("reset");
Table.create("AllResults");


Dir = getDirectory("Choose an Input Directory ");
list = getFileList(Dir);
iterations=list.length;

addon = getString("Type identifier for csv files:", " ");

ComDetDir=Dir+"ComDet"+File.separator;
if (!File.exists(ComDetDir)) File.makeDirectory(ComDetDir);

for (i=0; i<iterations; i++) {
	if (endsWith(list[i], ".tif") || endsWith(list[i], ".nd2") || endsWith(list[i], ".ome")){
		run("Bio-Formats Importer", "open="+Dir+list[i]+" autoscale color_mode=Grayscale concatenate_series open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
		image_window=getImageID();
		run("Duplicate...", "duplicate channels=2");
		fluorChannel=getImageID();
		setMinAndMax(100, 250);
		waitForUser("Evaluate number of cells for analysis");
		cells = getNumber("Number of Cells", 1);
		if (cells > 0){
			run("Z Project...", "projection=[Sum Slices] all");
			run("Z Project...", "projection=[Sum Slices]");
			run("Enhance Contrast", "saturated=0.35");
			projection=getImageID();
			ComDet(fluorChannel,projection,cells);
		}
		else{
			run("Close All");
		}
	}
}
selectWindow("AllResults");
saveAs("results", ComDetDir+addon+"_AllResults.csv");
AllResultsWindow=getInfo("window.title");
close(AllResultsWindow);
run("Close All");

function ComDet(image, projection, cellROIs){
	roiManager("reset");
	setTool("freehand");
	Filename = list[i];
	FileTitle = split(list[i], ".");
	for(n=0; n<cellROIs; n++){
		selectImage(projection);
		waitForUser("Select boundary box");
		roiManager("Add");
		selectImage(image);
		run("Z Project...", "projection=[Max Intensity] all");
		run("Grays");
		run("Restore Selection");
		setMinAndMax(100, 250);
		waitForUser("Adjust boundary");
		m = n+1;
		CellName = FileTitle[0]+"_cell"+m;
		ROIname = ComDetDir+File.separator+FileTitle[0]+"_cell"+m;
		roiManager("Save", ROIname+"_analysis_cellROI.zip");
		waitForUser("Scroll to 1 frame prior to anaphase-onset");
		MinRange=getSliceNumber();
		waitForUser("Scroll to last relevant frame");
		MaxRange=getSliceNumber();
		run("Duplicate...", "duplicate range="+MinRange+"-"+MaxRange+"");
		uncorrName = getImageID();
		saveAs("Tiff", uncorrName+"_crop");
		//preCyto=getImageID();
		//CytoBackground(preCyto);
		run("Subtract Background...", "rolling=5 stack");
		//ComDet v.0.4.2

//SPAST
//		run("Detect Particles", "include segment ch1a=2 ch1s=5 add=[All detections]");
////CHMP4B
	run("Detect Particles", "include segment ch1a=2 ch1s=10 add=[All detections]");
		detections=getImageID();
		doCommand("Start Animation [\\]");
		midscreen = 0.4*screenWidth;
		Dialog.create("Detection Evaluation");
		items = newArray("Yes", "No");
		Dialog.addRadioButtonGroup("Is this detection accurate?", items, 1, 2, "One");
		Dialog.setLocation(midscreen,0);
		Dialog.show();
		EvaluateRotation = Dialog.getRadioButton();
		if (EvaluateRotation=="Yes"){
			selectWindow("Summary");
			comdetParticles=Table.getColumn("Number_of_Particles");
			selectWindow("AllResults");
			Table.setColumn(CellName, comdetParticles);
			Table.update;
			selectWindow("Results");
			saveAs("Text", ROIname+"_results.csv");
			ResultsWindow=getInfo("window.title");
			close(ResultsWindow);
			selectWindow("Summary");
			saveAs("Text", ROIname+"_ComDetSummary.csv");
			SummaryWindow=getInfo("window.title");
			close(SummaryWindow);
			roiManager("reset");
		}else if (EvaluateRotation=="No"){
			close("Results");
			close("Summary");
			close(detections);
			selectWindow("AllResults");
			error = newArray("Error");
			Table.setColumn(CellName, error);
			Table.update;
			roiManager("reset");
			}
		}
	run("Close All");
}

function CytoBackground(image, Filename){
	//middle = round(slices/2);
	//Stack.setSlice(middle);
	waitForUser("Evaluate number of cells for analysis");
	cells = getNumber("Number of Cells", 1);
	if (cells > 0){
		for(n=0; n<cells; n++){
			m = n+1;
			ROIname = Filename+"_cell"+m;
			waitForUser("Scroll to 1 frame prior to anaphase-onset");
			setTool("point");
			waitForUser("Select cytosolic background region");
			roiManager("Add");
			Roi.getCoordinates(xpoints, ypoints);
			xpoint=(xpoints[0]-20);
			ypoint=ypoints[0]-20;
			makeRectangle(xpoint, ypoint, 40, 40);
			midscreen = 0.4*screenWidth;
			Dialog.create("Background Evaluation");
			items = newArray("Yes", "No");
			Dialog.addRadioButtonGroup("Is this region correct?", items, 1, 2, "One");
			Dialog.setLocation(midscreen,0);
			Dialog.show();
			EvaluateRotation = Dialog.getRadioButton();
			if (EvaluateRotation=="Yes"){
				run("Measure");
				Avg = getResult('Mean');
				XtopL = getResult('BX');
				YtopL = getResult('BY');
				width = getResult('Width');
				height = getResult('Height');
				slice = getResult('Slice');
				myTable(ROIname,Avg,XtopL,YtopL,width,height,slice);
			}else if (EvaluateRotation=="No"){
				CytoBackground(image, Filename);
			}
		}
	}
	close();
	close();
}

function myTable(a,b,c,d,e,f,g){
	title1="CytoBackground";
	title2="["+title1+"]";
	if (isOpen(title1)){
   		print(title2, a+"\t"+b+"\t"+c+"\t"+d+"\t"+e+"\t"+f+"\t"+g);
		}
	else{
   		run("Table...", "name="+title2+" width=600 height=400");
   		print(title2, "\\Headings:Filename\tCytosolic_Background\tBX\tBY\tWidth\tHeight\tZ_Slice");
   		print(title2, a+"\t"+b+"\t"+c+"\t"+d+"\t"+e+"\t"+f+"\t"+g);
	}
}
