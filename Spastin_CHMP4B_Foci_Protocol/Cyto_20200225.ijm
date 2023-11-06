//Cytoplasmic Background Intensity
//20200225
//Megan E. Kelley kelley.e.megan@gmail.com

run("Close All");
run("Bio-Formats Macro Extensions");
run("Brightness/Contrast...");
run("Set Measurements...", "area mean min bounding stack redirect=None decimal=3");
roiManager("reset");

Dir = getDirectory("Choose an Input Directory ");
list = getFileList(Dir);
iterations=list.length;

for (i=0; i<iterations; i++) {
	if (endsWith(list[i], ".tif") || endsWith(list[i], ".nd2") || endsWith(list[i], ".ome")){
		run("Bio-Formats Importer", "open="+Dir+list[i]+" autoscale color_mode=Grayscale concatenate_series open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
		image_window=getImageID();
		Filename = list[i];
		Maximize_TopLeft(image_window);
		selectImage(image_window);
		run("Duplicate...", "duplicate channels=2");
		run("Z Project...", "projection=[Max Intensity] all");
		run("Grays");
		setMinAndMax(100, 250);
		fluorChannel=getImageID();
		CytoBackground(image_window, Filename);
		close();
	}
}

selectWindow("CytoBackground");
saveAs("results", Dir+File.separator+"CytoBackground.csv");
run("Close");
selectWindow("Results");
run("Close");

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

function Maximize_TopLeft(input_window){
input_window = getTitle();
H=0.95*screenHeight;
W=0.7*screenWidth;
selectImage(input_window);
setLocation(0, 0, W, H);
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
