///Macro by Abby Kroken, Loyola University Chicago, 2023


////PLEASE SPECIFY CHANNELS HERE BEFORE STARTING:

DAPIchannel = "1"  //channel to detect nuclei
p65channel = "2" //channel to measure, e.g., p65 for NFkB activation
	
//cleanup	
	close("\\Others");
	run("Set Measurements...", "mean redirect=None decimal=3");	
	roiManager("reset")
	run("Clear Results");

run("Bio-Formats Macro Extensions");

//User selects source directory and output directory.
source_dir = getDirectory("Select source directory")
output_dir = getDirectory("Select output directory")
file_list = getFileList(source_dir);

setBatchMode(true);

for(a=0; a<file_list.length; a++){
	if(endsWith(file_list[a],  ".nd2")){
		run("Bio-Formats Importer", "open=[" + source_dir + file_list[a] + "] autoscale color_mode=Composite concatenate_series open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

			//make a nuclei mask from DAPI channel
			raw_image = getImageID();
			selectImage(raw_image);
			run("Select None");
			run("Duplicate...", "title=dapi duplicate channels=DAPIchannel"); ///user defined variable above
			run("Mean...", "radius=3 stack");
			run("Grays");
			//setAutoThreshold("Minimum dark");////needs manual threshold
			
				setAutoThreshold("Li dark stack");
				//waitForUser("Please adjust threshold");   //optional manual threshold adjustment here, remove comment to use.
				
				
			run("Convert to Mask");
			run("Watershed"); /////comment this out if nuclei are bean-shaped or over-separated
			run("Analyze Particles...", "size=1000-Infinity pixel show=Masks exclude"); //// size filter, may be adjusted if smaller nuclei are present
			rename("dapi size filter");
			run("Analyze Particles...", "size=40-Infinity pixel show=[Count Masks] exclude add"); //get ROIs here
			rename("Count Masks of dapi");
			run("Clear Results");
		
			//select all each ROI in manager and measure mean intensity
			count = roiManager("count");
			array = newArray(count);
		 	 for (i=0; i<array.length; i++) {
		      	array[i] = i;
		      	selectImage(raw_image);
		      	Stack.setChannel(p65channel); //////ensure this is the p65 channel
		 		roiManager("Select", i);
				getStatistics(area, mean, min, max, std, histogram);
				setResult("Nucleus Mean", nResults, mean);
				updateResults();
		 		 }
		
			//clear ROI manager
			roiManager("reset")
		
		
		//grow nuclei and re-separate using marker controlled watershed
			selectWindow("Count Masks of dapi");
			run("Duplicate...", "title=binary_mask");
			setThreshold(1, 255);
			run("Convert to Mask");
			
			selectWindow("Count Masks of dapi");
			run("Duplicate...", "title=distmap");
			setThreshold(1, 255);
			run("Convert to Mask");
			run("Invert");
			run("Distance Map");
			run("Invert");
			setThreshold(235, 255);/////original was 225,255 for a 40x image. increase the first number for smaller cytoplasm 
			run("Convert to Mask");
			
			run("Marker-controlled Watershed", "input=distmap marker=[Count Masks of dapi] mask=distmap compactness=0 binary calculate use");
		
			run("8-bit");
			setThreshold(1, 255, "raw");
			//setThreshold(1, 255);
			run("Convert to Mask");
		
			run("Analyze Particles...", "size=5-Infinity pixel show=[Count Masks] add");
		
			
		//select all lines in ROI manager and measure means of grown mask
			selectImage(raw_image);
			Stack.setDisplayMode("grayscale");
			Stack.setChannel(p65channel); ////ensure this is the p65 channel
			
			count = roiManager("count");
			array = newArray(count);
		 	 for (i=0; i<array.length; i++) {
		     	array[i] = i;
				roiManager("Select", i);
				selectImage(raw_image);
		      	Stack.setChannel(2); //////ensure this is the p65 channel
				getStatistics(area, mean, min, max, std, histogram);
				setResult("Cell Mean", nResults-count+i, mean);
				updateResults();
				
				//do math for ratio
				nuc = getResult("Nucleus Mean", nResults-count+i);
				cell = getResult("Cell Mean", nResults-count+i);
				setResult("Nucleus / cell ratio", nResults-count+i, nuc/cell);
		 	 }
		 	
		 	//save results with file name 
		 	selectWindow("Results");
			saveAs("Results", output_dir + file_list[a] + ".csv");
			
			//merge nucleus map and cytoplasm maps for checking integrity
			run("Merge Channels...", "c1=[Count Masks of dapi] c2=[Count Masks of distmap-watershed] create");
			saveAs("Tiff", output_dir + file_list[a] + ".tif");
			
		 	 //cleanup
		 	 	roiManager("reset")
				run("Clear Results");
				close("*");
				
	}
}

