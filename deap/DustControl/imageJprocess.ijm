// Strange thing: the pictures must be put in a folder without any other files except images! otherwise it can not find pictures!

dirInput = getDirectory("~/MyGithub/microscopeDarkside/pic/");
dirOutput = getDirectory("~/MyGithub/microscopeDarkside/resultsImageJ/");
list = getFileList(dirInput);
list = Array.sort(list);
setBatchMode(true);
suffix = '.png'
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
    	open(dirInput+list[i]);
        print(list[i]);
    	run("8-bit");
        run("Gaussian Blur...", "sigma=2");
        run("Subtract Background...","rolling=80.0 light background");
	run("Set Measurements...", "area mean standard min centroid perimeter fit feret's");
    	setAutoThreshold("Default");
        // run("Analyze Particles...", "  show=Overlay clear");
    	//run("Threshold...");
    	getThreshold(lower, upper);

        //Exclude on Edge ROI
        // w = getWidth();
        // h = getHeight();
        // Define the ROI
       //roi = new Rectangle(4712, 0, 4790-4712, 45-0); //x (from left), y (from top), width, height
       // Set the ROI
       //setRoi(roi);

       // Do something with the ROI, for example, measure the mean gray value
       // mean = getStatistics().mean;
       // print("Mean gray value in ROI: " + mean);

        // run("Find Edges");
        //NOTE: the enhance contrast needs to be tested before using, usually turn-off. Too many noises. However it can boost the particle detection.
        // run("Enhance Contrast...", "saturated=0.02 normalize");
        // setAutoThreshold("MaxEntropy");
    	if (upper>90) 
        {
           setAutoThreshold("MaxEntropy");
           //setThreshold(90, 220); //hard-coded number
           //run("Convert to Mask");
           //run("Close");
           run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 exclude");// show=Nothing display clear summarize");

           //run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 display");// show=Nothing display clear summarize");

           saveAs("Results", dirOutput+File.nameWithoutExtension+"_Results.csv");
           //saveAs("Results", dirOutput+String.pad(list[i],4)+"_Results.csv");

           // run("Flatten");
           saveAs("JPG", dirOutput+list[i]);//.JPG has smaller file size
    	}
    	else {
    		setAutoThreshold("Default");
    		run("Convert to Mask");
    		run("Close");
    		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00");// show=Nothing display clear summarize");
    	}
    close();
}
run("Close All");


// This tool is a wand tool that also runs the Measure command

//    macro "Wand Auto-Measure Tool -C00b-Lee22-o2244" {
//        requires("1.30k");
//        getCursorLoc(x, y, z, flags);
//        doWand(x,y);
//        if (selectionType!=0)
//            run("Measure");
//    }



//  saveSettings;
//  tolerance = 0.0000001;
//  labels = newArray("Area","Mean","StdDev","Min","Max", "X","Y",
//    "XM","YM","Perim.","Major","Angle", "Circ.","Feret");
//  run("AuPbSn 40 (56K)");
//  run("Set Measurements...", "area mean standard min centroid perimeter fit feret's");
//  setAutoThreshold("Huang");
//  run("Analyze Particles...", "  show=Overlay clear");
//  //resetThreshold;
//  for (i=0; i<nResults; i++) {
//     Overlay.activateSelection(i);
//     for (j=0; j<labels.length; j++) {
//        label = labels[j];
//        List.setMeasurements("limit");
//        v1 = getResult(label,i);
//        v2 = List.getValue(label);
//        if (abs(v1-v2)>tolerance || (!isNaN(v1)&&isNaN(v2)) )
//           print(label+"["+i+"]: "+v1+"  "+v2);
//
//     }
//  }
//  restoreSettings;
