dirInput = getDirectory("/media/vphys/2TB\ FireCuda/DarkSide-DEAP/microScopeStudy/10Nov2022/");
dirOutput = getDirectory("/media/vphys/2TB\ FireCuda/DarkSide-DEAP/microScopeStudy/10Nov2022/test11/");
list = getFileList(dirInput);
setBatchMode(true);
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
    	open(dirInput+list[i]);
    	run("8-bit");
        run("Gaussian Blur...", "sigma=2");
        run("Subtract Background...","rolling=100.0 light background");
	run("Set Measurements...", "area mean standard min centroid perimeter fit feret's");
        setAutoThreshold("Huang");
        run("Analyze Particles...", "  show=Overlay clear");



	saveAs("Results", "/media/vphys/My Passport/microScopePictures/15Nov2022/testResults/Results21.csv");
    	setAutoThreshold("Default");
    	//run("Threshold...");
    	getThreshold(lower, upper);
    	if (upper>90) {
    		setThreshold(0, 90);
    		run("Convert to Mask");
    		run("Close");
    		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00");// show=Nothing display clear summarize");
   
    	}
    	else {
    		setAutoThreshold("Default");
    		run("Convert to Mask");
    		run("Close");
    		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00");// show=Nothing display clear summarize");
   
    	}
    saveAs("TIFF", dirOutput+list[i]);
    close();

}
