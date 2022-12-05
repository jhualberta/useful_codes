dirInput = getDirectory("C:\Users\POOJA\lab work\experiments\Shankhamala\imagej batch\ ");
dirOutput = getDirectory("C:\Users\POOJA\lab work\experiments\Shankhamala\output\ ");
list = getFileList(dirInput);
setBatchMode(true);
for (i=0; i<list.length; i++) {
	
	showProgress(i+1, list.length);
    	open(dirInput+list[i]);
	
    	run("8-bit");
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
