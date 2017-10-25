/*
 * Tools for SIM STORM Correlation
 * 
 * DESCRIPTION:
 * ------------
 *  This ActionBar enables the registration of SIM/STORM images
 *  (Or any kind of imaging modality where identical bead images 
 *  can be acquired.
 *  Please see the documentation in the publication
 *  https://www.ncbi.nlm.nih.gov/pubmed/28924661
 *  
 * INSTALLATION:
 * -------------
 * 0. Dependencies: This plugin depends from multiple packages
 * 		- ThunderSTORM (If you want to use the shortcuts provided, optional)
 *	 		https://zitmen.github.io/thunderstorm/
 *	 		
 * 		- ActionBar (From the IBMP-CNRS Update Site)
 * 			http://imagejdocu.tudor.lu/doku.php?id=plugin:utilities:action_bar:start
 * 			
 * 		- The PTBIOP Update site
 * 		
 * 1. Copy this file to your Fiji installation under
 *  plugins > ActionBar
 *  eg. C:\Fiji\plugins\ActionBar\BIOP_SIM_STORM_Tools.ijm
 *  
 * 2. After restarting Fiji, you will find it under Plugins > ActionBar
 * 
 */

// Install the BIOP Library
call("BIOP_LibInstaller.installLibrary", "BIOP"+File.separator+"BIOPLib.ijm");

bar_name = "BIOP SIM STORM Tools";
bar_file = replace(bar_name, " ", "_")+".ijm";

runFrom = "/plugins/ActionBar/"+bar_file;

if(isOpen(bar_name)) {
	run("Close AB", bar_name);
}

run("Action Bar",runFrom);
exit();


<codeLibrary>

function toolName() {
	return "SIM STORM Tools";
}

function FWHMParameters(){

	names = newArray("Mean Filter Radius [px]", "Maxima Finder Noise Tolerance");
	types = newArray("n", "n");
	defaults = newArray(1, 10);
	promptParameters(names, types, defaults);

}
/*
 * Measures FWHM in XYZ of for multichannel stacks of beads
 */
function FWHMXYZMeasure() {
	meanRadius = parseFloat(getDataD("Mean Filter Radius [px]", 1)); // Mean Radius
	noise = parseFloat(getDataD("Maxima Finder Noise Tolerance", 10)); // Mean Radius
	name = getTitle();
	getDimensions(nX,nY,nC,nZ,nT);
	getVoxelSize(vx,vy,vz,U);
	
	
	// Suggested by romain, normalize all channels individually, then average them
	// Max projection
	run("Z Project...", "projection=[Max Intensity]");
	zProj = getTitle();
	
	// Suggested by Romain, normalize all channels individually, then average them
	run("32-bit");
	roiManager("Reset");
	cmin = newArray(nSlices);
	cmax = newArray(nSlices);
	
	for (i=1; i<=nSlices; i++) {
		setSlice(i);
		getStatistics(area, mean, tempmin, tempmax);
		cmin[i-1] = tempmin;
		cmax[i-1] = tempmax;
		run("Macro...", "code=v=(v-"+tempmin+")/("+tempmax+"-"+tempmin+")*100");
	}

	// Average them together
	run("Z Project...", "projection=[Average Intensity]");
	avgProj = getTitle();

	
	// Do a happy find maxima
	run("Mean...", "radius="+meanRadius+" stack");
	run("Find Maxima...", "noise="+noise+" exclude output=[Point Selection]");
	getSelectionCoordinates(xpoints,ypoints);
	n = xpoints.length;

	// Filter by intensity

	avgIntFiltered  = newArray(n);	// Instensity values
	xpointsFiltered = newArray(n);	// Kept X positions
	ypointsFiltered = newArray(n);	// Kept Y positions
	
	minSpotI = parseFloat(getDataD("Minimum Relative Instensity [0,1]", 0.7)); // Min spot intensity value


	// Filter the values
	k=0;		// Counter for positions above threshold
	for(a=0; a<n; a++) {
		spotI = getPixel(xpoints[a], ypoints[a]); // Get value, remember that it is averaged before
		if(spotI > minSpotI) {
			avgIntFiltered[k] = spotI;
			xpointsFiltered[k] = xpoints[a];
			ypointsFiltered[k] = ypoints[a];
			k++;
		}
	}
	
	// Finally clean up the arrays
	avgIntFiltered  = Array.trim(avgIntFiltered, k); 
	xpointsFiltered = Array.trim(xpointsFiltered, k); 
	ypointsFiltered = Array.trim(ypointsFiltered, k); 

	makeSelection("point", xpointsFiltered, ypointsFiltered);
	roiManager("Add");
	close(avgProj);

	// Duplicate for Z profile
	selectImage(name);
	//Duplicate and average
	run("Duplicate...", "title=[Averaged sigma "+meanRadius+" - "+name+"] duplicate");
	run("Mean...", "radius="+meanRadius+" stack");
	dup = getTitle();
	
	// Fit 2d Gaussian on Z projection, unfiltered
	selectImage(zProj);
	n = xpointsFiltered.length;
	
	zData = newArray(nZ);
	xes = Array.getSequence(nZ);
	
	// For each point
	for(i=0; i<n; i++) {
		// For Each Color, fit in XY
		for (c=1; c<=nC; c++) {
			selectImage(zProj);
			setSlice(c);
			makeRectangle(xpointsFiltered[i]-5, ypointsFiltered[i]-5,11,11);
			run("Gaussian Fit", "smoothing=1 box_size=2 background=5 min_height=0 fraction_above_background=0 min_width=0 top_n=0 border=2 fit_function=Free fit_background fit_criteria=[Least-squared error] max_iterations=100 significant_digits=4 coord_delta=0.0100 single_fit single_region_size=5 initial_stddev=1.000");
		}
	
	saveAndAppendToResults(name,i+1, nC);
	prepareTable("PSF Analysis");
	r=nResults;
	// Fit in Z again for each color
	for (c=1; c<=nC; c++) {
		selectImage(dup);
		xPos = getResult("Position X", r-c);
		yPos = getResult("Position Y", r-c);
		chan = getResult("Channel", r-c);
		xsd  = getResult("X SD",  r-c);				
		ysd  = getResult("Y SD",  r-c);				
		xyer = getResult("XY Error", r-c);
		
	    //makePoint(xPos,yPos);
	    //roiManager("Add");    
		for(s=1; s<=nZ;s++) {
			Stack.setPosition(chan,s,1);

			zData[s-1] = getPixel(xPos, yPos);
		}
			
		// Now fit in Z
		Fit.doFit("Gaussian", xes, zData);
		Fit.plot();
		rename(name + " - Fit Detection "+(i+1)+" Channel"+c);		
		A	= Fit.p(0); 	// baseLine 
		B 	= Fit.p(1);		// 
		C   = Fit.p(2);		// center
		D 	= Fit.p(3);		// 
		R   = Fit.rSquared;
		// Moments
		skew = getSkew(zData);
		
		setResult("Position X", r-c, xPos*vx);
		setResult("Position Y", r-c,yPos*vy);
		setResult("Position Z", r-c, C*vz);
		setResult("X SD",  r-c,xsd);				
		setResult("Y SD",  r-c,ysd);
		setResult("XY Error", r-c, xyer);
		setResult("Z FWHM", r-c, 2*sqrt(2*log(2))*D*vz);
		setResult("Z rSquared", r-c, R);
		setResult("Skewness", r-c, skew);
		setResult("Distance From FOV Center", r-c, dist(xpointsFiltered[i], ypointsFiltered[i], nX/2, nY/2));

	}
	closeTable("PSF Analysis");
	}

	run("Images to Stack", "name=["+name + "Fitted PSFs] title=[Fit Detection ] use");
	saveCurrentImage();
	// focus on initial image
	selectImage(name);
	
	// Retake all things and calculate delta channels
	prepareTable("PSF Analysis");
	
	deltacol = newArray(n*nC*(nC-1));
	deltaX = newArray(n*nC*(nC-1));
	deltaY = newArray(n*nC*(nC-1));
	deltaZ = newArray(n*nC*(nC-1));
	image  = newArray(n*nC*(nC-1));
	pos    = newArray(n*nC*(nC-1));
	col    = newArray(n*nC*(nC-1));
	det    = newArray(n*nC*(nC-1));
	
	k=0;
	for(i=0; i<n; i++) {
		for(c=1; c<=nC; c++) {
			for(c2=c+1; c2<=nC; c2++) {
				chan = getResultString("Channel", i*nC+c-1);
				img  = getResultString("Image", i*nC+c-1);
				p    = getResultString("Position", i*nC+c-1);
				co   = getResultString("Color", i*nC+c-1);
				de   = getResultString("Detection", i*nC+c-1);
				xPos  = getResult("Position X", i*nC+c-1);
				yPos  = getResult("Position Y", i*nC+c-1);
				zPos  = getResult("Position Z", i*nC+c-1);
				
				chan2 = getResult("Channel", i*nC+c2-1);
				xPos2 = getResult("Position X", i*nC+c2-1);
				yPos2 = getResult("Position Y", i*nC+c2-1);
				zPos2 = getResult("Position Z", i*nC+c2-1);
				
				deltacol[k] = "Delta C"+chan+"-C"+chan2;
				print(i*nC+c-1, chan);

				
				deltaX[k] = xPos - xPos2;
				deltaY[k] = yPos - yPos2;
				deltaZ[k] = zPos - zPos2;
				image[k]  = img;
				pos[k]    = p;
				col[k]    = co;
				det[k]    = de;
				k++;
			}
		}
	}
	
	closeTable("PSF Analysis");
	prepareTable("Delta XYZ Analysis");
	nR=nResults;
	k=0;
	Array.print(image);
	for(i=0; i<n; i++) {
		for(c=1; c<=nC; c++) {
			for(c2=c+1; c2<=nC; c2++) {
				setResult("Point", nR+i, i+1);
				setResult("Image", nR+i, image[k]);
				setResult("Position", nR+i,pos[k]);
				setResult("Color", nR+i, col[k]);
				setResult("Detection", nR+i, det[k]);
				setResult(deltacol[k]+"-X", nR+i, deltaX[k]);
				setResult(deltacol[k]+"-Y", nR+i, deltaY[k]);
				setResult(deltacol[k]+"-Z", nR+i, deltaZ[k]);
				
				k++;
			}
		}
	}
	closeTable("Delta XYZ Analysis");
}

function saveAndAppendToResults(name,d, nChan) {
	
	selectWindow("Fit Results");
	saveAs("Results", "C:\\Fiji.app\\Temp.csv");
	run("Close");
	open("C:\\Fiji.app\\Temp.csv");

	
	// Get the corr ring position for the image
	posnum = substring(name, indexOf(name, "-")+1,lastIndexOf(name,"-"));
	pos = substring(posnum,0,3);
	col = substring(posnum,3);
		
	
	n = nResults;
	if (n==nChan) { // Make sure that we detected something on each channel before using it.
		xpoints = newArray(n);
		ypoints = newArray(n);
		angles  = newArray(n);
		xsd     = newArray(n);
		ysd     = newArray(n);
		error   = newArray(n);
		chan    = newArray(n);
		
		for(i=0; i<n; i++) {
			xpoints[i] = getResult("X",i);
			ypoints[i] = getResult("Y",i);
			angles[i]  = getResult("Angle",i);
			xsd[i]     = getResult("X SD",i);
			ysd[i]     = getResult("Y SD",i);
			error[i]   = getResult("Error",i);
			chan[i]    = i+1;
		}
	
		prepareTable("PSF Analysis");
		r = nResults;
		for(i=0; i<n; i++) {
			setResult("Image", r+i, name);
			setResult("Channel", r+i, chan[i]);
			setResult("Position", r+i, pos);
			setResult("Color", r+i, col);
			setResult("Detection", r+i, d);
			setResult("Angle", r+i,angles[i]);
			setResult("X SD", r+i,xsd[i]);				
			setResult("Y SD", r+i,ysd[i]);
			setResult("XY Error", r+i, error[i]);
			setResult("Position X", r+i,xpoints[i]);
			setResult("Position Y", r+i, ypoints[i]);				
		}
		closeTable("PSF Analysis");
	}
}

function dist(x1,y1,x2,y2) {
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

function getSkew(data) {
	momentN = 0;
	momentD = 0;
	n = data.length;
	Array.getStatistics(data, min, max, mean, stdDev);
	for(i=0; i<n; i++) {
		momentN+= pow(data[i]-mean, 3);
		momentD+= pow(data[i]-mean, 2);
		
	}
	momentN *= 1/n;
	momentD *= 1/(n-1);

	
	
	return sqrt(pow(momentD/momentN,3));
}


function MLSSettings() {
	names = newArray("SIM Gaussian Filter Radius px", "SIM Maxima Finder Noise Tolerance", "STORM Gaussian Filter Radius px", "STORM Maxima Finder Noise Tolerance", "Points Distance Threshold px");
	types = newArray("n", "n", "n", "n", "n");
	defaults = newArray(3, 10, 3, 10, 5);
	promptParameters(names, types, defaults);
}

function MLSGetPoints(sigma_SIM, tolerance_SIM, sigma_STORM, tolerance_STORM, corresp_thr ) {
	
	// Get folder of images
	main_dir = getDirectory("Provide folder containing 'SIM' and 'STORM' acquisition folders");

	sim_dir   = main_dir+File.separator+"SIM"+File.separator;
	storm_dir = main_dir+File.separator+"STORM"+File.separator;
	
	files_sim   = getFileList(sim_dir);
	files_storm = getFileList(storm_dir);
	k=1;

	sim_image = "SIM Stack";
	storm_image = "STORM Stack";

	sim_image_large = "SIM Stack Scaled";
	storm_image_large = "STORM Stack Scaled";
	
	
for (i=0; i< files_sim.length; i++ ) {
	if (isImage(files_sim[i]) && isImage(files_storm[i]) ) {
		
		run("Bio-Formats Windowless Importer", "open=["+sim_dir+files_sim[i]+"]");
		rename("SIM "+k);
		
		run("Bio-Formats Windowless Importer", "open=["+storm_dir+files_storm[i]+"]");
		rename("STORM "+k);
		k++;
		
	}
}

	// Make a stack for each one
	run("Images to Stack", "name=["+sim_image+"] title=SIM use");
	run("Images to Stack", "name=["+storm_image+"] title=STORM use");

	// Prepare for projection SIM
	selectImage(sim_image);
	
	run("Scale...", "x=- y=- width=1024 height=1024 interpolation=Bilinear average process create");
	rename(sim_image_large);
	run("Z Project...", "projection=[Max Intensity]");
	rename("Blurred - "+sim_image+" Projected");

	run("Gaussian Blur...", "sigma="+sigma_SIM);
	
	// Prepare for projection STORM
	selectImage(storm_image);
	
	run("Scale...", "x=- y=- width=1024 height=1024 interpolation=Bilinear average process create");
	rename(storm_image_large);
	
	run("Z Project...", "projection=[Max Intensity]");
	rename("Blurred - "+storm_image+" Projected");

	run("Gaussian Blur...", "sigma="+sigma_STORM);

	// Do SIFT
	affineMat = getDataArray("Initial Affine Transform", ",");

	if(lengthOf(affineMat)< 5) {
		// Run SIFT Target is Blurred STORM Image
		run("Extract SIFT Correspondences", "source_image=[Blurred - "+sim_image+" Projected] target_image=[Blurred - "+storm_image+" Projected] initial_gaussian_blur=2 steps_per_scale_octave=12 minimum_image_size=32 maximum_image_size=512 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.92 filter maximal_alignment_error=50 minimal_inlier_ratio=0.05 minimal_number_of_inliers=7 expected_transformation=Affine");
		
		// Grab the Affine Transform from the Log Window
		affineMat= getMatrixFromLog();
		
		setDataArray("Initial Affine Transform", affineMat, ",");
	} else {
		for(i=0; i<affineMat.length;i++) {
			affineMat[i] = parseFloat(affineMat[i]);
		}
	}
	
	// Now work on the stack and find the nearest neighbors by transforming the SIM coordinates
	selectImage(sim_image);
	getDimensions(dx,dy,dc,dz,dt);

	for(i=0; i<dz;i++) {
		selectImage(sim_image_large);
		setSlice(i+1);
		run("Find Maxima...", "noise="+tolerance_SIM+" output=[Point Selection]");
		getSelectionCoordinates(x_sim, y_sim);
		getSelectionCoordinates(x_simt, y_simt);

		
		selectImage(storm_image_large);
		setSlice(i+1);
		run("Find Maxima...", "noise="+tolerance_STORM+" output=[Point Selection]");
		getSelectionCoordinates(x_storm, y_storm);

		// Transform SIM corrdinates
		applyArray(x_simt,y_simt, affineMat);

		// Find nearest neighbors and append
		idx = getCorrespondences(x_simt, y_simt, x_storm, y_storm, corresp_thr);
		
		// Now save the points in the right order
		orderPoints(idx, x_sim, y_sim, x_storm, y_storm);
	}
	
	
	
}


function getMatrixFromLog() {
	thelog = split(getInfo("log"),"\n");
	str = thelog[thelog.length-1];
	if(startsWith(str, "Estimated transformation model")) {
		arr = substring(str, indexOf(str,"[[")+2, lastIndexOf(str, "]]"));
		arr = split(arr,", ");
		for(i=0; i<arr.length;i++) {
			arr[i] = replace(arr[i],"[\\[\\]]","");
			arr[i] = parseFloat(arr[i]);
		}
		return arr;
	} else {
		return "";
	}
}

function applyArray(x,y, M) {
	for(i=0; i<x.length;i++) {
		xt = (M[0] * x[i]) + (M[1] * y[i]) + M[2];
		yt = (M[3] * x[i]) + (M[4] * y[i]) + M[5];
		x[i] = xt;
		y[i] = yt;
	}
}

function orderPoints(idx, sim_x, sim_y, storm_x, storm_y) {

	// If the Ordered SIM points ROI exists, append to it
	id = findRoiWithName("Ordered SIM Points");

	if (id > -1) {
		roiManager("Select", id);
		getSelectionCoordinates(new_sim_x, new_sim_y);
		
		roiManager("Select", id+1);
		getSelectionCoordinates(new_storm_x, new_storm_y);
		
	} else {
		new_sim_x   = newArray(0);
		new_sim_y   = newArray(0);
		new_storm_x = newArray(0);
		new_storm_y = newArray(0);
	}
	for (i=0; i<idx.length; i++) {
		if (!isNaN(idx[i])) {
			new_sim_x = Array.concat(new_sim_x, sim_x[i]);
			new_sim_y = Array.concat(new_sim_y, sim_y[i]);
			
			new_storm_x = Array.concat(new_storm_x, storm_x[idx[i]]);
			new_storm_y = Array.concat(new_storm_y, storm_y[idx[i]]);
		}
	}

	// Make new ROIs or append
	
	if (id > -1) {
		roiManager("Select", id);
		makeSelection("points", new_sim_x, new_sim_y);
		Roi.setName("Ordered SIM Points");
		roiManager("Update");

		roiManager("Select", id+1);
		makeSelection("points", new_storm_x, new_storm_y);
		Roi.setName("Ordered STORM Points");
		roiManager("Update");
	} else {
		
		makeSelection("points", new_sim_x, new_sim_y);
		Roi.setName("Ordered SIM Points");
		roiManager("Add");
		
		makeSelection("points", new_storm_x, new_storm_y);
		Roi.setName("Ordered STORM Points");
		roiManager("Add");
	}
}
	
function getCorrespondences(x1,y1,x2,y2, thr) {
	
	// Brownian motion setup, just accept for the closest
	// Works OK for low densities
	idx = newArray(0);
	for (i=0; i < x1.length; i++) {
		dMin = 1000;
		ix = NaN;
		for(j=0; j< x2.length; j++) {
			d = dist(x1[i],y1[i], x2[j], y2[j]);
			if ( d < thr && d < dMin ) {
				dMin = d;
				ix = j;
			}
		}
		idx = Array.concat(idx,ix);
	}
	return idx;
}


function reconstructSTORM(storm_file_path) {
	setBatchMode(true);
	// File name for new image name
	storm_file = substring(storm_file_path, lastIndexOf(storm_file_path, File.separator), lengthOf(storm_file_path)-4);
	run("Results... ", "open="+storm_file_path);

	imageSizeX = 42000 ; // nm
	imageSizeY = 42000 ;
	imageSizeZ = 800 ;
	
	voxelSize = 10; //nm
	
	imagePixelNumberX = imageSizeX / voxelSize ;
	imagePixelNumberY = imageSizeY / voxelSize ;
	imagePixelNumberZ = imageSizeZ / voxelSize ;

	storm_image_title = "Reconstructed from "+storm_file;
	newImage(storm_image_title, "32-bit black", imagePixelNumberX, imagePixelNumberY, imagePixelNumberZ);
	run("Properties...", "unit=micron pixel_width="+(voxelSize/1000)+" pixel_height="+(voxelSize/1000)+" voxel_depth="+(voxelSize/1000));

	totalRow = nResults;
	for (rowIndex = 0 ; rowIndex < totalRow ; rowIndex++){
		xRaw 		= getResult("Xc",rowIndex);
		yRaw 		= getResult("Yc",rowIndex);
		zRaw 		= getResult("Zc",rowIndex);
		intensity 	= getResult("I",rowIndex);
		//intensity 	= getResult("I",rowIndex);
		
		xCorr = round (xRaw/voxelSize);
		yCorr = round (yRaw/voxelSize);
		zCorr = round (zRaw/voxelSize);
		//print(xCorr+"-"+yCorr+"-"+rowIndex+"-"+zCorr);
		intensityCorr  = round (intensity/10); 

		print3DCore(xCorr, yCorr, zCorr,intensityCorr);	// only core of the spot
		
	}
	updateDisplay();
	run("Z Project...", "projection=[Max Intensity]");
	run("Brightness/Contrast...");
	// Add Blur if wanted
	
	setBatchMode(false);
}
	

//////////////////////////////////////////////////////////////////////////////// Functions

function print3DCore(xCorr, yCorr, zCorr,intensity){
	Stack.setSlice(zCorr);
	addIntensityToPixel(xCorr, yCorr, intensity);
}

function addIntensityToPixel(xCorr, yCorr, intensity){
	currentIntensity = getPixel(xCorr, yCorr);
	actualIntensity = intensity + currentIntensity;
	setPixel(xCorr, yCorr, actualIntensity);
}


</codeLibrary>

<text><html><font size=2.5 color=#0C2981>Parameters
<line>
<button>
label=Save Parameters
icon=noicon
arg=<macro>
saveParameters();
</macro>

<button>
label=Load Parameters
icon=noicon
arg=<macro>
loadParameters();
</macro>
</line>

<text><html><font size=2.5 color=#0C2981>Beads FWHM XYZ
<line>
<button>
label=Detla Processing Parameters
icon=noicon
arg=<macro>
FWHMParameters();
</macro>
</line>

<line>
<button>
label= current Image
icon=noicon
arg=<macro>
FWHMXYZMeasure();
</macro>
<button>
label= Folder
icon=noicon
arg=<macro>
dir = getDirectory("Gimme");
setData("Image Folder", dir);
setBatchMode(true);
nI = getNumberImages();
for(i=0; i<nI;i++) {
	openImage(i);
	FWHMXYZMeasure();
	run("Close All");
}
setBatchMode(false);
showMessage("DODONEEEEEE");
</macro>
</line>

<text><html><font size=2.5 color=#0C2981>ThunderSTORM Shortcuts
<line>
<button>
label=Astigmatism Calibration
icon=noicon
arg=<macro>
run("Cylindrical lens calibration");
</macro>
</line>

<line>
<button>
label=Run Analysis
icon=noicon
arg=<macro>
run("Run analysis");
</macro>
</line>

<text><html><font size=2.5 color=#0C2981>Warp SIM STORM
<line>
<button>
label=Working Parameters
icon=noicon
arg=<macro>
MLSSettings();
</macro>
</line>

<line>
<button>
label=Test Registration Detection
icon=noicon
arg=<macro>
waitForUser("Open a SIM or STORM bead image");
image_name = getTitle();
if(matches(image_name,".*SIM.*")) {
	blur = getDataD("SIM Gaussian Filter Radius px", 2);
} else {
	blur = getDataD("STORM Gaussian Filter Radius px", 2);
}

run("Scale...", "x=- y=- width=1024 height=1024 interpolation=Bilinear average process create");
run("Gaussian Blur...", "sigma="+blur);
showMessage("After closing this message, test the Find Maxima value that\nwill best match the bead detection");
run("Find Maxima...");
 MLSSettings();
run("Close All");

</macro>

</line>

<line>
<button>
label=Create Registration Points
icon=noicon
arg=<macro>

sigma_SIM       = getData("SIM Gaussian Filter Radius px");
tolerance_SIM   = getData("SIM Maxima Finder Noise Tolerance");
sigma_STORM     = getData("STORM Gaussian Filter Radius px");
tolerance_STORM = getData("STORM Maxima Finder Noise Tolerance");
corresp_thr     = getData("Points Distance Threshold px");

MLSGetPoints( sigma_SIM, tolerance_SIM, sigma_STORM, tolerance_STORM, corresp_thr );
run("Close All");

</macro>

</line>

<line>
<button>
label=Register SIM-STORM
icon=noicon
arg=<macro>
// Get ROIs
if(roiManager("Count") != 2) {
	roiManager("Reset");
	roiFile = File.openDialog("Choose ROI zip File With Matching point ROIs");
	roiManager("Open", roiFile);
}

images = newArray(nImages);
sim_image = "";
storm_image = "";
for(i=0; i<nImages; i++) {
	selectImage(i+1);
	images[i] = getTitle();
	if(matches(images[i], ".*SIM.*"))   { sim_image   = images[i]; }
	if(matches(images[i], ".*STORM.*"))	{ storm_image = images[i]; }
}

if(sim_image == "") sim_image = images[0];
if(storm_image == "") storm_image = images[1];

Dialog.create("Register SIM and STORM");
Dialog.addChoice("SIM Image", images,sim_image);
Dialog.addChoice("STORM Image", images, storm_image);

Dialog.show();

sim_image   = Dialog.getChoice();
storm_image = Dialog.getChoice();

// Adjust the ROIs to the sizes of the images
sim_roi   = findRoiWithName(".*SIM.*");
storm_roi = findRoiWithName(".*STORM.*");

// SIM Image
selectImage(sim_image);
roiManager("Select", sim_roi);
getSelectionCoordinates(x,y);
getDimensions(w,h,c,z,t);
for(i=0; i<x.length; i++) {
	x[i] = x[i] * w/1024;
	y[i] = y[i] * h/1024;
}
makeSelection("points", x,y);

// STORM Image
selectImage(storm_image);
roiManager("Select", storm_roi);
getSelectionCoordinates(x,y);
getDimensions(w,h,c,z,t);
for(i=0; i<x.length; i++) {
	x[i] = x[i] * w/1024;
	y[i] = y[i] * h/1024;
}
makeSelection("points", x,y);


run("Moving Least Squares", "source=["+sim_image+"] target=["+storm_image+"] method=Affine alpha=1.000 gridsize=0.000");
run("Enhance Contrast", "saturated=0.35");
rename("Registered - "+sim_image);
run("Merge Channels...", "c1=[Registered - "+sim_image+"] c2=["+storm_image+"] keep create");
selectImage("Composite");
rename(sim_image+" Registered to "+storm_image);
run("Green");
setSlice(2);
run("Magenta");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", "");
</macro>

</line>
<text><html><font size=2.5 color=#0C2981>Reconstruct STORM
<line>
<button>
label=From NIS-Elements Export
icon=noicon
arg=<macro>
storm_file_path = File.openDialog("txt file with point coordinates");
reconstructSTORM(storm_file_path) 
</macro>
</line>


<line>
<button>
label=Close Others
icon=noicon
arg=<macro>
close("\\Others");
</macro>

</line>
