#@ UpdateService updateService
#@ String (visibility=MESSAGE, value="<html><i>Script written by Yannick Krempp, CIF, 2020</i></html>", required=false) msg0
#@ String (visibility=MESSAGE, value="<html><br/><br/><p><b>Purpose of this script:</b></p><li>1) Identifies nuclei with Stardist2D</li><li>2) Identifies NCC-positive cells defined by <code>NCC_threshold</code></li><li>3) From NCC-positive cells, counts MR-positive cells defined by <code>MR_threshold</code></li></ul></html>", required=false) msg1
#@ Integer (label="Band thickness", value=10, style="slider", min=0, max=15, stepSize=1, description="Thickness of the band around the nuclei in pixels") band_thickness
#@ Integer (label="NCC threshold", value=2000, description="Intensity threshold for NCC") ncc_threshold
#@ Integer (label="MR threshold", value=1000, description="Intensity threshold for MR") mr_threshold
#@ Integer (label="Nucleus area", value=200, description="<html>Minimum area for the nucleus in pixels<sup>2</sup></html>") min_nucleus_area


from ij import IJ
from ij.plugin.frame import RoiManager
from ij.plugin.filter import Analyzer
from ij import WindowManager as wm

RM = RoiManager()        # we create an instance of the RoiManager class
rm = RM.getRoiManager()  # "activate" the RoiManager otherwise it can behave strangely 

GRP_AREA_LARGE_ENOUGH = 0
GRP_NCC_POS = 1
GRP_NCC_NEG = 2
GRP_MR_POS = 3
GRP_MR_NEG = 4

cells_idx = []
mr_pos_idx = []
ncc_pos_idx = []
mr_ncc_pos_idx = []
  
def initialize():
	#Step 0: cleanup and sanity check
	IJ.log("Parameters used:")
	IJ.log("thickness: "+str(band_thickness))
	IJ.log("area: "+str(min_nucleus_area))
	IJ.log("NCC threshold: "+str(ncc_threshold))
	IJ.log("MR threshold: "+str(mr_threshold))
	IJ.log("----")
	
	IJ.log("Initialize started...")
	rm.runCommand("reset")
	IJ.run("Clear Results")
	IJ.run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel")
	IJ.log("Initialize done.")

	#Step1: Getting image information
	ImageTitle = wm.getCurrentImage().getTitle()
	IJ.run("Split Channels");
	#Channel_1 = ImageTitle+" (red)"
	#Channel_2 = ImageTitle+" (green)"
	#Channel_3 = ImageTitle+" (blue)"
	Channel_1 = "C1-"+ImageTitle
	Channel_2 = "C2-"+ImageTitle
	Channel_3 = "C3-"+ImageTitle
	
	wm.getImage(Channel_3).setTitle("MR")
	wm.getImage(Channel_2).setTitle("NCC")
	wm.getImage(Channel_1).setTitle("Nuclei")

def detectNuclei(channel):
	IJ.log("Detecting nuclei with Stardist2D...")
	wm.getImage(channel)
	IJ.run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':"+channel+", 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'Both', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	wm.getImage(channel)
	rm.runCommand("Show All without labels")
	IJ.log("Nuclei detected.")
	IJ.log("Number of objects detected by Stardist: "+str(rm.getCount()))
	IJ.run("Tile", "")


def filter_small_nuclei():

	area_idx = []
	stardist_rois = rm.getCount()
	IJ.run("Set Measurements...", "area mean redirect=NCC decimal=3")
	IJ.selectWindow("Nuclei")
	rm.runCommand("Show All without labels")
	IJ.selectWindow("NCC")
	
	for nuclei in range(0,stardist_rois):
		roi = rm.select(nuclei)
		if IJ.getValue(wm.getImage("NCC"), "Area") > min_nucleus_area:
			area_idx.append(nuclei)
	IJ.log("Number of cells: "+str(len(area_idx)))
	return area_idx

def detect_MR_pos(nuclei):

	IJ.selectWindow("MR")
	IJ.run("Remove Overlay", "")
	IJ.run("Overlay Options...", "stroke=red width=0 fill=red set")
	IJ.run("Set Measurements...", "area mean redirect=MR decimal=3")
	pos_idx = []
	neg_idx = []
	
	for idx in nuclei:
		rm.select(idx)
		MR_mean = IJ.getValue(wm.getImage("MR"), "Mean")
		if MR_mean > mr_threshold:
			pos_idx.append(idx)
			IJ.run(IJ.getImage(), "Add Selection...", "")
		else:
			neg_idx.append(idx)
			
	IJ.run("Show Overlay", "")
	IJ.log("Number of MR+ cells: "+str(len(pos_idx)))
	IJ.log("Number of MR- cells: "+str(len(neg_idx)))
	return pos_idx

def detect_NCC_pos(nuclei, thickness):
			
	IJ.selectWindow("NCC")
	IJ.run("Remove Overlay", "")
	IJ.run("Overlay Options...", "stroke=green width=0 fill=green set")
	IJ.run("Set Measurements...", "area mean redirect=NCC decimal=3");
	pos_idx = []
	neg_idx = []
	
	for idx in nuclei:
		
		ncc_roi = rm.select(idx)
		imp = IJ.getImage()
		IJ.run("Make Band...", "band="+str(thickness))
		#rm.runCommand("Update")
		if IJ.getValue(wm.getImage("NCC"), "Mean") > ncc_threshold:
			pos_idx.append(idx)
			IJ.run(imp, "Add Selection...", "")
			IJ.run(imp, "Select None", "")
		else:
			neg_idx.append(idx)
			
	IJ.run("Show Overlay", "")
	IJ.log("Number of NCC+ cells: "+str(len(pos_idx)))
	IJ.log("Number of NCC- cells: "+str(len(neg_idx)))
	return pos_idx

def displayResults(imp, mr_pos, ncc_pos):
	imp.show()
	mr_ncc_pos_idx = []
	
	IJ.run("Overlay Options...", "stroke=red width=1 fill=red set")
	for idx_mr in mr_pos:
		rm.select(idx_mr)
		IJ.run(imp, "Add Selection...", "")
		
	IJ.run("Overlay Options...", "stroke=green width=2 fill=[]")
	for idx_ncc in ncc_pos:
		rm.select(idx_ncc)
		if idx_ncc in mr_pos:
			IJ.run("Overlay Options...", "stroke=orange width=2 fill=yellow")
			IJ.run(imp, "Add Selection...", "")
			mr_ncc_pos_idx.append(idx_ncc)
		else:
			IJ.run("Overlay Options...", "stroke=green width=2 fill=[]")
			IJ.run(imp, "Add Selection...", "")
		
	IJ.run("Show Overlay", "")
	IJ.log("Number of MR+ AND NCC+ cells: "+str(len(mr_ncc_pos_idx)))

def cleanup():
	cells_idx = []
	mr_pos_idx = []
	ncc_pos_idx = []
	mr_ncc_pos_idx = []
	rm.reset();

if not updateService.getUpdateSite("StarDist").isActive():
	IJ.error("StarDist plugin required ! Please activate the StarDist update site.")
	IJ.exit() 

else:
	
	original_image = wm.getCurrentImage().duplicate()
	original_image.setTitle("Original")
	
	initialize()
	detectNuclei("Nuclei")
	cells_idx = filter_small_nuclei()
	mr_pos_idx = detect_MR_pos(cells_idx)
	ncc_pos_idx = detect_NCC_pos(cells_idx, band_thickness)
	displayResults(original_image, mr_pos_idx, ncc_pos_idx)
	cleanup()



