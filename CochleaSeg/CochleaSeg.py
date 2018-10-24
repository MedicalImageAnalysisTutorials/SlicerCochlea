
#======================================================================================
#  3D Slicer [1] plugin that uses elastix toolbox [2] Plugin for Automatic Cochlea    # 
#  Image Segmentation [3]. More info can be found at [4].                              #
#  Sample cochlea datasets can be downloaded using Slicer Datastore module            #
#                                                                                     #
#  Contributers:                                                                      #   
#      - Christopher L. Guy,   guycl@vcu.edu              : Original source code.     #
#      - Ibraheem Al-Dhamari,  idhamari@uni-koblenz.de    : Plugin design.            #
#      - Michel Peltriaux,     mpeltriaux@uni-koblenz.de  : Programming & testing.    #
#      - Anna Gessler,         agessler@uni-koblenz.de    : Programming & testing.    #
#      - Jasper Grimmig        jgrimmig@uni-koblenz.de    : Programming & testing.    #
#      - Pepe Eulzer           eulzer@uni-koblenz.de      : Programming & testing.    #  
#  [1] http://elastix.isi.uu.nl                                                       #
#  [2] http://elastix.isi.uu.nl                                                       #
#  [3] Al-Dhamari et al.,(2018), Automatic Cochlear Length and Volume Size Estimation #
#       Accepted in: First  International Workshop on Context-Aware Operating         #
#       Theater OR 2, MICCAI 2018, Granada Spain.                                     #
#  [4] https://mtixnat.uni-koblenz.de                                                 #
#                                                                                     #
#  Updated: 7.9.2018                                                                  #    
#                                                                                     #  
#======================================================================================

try:
    # for Python2
    import Tkinter   ## capitalized T in Tkinter 
except ImportError:
    # for Python3
    import tkinter   ## lowercase 't' in tkinter here

import tkFileDialog
import os, re , datetime, time ,shutil, unittest, logging, zipfile, urllib2, stat,  inspect
import sitkUtils, sys ,math, platform  
import  numpy as np, SimpleITK as sitk
import vtkSegmentationCorePython as vtkSegmentationCore
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *   
from copy import deepcopy
from collections import defaultdict
from os.path import expanduser
from os.path import isfile
from os.path import basename
from PythonQt import BoolResult
from shutil import copyfile

#TODO:
# - install requirements automatically     
# - Checking if all above are needed 
# - Cleaning, optimizing, commenting.  
# - Using smaller size binaries or the SlierElastix binaries.   
# - Visualizing the interimediate steps. 
# - Testing in both Windows and Linux. 
# - Supporting DICOM. 
# - Supporting illegal filename.  
# - Use updated paths
# 
#  
# Terminology
#  img         : ITK image 
#  imgNode     : Slicer Node
#  imgPath     : wholePath + Filename
#  imgFnm      : Filename without the path and the extension
#  imgFileName : Filename without the path

#===================================================================
#                           Main Class
#===================================================================
class CochleaSeg(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Cochlea Segmentation"
        parent.categories = ["VisSimTools"]
        parent.dependencies = []
        parent.contributors = [
                               "Christopher Guy",
                               "Ibraheem Al-Dhamari",
                               "Michel Peltriauxe",
                               "Anna Gessler",
                               "Jasper Grimmig",
                               "Pepe Eulzer"  
         ]
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        parent.acknowledgementText = " This work is sponsored by Cochlear as part of COMBS project "
        self.parent = parent

#===================================================================
#                           Main Widget
#=================================================================== 
class CochleaSegWidget(ScriptedLoadableModuleWidget):
    #----------------------------------------   Initialization 
    def setup(self):
        print(" ")
        print("=======================================================")   
        print("   Automatic Cochlea Image Segmentation               ")
        print("=======================================================")           
        
        # to avoid conflict with elastix ITK
        os.environ['ITK_AUTOLOAD_PATH'] = ' '

        ScriptedLoadableModuleWidget.setup(self)

        # Create flags, initial values for checkboxes
        self.sideActivated       = False    # use left side by default   
         
        # Declare global labels   
        self.timeLbl = qt.QLabel("                 Time: 00:00")
        self.timeLbl.setFixedWidth(500)

        # Initialize GUI
        self.initMainPanel()
        #--------------------------------------------------------------------------------------------
        #                       Default Settings  
        #--------------------------------------------------------------------------------------------
        print("   Default Settings: ")        
        # Set default VisSIm location in the user home 
        # It can be a user-defined path when installed first time 
        self.vissimPath    = expanduser("~/VisSimTools")
        print("      VisSimTools folder: " + self.vissimPath)
        # TODO: provide an option to use SlicerElastix binaries        
        self.elastixBinPath    =  self.vissimPath + "/sw/elastix-4.9.0/bin/elastix"
        self.transformixBinPath = self.vissimPath + "/sw/elastix-4.9.0/bin/transformix"
        self.elxInvTransBinPath = self.vissimPath + "/sw/elastix-4.9.0/bin/elxInvertTransform"
        self.elastixWebLink =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/VoxfbJ1kHw0EAQ6/download")      
        self.othersWebLink  =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/TCLlSzwoGK5yX0v/download")   
        self.noOutput= " >> /dev/null"

        # Set default output directory
        self.outputPath = self.vissimPath+"/outputs"
        #self.outputPathLbl.setText("Output dir: " + self.outputPath)
        print("      Output folder : " + self.outputPath)
        
        # Set default parameter file
        self.parsPath = self.vissimPath +"/pars/parCochSeg.txt"
        print("      Parameter file: " + self.parsPath)
        # Set default parameter file
        self.modelPath = self.vissimPath +"/models/modelCochlea"
    
        #radius in mm around cochlea that will be cropped, e.g. 10 mm
        self.croppingLength = 10   
        print("      Cropping Length: " + str(self.croppingLength))
        # download size of VisSimTools folder in MB
        self.downSz= 160    
        self.winOS=0       
        # windows
        if platform.system()=='Windows':
           self.elastixBinPath    = self.elastixBinPath      + ".exe"
           self.transformixBinPath =self.transformixBinPath  + ".exe"
           self.elxInvTransBinPath = self.elxInvTransBinPath + ".exe"  
           self.elastixWebLink =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/TAc8toxaajSdfy7/download")   
           self.noOutput= " > nul"   
           self.winOS=1    
           self.downSz= 500    
        #endif
        
        #check if VisSimTools folder is found 
        self.checkVisSimTools()
        #Resampling parameters 
        self.RSx = 0.125
        self.RSy = 0.125
        self.RSz = 0.125
    #--------------------------------------------------------------------------------------------
    #                        Main  Panel
    #--------------------------------------------------------------------------------------------
    def initMainPanel(self):
        # Create main collapsible Button 
        self.mainCollapsibleBtn = ctk.ctkCollapsibleButton()
        self.mainCollapsibleBtn.setStyleSheet("ctkCollapsibleButton { background-color: DarkSeaGreen  }")
        self.mainCollapsibleBtn.text = "ACIR: Automatic Cochlea Image segmentation"
        self.layout.addWidget(self.mainCollapsibleBtn)
        self.mainFormLayout = qt.QFormLayout(self.mainCollapsibleBtn)
       
        # Create input Volume Selector
        self.inputSelectorCoBx = slicer.qMRMLNodeComboBox()
        self.inputSelectorCoBx.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.inputSelectorCoBx.selectNodeUponCreation = True
        self.inputSelectorCoBx.addEnabled = False
        self.inputSelectorCoBx.removeEnabled = False
        self.inputSelectorCoBx.noneEnabled = False
        self.inputSelectorCoBx.showHidden = False
        self.inputSelectorCoBx.showChildNodeTypes = False
        self.inputSelectorCoBx.setMRMLScene( slicer.mrmlScene )
        self.inputSelectorCoBx.setToolTip("select the input image")
        self.mainFormLayout.addRow("Input image: ", self.inputSelectorCoBx)
     

        # Create Fiducial selector button. The fixed input volume will be cropped around the chosen Fiducial.
        # The lambda function is useful to send which button pressed the function onFiducialButton
        # instead of having two functions 
        self.inputPoint = [0,0,0]
        self.CPos = [0,0,0,0] # location of the cochlea
        self.inputPointEdt = qt.QLineEdit()
        self.inputPointEdt.setFixedHeight(40)

        #self.inputPointEdt.setReadOnly(False)
        self.inputPointEdt.setText(str(self.inputPoint))
        self.inputFiducialBtn = qt.QPushButton("Pick cochlea location in input image    ")
        self.inputFiducialBtn.setFixedHeight(40)
        self.inputFiducialBtn.setToolTip("Pick the input fiducial point that will be the center of the cropped image")
        self.inputFiducialBtn.connect('clicked(bool)', lambda: self.inputFiducialBtnClick("input"))
        self.mainFormLayout.addRow( self.inputFiducialBtn, self.inputPointEdt)    
                     
        # Create and link Button to run segmentation
        self.runBtn = qt.QPushButton("Run")
        self.runBtn.setFixedHeight(50)
        self.runBtn.setFixedWidth (250)
        self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
        self.runBtn.toolTip = ('How to use:' ' Load at least two images into Slicer. Pick cochlea locations using the buttons and the Slicer Fiducial tool ')
        self.runBtn.connect('clicked(bool)', self.runBtnClick)
        self.mainFormLayout.addRow(self.runBtn, self.timeLbl)

        # Add check box for right ear side
        self.sideChkBox = qt.QCheckBox()
        self.sideChkBox.text = "Right side cochlea"
        self.sideChkBox.stateChanged.connect(self.sideChkBoxChange)
        self.mainFormLayout.addRow(self.sideChkBox)
 
        # Add time label
        #self.mainFormLayout.addRow(self.timeLbl)
        self.timeLbl.setText("                 Time: 00:00")
        
        # Create and link Btn to update measuerments 
        self.updateLengthBtn = qt.QPushButton("Update Length")
        self.updateLengthBtn.setFixedHeight(40)
        self.updateLengthBtn.setFixedWidth(250)        
        self.updateLengthBtn.toolTip = ('How to use:' ' Run segmentation first. ')
        self.updateLengthBtn.connect('clicked(bool)', self.updateLengthBtnClick)
        self.mainFormLayout.addRow(self.updateLengthBtn )
    #--------------------------------------------------------------------------------------------
    #                        Locating the cochlea: Fiducial placement
    #--------------------------------------------------------------------------------------------
    def inputFiducialBtnClick(self, volumeType):
        # Create Fiducial Node for the cochlea location  
        if (volumeType=="input"):
            print(" ..... getting cochlea location in the input image")  
            # Reset global point label
            self.inputPoint = [0,0,0]
            self.inputPointEdt.setText("[0, 0, 0]")
            # Check if a volume is selected
            if not self.inputSelectorCoBx.currentNode():
                print >> sys.stderr, "You need to pick a input volume first."
                return -1
            try:  
                self.inputNode = self.inputSelectorCoBx.currentNode()
            except:
                self.suppDICOM()
            #end try 
        
        #  Display Coronal during locating the cochlea
        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_cn = green_logic.GetSliceCompositeNode()
        green_cn.SetBackgroundVolumeID(self.inputNode.GetID())
        lm = slicer.app.layoutManager()
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpGreenSliceView)
        # Fit slice to window
        sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
        layoutManager = slicer.app.layoutManager()
        for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()

        # Remove old Fiducial nodes
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        for f in nodes:
            if ((f.GetName() == "CochleaLocation") ):
                    slicer.mrmlScene.RemoveNode(f)
        self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
        self.inputFiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
        self.inputFiducialNode.SetName("CochleaLocation")
        slicer.mrmlScene.AddNode(self.inputFiducialNode)

        # Start Fiducial Placement Mode in Slicer
        placeModePersistance = 0
        slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

        # Observe scene for updates
        self.inputFiducialNode.AddObserver(self.inputFiducialNode.MarkupAddedEvent, self.convRAS2IJK)
        self.inputFiducialNode.AddObserver(self.inputFiducialNode.MarkupRemovedEvent, self.convRAS2IJK)
        self.inputFiducialNode.AddObserver(self.inputFiducialNode.PointModifiedEvent, self.updateIJK)

        self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")

    #--------------------------------------------------------------------------------------------
    #                     RAS to  IJK  Calculation
    #--------------------------------------------------------------------------------------------
    # The fiducial point saved in RAS, we need to convert to IJK
    #  more info in our wiki 
    def convRAS2IJK(self, caller, event):
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        isInputImage = True
        for f in nodes:
            if (f.GetName() == "CochleaLocation"):
                    inputFiducial = self.inputFiducialNode
                    inputVolume = self.inputSelectorCoBx.currentNode()
         
        print 
        # Get RAS-to-IJK-conversion matrix of the volume node associated with the Fiducial
        RasToIjkMatrix = vtk.vtkMatrix4x4()
        inputVolume.GetRASToIJKMatrix(RasToIjkMatrix)

        # Call function to calculate IJK coordinates
        world = [0,0,0,0]
        inputFiducial.GetNthFiducialWorldCoordinates(0,world)
        ijkDoubleCoordinates = RasToIjkMatrix.MultiplyDoublePoint(world)
        ijkIntCoordinates = [0,0,0]
        for i in range(0,3):
            ijkIntCoordinates[i] = int(ijkDoubleCoordinates[i])

        # Save IJK coordinates in corresponding point variable and display them
        self.CPos = world
        if isInputImage:
            self.inputPoint = ijkIntCoordinates
            self.inputPointEdt.setText(str(ijkIntCoordinates))
            self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")
            print(" ..... cochlea location RAS: " + str(world))  
            print(" ..... cochlea location in the fixed image set to: " + str(ijkIntCoordinates))  

    # in case the user change the location 
    def updateIJK(self, caller, event):
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        isInputImage = True
        for f in nodes:
            if (f.GetName() == "CochleaLocation"):
               inputFiducial = self.inputFiducialNode
               inputVolume = self.inputSelectorCoBx.currentNode()
         
        # Get RAS-to-IJK-conversion matrix of the volume node associated with the Fiducial
        RasToIjkMatrix = vtk.vtkMatrix4x4()
        inputVolume.GetRASToIJKMatrix(RasToIjkMatrix)

        # Call function to calculate IJK coordinates
        world = [0,0,0,0]
        inputFiducial.GetNthFiducialWorldCoordinates(0,world)
        ijkDoubleCoordinates = RasToIjkMatrix.MultiplyDoublePoint(world)
        ijkIntCoordinates = [0,0,0]
        for i in range(0,3):
            ijkIntCoordinates[i] = int(ijkDoubleCoordinates[i])

        # Save IJK coordinates in corresponding point variable and display them
        self.inputPoint = ijkIntCoordinates
        self.inputPointEdt.setText(str(ijkIntCoordinates))
      
    #--------------------------------------------------------------------------------------------
    #                       Cropping Process  
    #--------------------------------------------------------------------------------------------
    # Using the location as a center point, we cropp around it using the defined cropLength 
    def doCropping(self, inputVolume, point):
        print("================= Begin cropping ... =====================")
        print("Cochlea location: " + str(point))

        spacing = inputVolume.GetSpacing()
        imgData = inputVolume.GetImageData()
        dimensions = imgData.GetDimensions()
        croppingBounds = [[0,0,0],[0,0,0]]
        
        size = [0,0,0]
        for i in range(0,3):
            size[i] = int((self.croppingLength/spacing[i])/2)

        # Calculate lower and upper cropping bounds
        lower = [0,0,0]
        for i in range(0,3):
            lower[i] = point[i] - size[i]

        upper = [0,0,0]
        for i in range(0,3):
            upper[i] = dimensions[i] - (point[i]+size[i])

        # Check if calculated boundaries exceed image dimensions
        for i in [lower,upper]:
            for j in range(0,3):
                if i[j] < 0:
                    i[j] = 0
                if i[j] > dimensions[j]:
                    i[j] = dimensions[j]

        croppingBounds = [lower,upper]

        # Call SimpleITK CropImageFilter
        print("Cropping with " + str(croppingBounds[0]) + " and " + str(croppingBounds[1]) + ".")
        inputImage = sitkUtils.PullVolumeFromSlicer(self.inputNode.GetID())
        cropper = sitkUtils.sitk.CropImageFilter()
        #this generates itk image
        croppedImage = cropper.Execute(inputImage, croppingBounds[0], croppingBounds[1])          
        self.inputCropPath = os.path.splitext(self.inputNode.GetStorageNode().GetFileName())[0] + "_crop.nrrd"
        # Make a node with cropped image 
        sitkUtils.PushVolumeToSlicer(croppedImage, None,  inputVolume.GetName() + "_crop", 'vtkMRMLScalarVolumeNode' )

        nodeName = str(inputVolume.GetName()) + "_crop"
        self.croppedNode = slicer.util.getNode(nodeName)
        print("self.inputCropPath : " + str(self.inputCropPath ))
                
        #-------------------------------------------------------
        # Resampling: this produces better looking models  
        #-------------------------------------------------------
       #Run slicer cli module: resample scalar volume
        params = {} 
        params['InputVolume']  = self.croppedNode
        params['OutputVolume'] = self.croppedNode #Resample the cropped image inplace
        params['outputPixelSpacing'] = str(self.RSx) + "," + str(self.RSy) + "," + str(self.RSz) 
        params['interpolationType'] = 'bspline'
        print("....... Resampling")
        slicer.cli.runSync(slicer.modules.resamplescalarvolume, None, params)

        # Save the resulted image to be used in elastix
        properties = {}
        properties["fileType"] = ".nrrd"
        self.inputCropPath = os.path.splitext(self.inputNode.GetStorageNode().GetFileName())[0] + "_crop_iso.nrrd"                                    
        print(" Cropping and resampling are done !!! ")
        # Save cropped image in directory of the original volume
        slicer.util.saveNode( self.croppedNode, self.inputCropPath)
        return self.inputCropPath
     #--------------------------------------------------------------------------------------------
    #                        Extra User Options
    #--------------------------------------------------------------------------------------------      
    def sideChkBoxChange(self):
        if  (self.sideChkBox.checked == True):
            self.sideActivated = True
        else:
            self.sideActivated = False
    #===========================================================================================
    #                       Segmentation Process 
    #--------------------------------------------------------------------------------------------        
    # This method tests the provided inputs before beginning the overlay procedure
    def runBtnClick(self):   
        self.runBtn.setText("...please wait")
        self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
        slicer.app.processEvents()

        # Get the cochlea location point
        sr = self.inputPointEdt.text.strip()
        sr = sr[1:-1]
        self.inputPoint = map(int, sr.split(','))

        #remove old files
        shutil.rmtree(self.outputPath) 
        os.mkdir(self.outputPath)      
         
        Styp="Dv"
        side="L" 

        if self.sideActivated: # right side
              side="R"               
             
        self.modelCropPath       = self.modelPath + "/Mdl"+Styp+side+"c.nrrd" 
        self.modelCropSegPath    = self.modelPath + "/Mdl"+Styp+side+"cS.nrrd" 
        self.modelCropImgStPath  = self.modelPath + "/Mdl"+Styp+side+"cSt.nrrd" 
        
        resTransPath = self.outputPath  + "/TransformParameters.0.txt"
        res0ImgPath  = self.outputPath  + "/result.0.nrrd"
        resImgPath   = self.outputPath  + "/result.nrrd"
             
        # check if the model is found
        if not isfile(self.modelCropPath): 
            print >> sys.stderr, "ERROR: model is not found"
            return False
        # endif

        try:  # is there a node loaded
            self.inputNode = self.inputSelectorCoBx.currentNode() 
        except AttributeError:
            print >> sys.stderr, "ERROR: No input volume, Bad Path or DICOM."
            self.suppDICOM()
        #end try 
        self.inputPath = self.inputNode.GetStorageNode().GetFileName()
        self.inputFnm  = basename(os.path.splitext(self.inputPath)[0])    
        # Remove old result
        rNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in rNodes:
            if f.GetName()[0:3]=='res':
                slicer.mrmlScene.RemoveNode(f)

        self.stm=time.time()
        print("time:" + str(self.stm))
        self.timeLbl.setText("                 Time: 00:00")
   
        if  self.inputPointEdt.text =="[0, 0, 0]" :
            print("Error: select cochlea point")
            return False
        #endif          
        self.intputCropPath = self.doCropping(self.inputNode, self.inputPoint) 
                    
        print("=================== Segmentation =====================")            
        print           
        print ("************  Compute the Transform **********************")
        # register the cropped image  of the model  
        cmd =self.elastixBinPath + " -f " + self.intputCropPath +" -m "+ self.modelCropPath  + " -out " + self.outputPath  + " -p " + self.parsPath + self.noOutput
        print("Executing: " + cmd)
        cTI=os.system(cmd)
        errStr="elastix error at line 601, check the log files"
        self.chkElxER(cTI,errStr) # Check if errors happen during elastix execution
                
        print ("************  Transform The Segmentation **********************")
        # Apply the transformation to the segmentation image:
        Cmd = self.transformixBinPath + " -in " +self.modelCropSegPath + " -out " + self.outputPath  + " -tp " + resTransPath + self.noOutput
        print("Executing... " + str(Cmd))
        cTS=os.system(Cmd)
        errStr="Transformix error at line 611, check the log file"
        self.chkElxER(cTS,errStr) # Check if errors happen during elastix execution
       
        #rename the result file        
        self.resImgLabelPath = self.outputPath  +"/"+ self.inputNode.GetName()  +"-label.nrrd"
        self.resultFnm = basename(os.path.splitext(self.resImgLabelPath)[0])                 
        os.rename(resImgPath, self.resImgLabelPath)
                         
        print ("************  Transform The Points **********************")
        # transform the Scala Tympani Points for length Computation 
        Cmd = self.transformixBinPath + " -in " +self.modelCropImgStPath + " -out " + self.outputPath  + " -tp " + resTransPath + self.noOutput
        print("Executing... " + str(Cmd))
        cTP=os.system(Cmd)
        errStr="Transformix error, transform points at line 623, check the log file"
        self.chkElxER(cTP,errStr) # Check if errors happen during elastix execution

        #rename the result image file        
        self.resImgPtsPath = self.outputPath  +"/"+ self.inputNode.GetName()  + "-IPtsSt.nrrd"
        os.rename(resImgPath, self.resImgPtsPath )       
    
        # Load the image poins, calculat the length and crete the points          
        [success, self.resImgPtsNode] = slicer.util.loadVolume(self.resImgPtsPath, returnNode=True)    
       
        # Convert to points  
        self.image2points(self.resImgPtsNode) 
        
        # Display the result if no error
        # Clear cochlea location labels
        if  (cTS==0) and (cTP==0):
             #remove the temprary cropped file and node  
             os.remove(self.inputCropPath)     
             slicer.mrmlScene.RemoveNode(self.croppedNode )

             # Remove old Fiducial nodes
             nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
             for f in nodes:
                if ((f.GetName() == "CochleaLocation") ):
                    slicer.mrmlScene.RemoveNode(f)   
             
             # The result image file is loaded into the slicer scene.
             [success , self.resultNode]=slicer.util.loadLabelVolume(self.resImgLabelPath,returnNode=True) 
             nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
             for f in nodes:
                match = re.search(r"result_[0-9]*", f.GetName())
                if match:
                   self.resultFnm = f.GetName()
             self.resultNode = slicer.util.getNode(self.resultFnm)
             self.resultNode.SetName('result') 

             # Generate a .seg node
             self.segNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode')
             self.segNode.SetName(self.inputFnm+"S")
             slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(self.resultNode, self.segNode) 
             #self.segNode.SetName(self.inputFnm+"S")  #TODO check if this shoud be removed 
             sg=self.segNode.GetSegmentation()
             st=sg.GetSegment("7")
             st.SetName(self.inputFnm+"_St")
             sv=sg.GetSegment("300")
             sv.SetName(self.inputFnm+"_Sv")

             #remove the temprary labelmap  
             slicer.mrmlScene.RemoveNode(self.resultNode )
        
             # change the model type from vtk to stl 
             msn=slicer.vtkMRMLModelStorageNode()
             msn.SetDefaultWriteFileExtension('stl')
             slicer.mrmlScene.AddDefaultNode(msn)

             # get Cochlea measuerments
             self.getCochleaSize()
        else:
            print("error happened during segmentation ")

        # TODO: find a reduced code to replace this block  
        lm = slicer.app.layoutManager()
        r_logic = lm.sliceWidget("Red").sliceLogic()
        r_cn = r_logic.GetSliceCompositeNode()
        r_cn.SetBackgroundVolumeID(self.inputNode.GetID())
        y_logic = lm.sliceWidget("Yellow").sliceLogic()
        y_cn = y_logic.GetSliceCompositeNode()
        y_cn.SetBackgroundVolumeID(self.inputNode.GetID())
        g_logic = lm.sliceWidget("Green").sliceLogic()
        g_cn = g_logic.GetSliceCompositeNode()
        g_cn.SetBackgroundVolumeID(self.inputNode.GetID())
        
        print("================= Cochlea analysis is complete  =====================")
       
        # Update time label
        self.etm=time.time()
        tm=self.etm - self.stm
        self.timeLbl.setText("Time: "+str(tm)+"  seconds")
        
        self.runBtn.setText("Run")
        self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
        slicer.app.processEvents()
 
    #--------------------------------------------------------------------------------------------
    #                       Check Elastix error
    #--------------------------------------------------------------------------------------------
    # This method checks if errors happen during elastix execution
    def chkElxER(self,c, s):
        if c>0:
           qt.QMessageBox.critical(slicer.util.mainWindow(),'segmentation', s)
           print(s)  
           return False
        else: 
            print("done !!!")

    #--------------------------------------------------------------------------------------------
    #                        Calculate length and volume of scalas
    #--------------------------------------------------------------------------------------------
    def getCochleaSize(self):
        import SegmentStatistics
        # Update time label
        lstm=time.time()

        print("================= Begin Size Calculation ... =====================")
        volInput = self.inputSelectorCoBx.currentNode() 
        volSeg   = slicer.util.getNode((volInput.GetName()+"S"))
        
        segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
        segStatLogic.getParameterNode().SetParameter("Segmentation", volSeg.GetID())
        segStatLogic.getParameterNode().SetParameter("ScalarVolume", volInput.GetID())
        segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.enabled","False")
        segStatLogic.getParameterNode().SetParameter("ScalarVolumeSegmentStatisticsPlugin.voxel_count.enabled","False")
        segStatLogic.computeStatistics()

        # Remove old table
        nodes = slicer.util.getNodesByClass('vtkMRMLTableNode')
        for f in nodes:
            if ((f.GetName() == "CochleaStatistics") ):
                    slicer.mrmlScene.RemoveNode(f)
                    
        #make a new table 
        self.resultsTableNode = slicer.vtkMRMLTableNode()
        self.resultsTableNode.SetName(self.inputFnm+"Tbl") #TODO: changet he name to the image+"_t"         
        segStatLogic.exportToTable(self.resultsTableNode)      

        for i in range (0,7):
            self.resultsTableNode.RemoveColumn(3)
           
        self.resultsTableNode.GetTable().GetColumn(1).SetName("Size (mm^3)")
        self.resultsTableNode.GetTable().GetColumn(2).SetName("Length (mm)")        
        self.resultsTableNode.SetCellText(0,0,"Scala Tympani")
        self.resultsTableNode.SetCellText(1,0,"Scala Vestibuli")
        self.resultsTableNode.SetCellText(0,2,str(self.StLength))           
        self.resultsTableNode.SetCellText(1,2,"0")   
        slicer.mrmlScene.AddNode(self.resultsTableNode)       
        segStatLogic.showTable(self.resultsTableNode)                
        
        # generate and show 3D model
        volSeg.CreateClosedSurfaceRepresentation()
        print("Calculation is done !!! ")
        
    # Convert the transformed image to points
    # TODO: replace when points transform works
#------------------------------------------------------
#         image2points
#------------------------------------------------------
    def image2points(self,inputImg):          

        # to avoid loading system ITK libs
        startTime =time.time()     
        print("====================================================")            
        print("=                Image to Points                  =")      
        print("====================================================") 
        # clone an image      
        #print(inputImg.GetName())   
        nimg = sitkUtils.PullVolumeFromSlicer(inputImg.GetID())
        sz= inputImg.GetImageData().GetDimensions()
        tmpImgArray = slicer.util.array(inputImg.GetID())
        nimgMax = tmpImgArray.max()
        nimgMin = tmpImgArray.min()
          
        b= zip( *np.where(tmpImgArray > (nimgMax/2 ) ) )
        NoPts = len(b)
        ptsIJK   =np.zeros((NoPts,4))        
        ptsIJKtmp=np.zeros((NoPts,4))        
        print("Number of points imported: " + str(NoPts))                  
        for j in range (0,NoPts):
               x= b[j][2] ; y= b[j][1] ; z= b[j][0]            
               ptsIJK[j][0:3] =[ x ,y, z ]
               ptsIJKtmp[j][0:3]   =[ x ,y, z ]
               ptsIJKtmp[j][3]     = tmpImgArray[z][y][x]      
                
        ptsIJKtmp=  sorted(ptsIJKtmp, key = lambda t: t[-1])     
        
        for j in range (0,NoPts):
            ptsIJK[j][0:3] = ptsIJKtmp[j][0:3] 

        print("  Convert the points from IJK to RAS               ")      
        print("---------------------------------------------------")            
        # Convert to RAS
        ptsRAS=  self.ptsIJK2RAS(ptsIJK,inputImg)
        
        print("  Create new points from the input image           ")      
        # create fiducial points 
        if NoPts >= 0 :
           
           mrk=slicer.modules.markups.logic()
           mrk.SetDefaultMarkupsDisplayNodeColor(9)
           mrk.SetDefaultMarkupsDisplayNodeGlyphScale(1.5)   
           mrk.SetDefaultMarkupsDisplayNodeTextScale(0.20)            
           
           self.markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
           self.markupsNode.CreateDefaultDisplayNodes()
           self.markupsNode.SetName("P")  
           for j in range (0, NoPts):
               x = ptsRAS[j][0]
               y = ptsRAS[j][1]
               z = ptsRAS[j][2]
               #print(self.ptsRAS[j][0:3])
               self.markupsNode.AddFiducial(x,y,z)
           #endfor
        #endif
                
        # get the file name at the first of the plugin
        self.markupsNode.SetName(self.inputFnm + "StPts")
        # Remove image points 
        slicer.mrmlScene.RemoveNode(self.resImgPtsNode)        
#------------------------------------------------------
#                  IJK to RAS  
#------------------------------------------------------   
    def ptsIJK2RAS(self,ptsIJK,inputImg): # input a markup points set                   
        # Convert thepoints from IJK to RAS
        ijk2rasM = vtk.vtkMatrix4x4()
        inputImg.GetIJKToRASMatrix(ijk2rasM)
        ptsRAS=np.zeros((len(ptsIJK),3))
        NoPts=len(ptsIJK)
        PtsLength=0 
        for i in range(0, NoPts):  
            ijk= ptsIJK[i][:]
            ijkv=[ijk[0],ijk[1],ijk[2],1]             
            rasPts=ijk2rasM.MultiplyDoublePoint(ijkv)
            ptsRAS[i]=rasPts[0:3]
            if (i>0) : 
                x1 =ptsRAS[i-1][0] ; y1= ptsRAS[i-1][1]; z1=ptsRAS[i-1][2]                
                x2 =rasPts[0]           ; y2= rasPts[1]          ; z2=rasPts[2]                 
                PtsLength =PtsLength +  math.sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
        #print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        self.StLength = PtsLength
        print(" length = " + str(PtsLength))
        #print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        return  ptsRAS       
    
    def updateLengthBtnClick (self): # input a markup points set                   
        # Convert thepoints from IJK to RAS
        
        ptsRASNode =  self.markupsNode
        NoPts= ptsRASNode.GetNumberOfFiducials()
        ptsRAS=np.zeros((NoPts,3))
        PtsLength=0 
        for i in range(NoPts):
            ptsRASNode.GetNthFiducialPosition(i, ptsRAS[i])
            if (i>0) : 
                x1 =ptsRAS[i-1][0] ; y1= ptsRAS[i-1][1]; z1=ptsRAS[i-1][2]                
                x2 =ptsRAS[i][0] ; y2= ptsRAS[i][1]; z2=ptsRAS[i][2]                
                PtsLength =PtsLength +  math.sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
        #print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        self.StLength = PtsLength
        print(" length = " + str(PtsLength))
       
        self.resultsTableNode.GetTable().GetColumn(1).SetName("Size (mm^3)")
        self.resultsTableNode.GetTable().GetColumn(2).SetName("Length (mm)")        
        self.resultsTableNode.SetCellText(0,0,"Scala Tympani")
        self.resultsTableNode.SetCellText(1,0,"Scala Vestibuli")
        self.resultsTableNode.SetCellText(0,2,str(self.StLength))           
        self.resultsTableNode.SetCellText(1,2,"0")   
        slicer.mrmlScene.AddNode(self.resultsTableNode)       
        print("Length is updated !!! ")
        
    def checkVisSimTools(self):
        # TODO: optimise this part to download only the missing files        
        # Check if elastix exist or download it 
        if isfile(self.elastixBinPath.strip()): 
           print("elastix binaries are found in " + self.elastixBinPath )
        else: 
            print("elastix binaries are missing, trying to download ... ")
            msg = qt.QMessageBox()
            msg.setIcon(qt.QMessageBox.Information)
            msg.setText("elastix binaries are missing!")
            msg.setInformativeText("VisSimTools elastix binaries will be downloaded, this may take some time, please wait!")
            msg.setWindowTitle("VisSimTools")
            msg.exec_()
            try:                               
                print("Downloading VisSimTools elastix ...")
                #cmd=" wget --no-check-certificate ""https://mtixnat.uni-koblenz.de/owncloud/index.php/s/3bYztVkSrJxdpDz/download"" -O ~/VisSimToolsTmp.zip"               
                vissimZip = expanduser("~/VisSimToolsTmp.zip")
                with open(vissimZip ,'wb') as f:
                     uFile = urllib2.urlopen(self.elastixWebLink)              
                     chunk = 10024096
                     while 1:
                           data = uFile.read(chunk)
                           f.write(data)                   
                           if not data:
                              f.close()                               
                              print "done!"
                              break
                           #endIf
                           print "Reading ...  %s bytes"%len(data) 
                     #endWhile                               
                print("Extracting to user home ")
                zip_ref = zipfile.ZipFile(vissimZip, 'r')
                zip_ref.extractall(expanduser("~/"))
                zip_ref.close()  
                #remove the downloaded zip file     
                os.remove(vissimZip)                                            
            except Exception as e:
                  print("Error: can not download and extract VisSimTools Elastix ...")
                  print(e)   
                  return -1
            #end try-except 
        #endif
        # check if other files exist
        if isfile(self.parsPath.strip()): 
           print("Other files are found !" )
        else: 
            print("Other files are  missing, trying to download ... ")
            msg = qt.QMessageBox()
            msg.setIcon(qt.QMessageBox.Information)
            msg.setText("Other files are missing!")
            msg.setInformativeText("VisSimTools other files will be downloaded, this may take some time, please wait!")
            msg.setWindowTitle("VisSimTools")
            msg.exec_()
            try:                               
                print("Downloading VisSimTools others ...")
                vissimZip = expanduser("~/VisSimToolsTmp.zip")
                with open(vissimZip ,'wb') as f:
                     uFile = urllib2.urlopen(self.othersWebLink)              
                     chunk = 10024096
                     while 1:
                           data = uFile.read(chunk)
                           f.write(data)                   
                           if not data:
                              f.close()                               
                              print "done!"
                              break
                           #endIf
                           print "Reading ...  %s bytes"%len(data) 
                     #endWhile                               
                print("Extracting to user home ")
                zip_ref = zipfile.ZipFile(vissimZip, 'r')
                zip_ref.extractall(expanduser("~/"))
                zip_ref.close()  
                #remove the downloaded zip file     
                os.remove(vissimZip)   
                # change permission of bin folder for Linux
                if self.winOS==0:   
                   print("Making binaries executable for Linux ")
                   md=  stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH |stat.S_IXGRP |stat.S_IXOTH
                   os.chmod(self.elastixBinPath.strip()    ,  md)
                   os.chmod(self.transformixBinPath.strip(),  md)
                   os.chmod(self.elxInvTransBinPath.strip(),  md)
                #endif 
                msg.setInformativeText("VisSimTools folder is downloaded and ready to use!")
                msg.exec_()                      
                                          
            except Exception as e:
                  print("Error: can not download and extract VisSimTools ...")
                  print(e)   
                  return -1
            #end try-except  
                              
