#======================================================================================
#  3D Slicer [1] plugin that uses elastix toolbox [2] Plugin for Automatic Cochlea    # 
#  Image Registration (ACIR) [3]. More info can be found at [4]                       #
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
#  [3] Al-Dhamari et al., (2017): ACIR: automatic cochlea image registration.         #
#      In: Proceedings SPIE Medical Imaging 2017: Image Processing;. SPIE. Bd.        #
#          10133. S. 10133p1-10133p5                                                  #
#  [4] https://mtixnat.uni-koblenz.de                                                 #
#                                                                                     #
#  Updated: 7.9.2018                                                                  #    
#                                                                                     #  
#======================================================================================

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
class CochleaReg(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Cochlea Registration"
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
        parent.helpText            = " This module uses ACIR method to auatomatically register cochlea images" 
        parent.acknowledgementText = " This work is sponsored by Cochlear as part of COMBS project "
        self.parent = parent

#===================================================================
#                           Main Widget
#===================================================================
class CochleaRegWidget(ScriptedLoadableModuleWidget):

    #----------------------------------------   Initialization 
    def setup(self):
        print(" ")
        print("=======================================================")   
        print("   Automatic Cochlea Image Registration                ")
        print("=======================================================")    

        # to avoid conflict between elastix itk and slicer itk 
        os.environ['ITK_AUTOLOAD_PATH'] = ' '
       
        ScriptedLoadableModuleWidget.setup(self)
      
        # Declare global labels   
        self.timeLbl = qt.QLabel("                 Time: 00:00")
        self.timeLbl.setFixedWidth(500)

        # Initialize GUI
        self.initRegistrationPanel()

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
        self.othersWebLink  =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/x6kts1R4f5RkDcN/download")   
          
        self.noOutput= " >> /dev/null"

        # Set default output directory
        self.outputPath = self.vissimPath+"/outputs"
        #self.outputPathLbl.setText("Output dir: " + self.outputPath)
        print("      Output folder : " + self.outputPath)
        
        # Set default parameter file
        self.parsPath = self.vissimPath +"/pars/parCochReg.txt"
        print("      Parameter file: " + self.parsPath)
    
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
           winOS=1    
           self.downSz= 500    
        #endif
        
        #check if VisSimTools folder is found 
        self.checkVisSimTools()
           

    #--------------------------------------------------------------------------------------------
    #                        Main Panel
    #--------------------------------------------------------------------------------------------
    def initRegistrationPanel(self):
        # Create collapsible Button for registration, transformix and invert transform
        self.registrationCollapsibleBtn = ctk.ctkCollapsibleButton()
        self.registrationCollapsibleBtn.setStyleSheet("ctkCollapsibleButton { background-color: DarkSeaGreen  }")
        self.registrationCollapsibleBtn.text = "ACIR: Automatic Cochlea Image Registration"
        self.layout.addWidget(self.registrationCollapsibleBtn)
        self.registrationFormLayout = qt.QFormLayout(self.registrationCollapsibleBtn)

        # Create fixed Volume Selector
        self.fixedSelectorCoBx                        = slicer.qMRMLNodeComboBox()
        self.fixedSelectorCoBx.nodeTypes              = ["vtkMRMLScalarVolumeNode"]
        self.fixedSelectorCoBx.selectNodeUponCreation = True
        self.fixedSelectorCoBx.addEnabled             = False
        self.fixedSelectorCoBx.removeEnabled          = False
        self.fixedSelectorCoBx.noneEnabled            = False
        self.fixedSelectorCoBx.showHidden             = False
        self.fixedSelectorCoBx.showChildNodeTypes     = False
        self.fixedSelectorCoBx.setMRMLScene( slicer.mrmlScene )
        self.fixedSelectorCoBx.setToolTip("Pick the fixed volume")
        self.registrationFormLayout.addRow("Fixed Volume: ", self.fixedSelectorCoBx)

        # Create moving Volume Selector
        self.movingSelectorCoBx                        = slicer.qMRMLNodeComboBox()
        self.movingSelectorCoBx.nodeTypes              = ["vtkMRMLScalarVolumeNode"]
        self.movingSelectorCoBx.selectNodeUponCreation = True
        self.movingSelectorCoBx.addEnabled             = False
        self.movingSelectorCoBx.removeEnabled          = False
        self.movingSelectorCoBx.noneEnabled            = False
        self.movingSelectorCoBx.showHidden             = False
        self.movingSelectorCoBx.showChildNodeTypes     = False
        self.movingSelectorCoBx.setMRMLScene( slicer.mrmlScene )
        self.movingSelectorCoBx.setToolTip("Pick the moving volume")
        self.registrationFormLayout.addRow("Moving Volume: ", self.movingSelectorCoBx)       

        # Create fixed Fiducial selector. The fixed input volume will be cropped around the chosen Fiducial.
        # The lambda function is useful to send which button pressed the function FiducialButtonClick
        # instead of having two functions 
        self.fixedPoint       = [0,0,0]
        self.fixedPointLbl    = qt.QLineEdit()
        self.fixedPointLbl.setReadOnly(True) # The point can only be edited by placing a new Fiducial
        self.fixedPointLbl.setText(str(self.fixedPoint))
        self.fixedFiducialBtn = qt.QPushButton("Pick cochlea location in fixed image    ")
        #self.fixedFiducialBtn.setFixedWidth(400)
        self.fixedFiducialBtn.setToolTip("Pick the input fiducial point that will be the center of the cropped image")
        self.fixedFiducialBtn.connect('clicked(bool)', lambda: self.FiducialButtonClick("fixed"))
        self.registrationFormLayout.addRow( self.fixedFiducialBtn, self.fixedPointLbl)

        # Create moving Fiducial selector. The moving input volume will be cropped around the chosen Fiducial.
        self.movingPoint       = [0,0,0]
        self.movingPointLbl    = qt.QLineEdit()
        self.movingPointLbl.setReadOnly(True) # The point can only be edited by placing a new Fiducial
        self.movingPointLbl.setText(str(self.movingPoint))
        self.movingFiducialBtn = qt.QPushButton("Pick cochlea location in moving image")
        #self.movingFiducialBtn.setFixedWidth(400)
        self.movingFiducialBtn.setToolTip("Pick the input fiducial point that will be the center of the cropped image")
        self.movingFiducialBtn.connect('clicked(bool)', lambda: self.FiducialButtonClick("moving"))
        self.registrationFormLayout.addRow(self.movingFiducialBtn, self.movingPointLbl)

        # Create and link Button to run the registration
        self.runBtn = qt.QPushButton("Run")
        self.runBtn.setFixedHeight(50)
        self.runBtn.setFixedWidth (250)
        self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
        self.runBtn.toolTip = ('How to use:' ' Load at least two images into Slicer. Pick cochlea locations using the buttons and the Slicer Fiducial tool ')
        self.runBtn.connect('clicked(bool)', self.runBtnClick)
        self.registrationFormLayout.addRow(self.runBtn, self.timeLbl)

        # Add time label
        self.timeLbl.setText("                 Time: 00:00")

        self.registrationFormLayout.addWidget(qt.QLabel("")) # Spacer
        self.layout.addStretch(1) # Collapsible button is held in place when collapsing/expanding.

    #--------------------------------------------------------------------------------------------
    #                        Locating the cochlea: Fiducial placement
    #--------------------------------------------------------------------------------------------               
    def FiducialButtonClick(self, volumeType):
       
        # Remove old Fiducial nodes
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        for f in nodes:
            if ((f.GetName() == "FixedImageFiducial") or (f.GetName() == "MovingImageFiducial")):
                    slicer.mrmlScene.RemoveNode(f)
       
        # Create Fiducial Node for the cochlea location in both images
        if (volumeType=="fixed"):
            print(" ..... getting cochlea location in the fixed image")  
            # Reset global point label
            self.fixedPoint = [0,0,0]
            self.fixedPointLbl.setText("[0, 0, 0]")
            # Check if a volume is selected
            if not self.fixedSelectorCoBx.currentNode():
                print >> sys.stderr, "ERROR: No fixed volume"
                return -1
 
            # TODO bug: if the fiducial placement mode is cancelled using the Slicer Fiducial button, the button remains colored
            self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
            self.movingFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")

            self.fixedFiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
            self.fixedFiducialNode.SetName("FixedImageFiducial")
            slicer.mrmlScene.AddNode(self.fixedFiducialNode)

            self.displayCoronalView(self.fixedSelectorCoBx.currentNode())

            # Start Fiducial Placement Mode in Slicer
            placeModePersistance = 0
            slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

            # Observe scene for updates
            self.fixedFiducialNode.AddObserver(self.fixedFiducialNode.MarkupAddedEvent, self.convRAS2IJK)
            self.fixedFiducialNode.AddObserver(self.fixedFiducialNode.MarkupRemovedEvent, self.convRAS2IJK)
        elif (volumeType=="moving"):
            print(" ..... getting cochlea location in the moving image")  
            # Reset global point label
            self.movingPoint = [0,0,0]
            self.movingPointLbl.setText("[0, 0, 0]")
            # Check if a volume is selected
            if not self.movingSelectorCoBx.currentNode():
                print >> sys.stderr, "ERROR: No moving volume, Bad Path or DICOM."
                return -1
            

            # TODO bug: if the fiducial placement mode is cancelled using the Slicer Fiducial button, the button remains colored
            self.movingFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
            self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")

            self.movingFiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
            self.movingFiducialNode.SetName("MovingImageFiducial")
            slicer.mrmlScene.AddNode(self.movingFiducialNode)

            self.displayCoronalView(self.movingSelectorCoBx.currentNode())

            # Start Fiducial Placement Mode in Slicer
            placeModePersistance = 0
            slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

            # Observe scene for updates
            self.movingFiducialNode.AddObserver(self.movingFiducialNode.MarkupAddedEvent, self.convRAS2IJK)
            self.movingFiducialNode.AddObserver(self.movingFiducialNode.MarkupRemovedEvent, self.convRAS2IJK)

    #--------------------------------------------------------------------------------------------
    #                        Display Coronal during locating the cochlea
    #--------------------------------------------------------------------------------------------
    def displayCoronalView(self, inputVolume):
        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_cn = green_logic.GetSliceCompositeNode()
        green_cn.SetBackgroundVolumeID(inputVolume.GetID())
        lm = slicer.app.layoutManager()
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpGreenSliceView)

        # Fit slice to window
        sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
        layoutManager = slicer.app.layoutManager()
        for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()

    #--------------------------------------------------------------------------------------------
    #                       IJK coordinates from Fiducial 
    #--------------------------------------------------------------------------------------------
    def convRAS2IJK(self, caller, event):
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        isFixedImage = True
        for f in nodes:
            if (f.GetName() == "FixedImageFiducial"):
                    inputFiducial = self.fixedFiducialNode
                    inputVolume = self.fixedSelectorCoBx.currentNode()
            elif (f.GetName() == "MovingImageFiducial"):
                    inputFiducial = self.movingFiducialNode
                    inputVolume = self.movingSelectorCoBx.currentNode()
                    isFixedImage = False

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
        if isFixedImage:
            self.fixedPoint = ijkIntCoordinates
            self.fixedPointLbl.setText(str(ijkIntCoordinates))
            self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")
            print(" ..... cochlea location in the fixed image set to: " + str(ijkIntCoordinates))  
        else:
            self.movingPoint = ijkIntCoordinates
            self.movingPointLbl.setText(str(ijkIntCoordinates))
            self.movingFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")
            print(" ..... cochlea location in the moving image set to: " + str(ijkIntCoordinates))  
            
  
    #--------------------------------------------------------------------------------------------
    #                       Cropping Process  
    #--------------------------------------------------------------------------------------------
    def calculateCroppingBounds(self, dimensions, spacing, point):
        self.croppingLength
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

        return [lower,upper]

    # Returns cropped volume
    def doCropping(self, inputVolume, point):
        print("================= Begin cropping ... =====================")
        print("Cochlea location: " + str(point))

        spacing = inputVolume.GetSpacing()
        imgData = inputVolume.GetImageData()
        dimensions = imgData.GetDimensions()
        croppingBounds = [[0,0,0],[0,0,0]]
        croppingBounds = self.calculateCroppingBounds(dimensions, spacing, point)

        # Call SimpleITK CropImageFilter
        print("Cropping with " + str(croppingBounds[0]) + " and " + str(croppingBounds[1]) + ".")
        inputImage = sitkUtils.PullVolumeFromSlicer(inputVolume.GetID())
        cropper = sitkUtils.sitk.CropImageFilter()
        croppedImage = cropper.Execute(inputImage, croppingBounds[0], croppingBounds[1])

        # Calculate file name and path of cropped image
        storageNode = inputVolume.GetStorageNode()
        pathWithoutExtension = os.path.splitext(storageNode.GetFileName())[0]
        savepath = pathWithoutExtension + "_crop.nrrd"

        # Save cropped image in directory of the original volume
        sitkUtils.PushVolumeToSlicer(croppedImage, None,  inputVolume.GetName() + "_crop",'vtkMRMLScalarVolumeNode')

        nodeName = str(inputVolume.GetName()) + "_crop"
        print("savepath" + str(savepath))
        slicer.util.saveNode(slicer.util.getNode(nodeName), savepath)

        print("================= Cropping finished =====================")
        return savepath

   
    #===========================================================================================
    #                       Registration Process 
    #--------------------------------------------------------------------------------------------
    # This method tests the provided inputs before beginning the overlay procedure
    def runBtnClick(self):
        self.runBtn.setText("...please wait")
        self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
        slicer.app.processEvents()     
        self.runBtn.setText("...please wait")
        self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
        slicer.app.processEvents() 
        self.runBtn.setText("...please wait")
        self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
        slicer.app.processEvents() 
        
        #self.removeOldResults()
        # remove old files in the output folder
        shutil.rmtree(self.outputPath) 
        os.mkdir(self.outputPath)   
        
        try:  # is there a node loaded
            self.fixedNode = self.fixedSelectorCoBx.currentNode() 
        except AttributeError:
                print >> sys.stderr, "ERROR: You need to pick a fixed input volume."
                return False
        try:
            self.fixedPath = self.fixedNode.GetStorageNode().GetFileName()
        except:                    
            self.fixedPath = self.fixedNode.GetStorageNode().GetFileName()
            print(self.fixedPath)
            
        try:
            self.movingNode = self.movingSelectorCoBx.currentNode()                       
        except AttributeError:
            print >> sys.stderr, "ERROR: You need to pick a moving input volume."
            return False
        
        try:
            self.movingPath = self.movingNode.GetStorageNode().GetFileName()
        except:                    
            self.suppDICOM("moving")
            self.movingPath = self.movingNode.GetStorageNode().GetFileName()
            print(self.movingPath)       

        
        # Remove old Nodes 
        rNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in rNodes:
            if f.GetName()[0:3]=='res':
                slicer.mrmlScene.RemoveNode(f)
        
        # Start time  
        self.stm=time.time()
        print("time:" + str(self.stm))
        self.timeLbl.setText("                 Time: 00:00")
                   
        # TODO: optimize this part if possible           
        # here we have 4 situations: 1-nothing  2-crop   3-swp   4-swp+crop 
        print("================= Registration =====================")
        
        transResPath  = self.outputPath + "/TransformParameters.0.txt"
        transInvPath  = self.outputPath + "/TF_Invert.txt"  
        transCropPath = transResPath
        self.resPath = self.outputPath + "/result.nrrd"          
        
        if not ( (self.movingPointLbl.text =="[0, 0, 0]") | (self.movingPointLbl.text =="[0, 0, 0]")) :
                   # return Path of the cropped images  
                   self.fixedCropPath   = self.doCropping(self.fixedNode, self.fixedPoint)
                   self.movingCropPath  = self.doCropping(self.movingNode, self.movingPoint)    

                   # register the cropped images  
                   Cmd = self.elastixBinPath  + " -f " + self.fixedCropPath +" -m "+ self.movingCropPath + " -out " + self.outputPath + " -p " + self.parsPath + self.noOutput
                   print("executing: " + Cmd)
                   c=os.system(Cmd)
                   errStr="elastix error ( crop ) at line 604, check the log file"
                   self.chkElxER(c,errStr) # Check if errors happen during elastix execution   
                               
                   # register the original large moving image using transformix  
                    # change the size in the transform file to be like the large one                  
                   transLargePath = self.modifyTransformFile(transCropPath)
                   
                   Cmd = self.transformixBinPath + " -in " + self.movingPath + " -out " + self.outputPath + " -tp " + transLargePath + self.noOutput
                   print("Executing... " + str(Cmd))
                   c=os.system(Cmd)
                   errStr="Transformix error (cropping) at line 616, check the log file"
                   self.chkElxER(c,errStr) # Check if errors happen during elastix execution 
                   
                   #remove the temprary cropped file and node  
                   os.remove(self.fixedCropPath)     
                   os.remove(self.movingCropPath)  
        else:
                    qt.QMessageBox.critical(slicer.util.mainWindow(),'Registration', 'Cochlea locations are missing')
                    print("Error cochlea locations are missing")
                    return False
        #endif
                                
        #Display the result if no error
        # Clear cochlea location labels
        print("Result saved as: " + self.resPath)
        if  c==0:
            # The images are loaded in Slicer, and the scenes are adjusted.
            self.displayResult(self.fixedPath,  self.movingPath, self.resPath)
        else:
            print("error happened during registration, no display ")   
         
       
        self.runBtn.setText("Run")
        self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    slicer.app.processEvents()     
  
    #--------------------------------------------------------------------------------------------
    #                       Check Elastix error
    #--------------------------------------------------------------------------------------------
    # This method checks if errors happen during elastix execution
    def chkElxER(self,c, s):
        if c>0:
           qt.QMessageBox.critical(slicer.util.mainWindow(),'Registration', s)
           print(s)  
           return False
        else: 
            print("done !!!")

    #--------------------------------------------------------------------------------------------
    #                      Modify Cropped Transform File
    #--------------------------------------------------------------------------------------------
    def modifyTransformFile(self, transCropPath):
        transOutPath = self.outputPath + "/TransformParameters.txt"
        
        # Get spacing, size and origin of the fixed volume
        mFid     = self.fixedNode.GetImageData()
        mDimensions = mFid.GetDimensions()
        mSpacing = self.fixedNode.GetSpacing()
        mOrigin  = self.fixedNode.GetOrigin()

        # TODO check this part 
        # try to solve AX2MR problem         
        mOrigin= (-1 * mOrigin[0] , -1 * mOrigin[1],  mOrigin[2] )

       # Get IJKToRAS direction matrix of the fixed volume and save its values into a list
        matrix = vtk.vtkMatrix4x4()
        self.fixedNode.GetIJKToRASDirectionMatrix(matrix)
        mE = []
        for i in range(0,3):
            for j in range(0,3):
                mE.append(str(matrix.GetElement(i,j)))
        mE[0]= -1 * mE[0]
        mE[5]= -1 * mE[0]

        # Open the cropped transform file and copy it to the new transform file.
        # Replace the lines containing size, spacing, origin and direction matrix.       
        f = open(transCropPath,'r')
        # remove contents
        modifiedTransFile = open(transOutPath, "w+")
        for line in f:
            if line.startswith("(Size"):
                modifiedTransFile.write("(Size " + str(mDimensions[0]) + " " + str(mDimensions[1]) + " " + str(mDimensions[2]) + ")\n")
            elif line.startswith("(Spacing"):
                modifiedTransFile.write("(Spacing " + str(mSpacing[0]) + " " + str(mSpacing[1]) + " " + str(mSpacing[2]) + ")\n")
            elif line.startswith("(Origin"):
                modifiedTransFile.write("(Origin " + str(mOrigin[0]) + " " + str(mOrigin[1]) + " " + str(mOrigin[2]) + ")\n")
            #elif line.startswith("(Direction"):
            #    modifiedTransFile.write("(Direction " + mE[0] + " " + mE[1] + " " + mE[2] + " " + mE[3] + " " + mE[4] + " " + mE[5] + " " + mE[6] + " " + mE[7] + " " + mE[8] + ")\n")
            else:
                modifiedTransFile.write(line) # no need for new line as it is included in the original file
        f.close()
        modifiedTransFile.close()
        return transOutPath

    #--------------------------------------------------------------------------------------------
    #                       Remove Old results and Reset
    #--------------------------------------------------------------------------------------------
    #TODO: add reset button
    # This method removes temporary files and old result files from the output directory and the Slicer scene.
    def removeOldResults(self):
     
        # Remove remaining Fiducial nodes in Slicer
        nodes = slicer.util.getNodesByClass("vtkMRMLMarkupsFiducialNode")
        for f in nodes:
            if ((f.GetName() == "FixedImageFiducial") or (f.GetName() == "MovingImageFiducial")):
                    slicer.mrmlScene.RemoveNode(f)
        # Remove cropped nodes and result nodes in Slicer
        nodes2 = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in nodes2:
            if (f.GetName() == "result"):
                    slicer.mrmlScene.RemoveNode(f)
            if ("_crop" in f.GetName()):
                    slicer.mrmlScene.RemoveNode(f)
        # Reset colors
        try:
            self.s1DisplayNode.SetAndObserveColorNodeID(self.Nd0)
            self.s2DisplayNode.SetAndObserveColorNodeID(self.Nd0)
            fixedNode = slicer.util.getNode(self.fixedName)
            s1 = fixedNode.GetScalarVolumeDisplayNode()
            s1.AutoWindowLevelOn() # TODO Bug - doesn't work?
        except AttributeError:
            pass

        self.timeLbl.setText("Time: 00:00") # Reset time label

    #--------------------------------------------------------------------------------------------
    #                       Display Results
    #--------------------------------------------------------------------------------------------
    # This method loads all registration images and adjusts their visualization properties.
    def displayResult(self, firstScan, secondScan, thirdScan):
        print("=================  Displaying Result =====================")
        # Remove old Fiducial nodes
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        for f in nodes:
            if ((f.GetName() == "FixedImageFiducial") or (f.GetName() == "MovingImageFiducial")):
                    slicer.mrmlScene.RemoveNode(f)
        # The result image file is loaded into the slicer scene.
        [success , resultNode ]=slicer.util.loadVolume(thirdScan, returnNode=True)  
        
        fixedFnm  = basename(os.path.splitext(firstScan)[0])   #(fixedList[len(fixedList) - 1].split('.')[0])
        movingFnm = basename(os.path.splitext(secondScan)[0])  #= (movingList[len(movingList) - 1].split('.')[0])
        resultFnm = basename(os.path.splitext(thirdScan)[0])   #(resultList[len(resultList) - 1].split('.')[0])     
        print(resultFnm)
        # Check if the result volume was loaded under the name 'result_<NUMBER>' by Slicer
        # This can happen when running registration multiple times in a single session
        # TODO The old result nodes are removed in removeOldResults, but this problem still occurs
        nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in nodes:
            match = re.search(r"result_[0-9]*", f.GetName())
            if match:
                resultFnm = f.GetName()
        resultNode = slicer.util.getNode(resultFnm)
        resultNode.SetName('result') 

        # The Slicer nodes for the images are obtained.
        fixedNode  = slicer.util.getNode(fixedFnm)
        movingNode = slicer.util.getNode(movingFnm)
        #resultNode = slicer.util.getNode('result')
 
        # The volumes are added to the Red (axial) slice of the scene.
        red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
        red_cn = red_logic.GetSliceCompositeNode()
        red_cn.SetBackgroundVolumeID(fixedNode.GetID())
        red_cn.SetForegroundVolumeID(resultNode.GetID())
        red_cn.SetForegroundOpacity(0.5)

        # The volumes are added to the Yellow (sagittal) slice of the scene.
        yellow_logic = slicer.app.layoutManager().sliceWidget("Yellow").sliceLogic()
        yellow_cn = yellow_logic.GetSliceCompositeNode()
        yellow_cn.SetBackgroundVolumeID(fixedNode.GetID())
        yellow_cn.SetForegroundVolumeID(resultNode.GetID())
        yellow_cn.SetForegroundOpacity(0.5)

        # The volumes are added to the Green (coronal) slice of the scene.
        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_cn = green_logic.GetSliceCompositeNode()
        green_cn.SetBackgroundVolumeID(fixedNode.GetID())
        green_cn.SetForegroundVolumeID(resultNode.GetID())
        green_cn.SetForegroundOpacity(0.5)

        # The layout is set to show only the Red slice by default.
        lm = slicer.app.layoutManager()
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)

        # The window level for each image is set to be the same value.
        scan1Node = fixedNode.GetScalarVolumeDisplayNode()
        scan2Node = resultNode.GetScalarVolumeDisplayNode()
        scan1Node.AutoWindowLevelOff()
        scan2Node.AutoWindowLevelOff()
        scan1Node.SetWindowLevel(1000, 400)
        scan2Node.SetWindowLevel(1000, 400)

        # Get lookup tables for colors
        self.Nd0 = slicer.modules.colors.logic().GetColorTableNodeID(1)   # original color table
        self.NdG = slicer.modules.colors.logic().GetColorTableNodeID(16)   # green color table
        self.NdM = slicer.modules.colors.logic().GetColorTableNodeID(20) # magenta color table

        # The lookup tables for each image are applied. Green is assigned to the fixed Image, while magenta is assigned to both the moving image and the transformed image.
        magentaNode = slicer.modules.colors.logic().GetColorTableNodeID(20)
        greenNode = slicer.modules.colors.logic().GetColorTableNodeID(16)
        self.s1DisplayNode = fixedNode.GetDisplayNode()
        self.s2DisplayNode = resultNode.GetDisplayNode()
        #resultDisplayNode = resultNode.GetDisplayNode()
        self.s1DisplayNode.SetAndObserveColorNodeID(greenNode)
        self.s2DisplayNode.SetAndObserveColorNodeID(magentaNode)
        #resultDisplayNode.SetAndObserveColorNodeID(magentaNode)

        # Fit slices to window
        sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
        layoutManager = slicer.app.layoutManager()
        for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()

        # Update time label
        self.etm=time.time()
        tm=self.etm - self.stm
        self.timeLbl.setText("Time: "+str(tm)+"  seconds")
        print("==============================================")
        print(" All tasks are done!")
        print("==============================================")        

    # Download VisSimTools folder if not found 
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
