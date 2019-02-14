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
#  [1] https://www.slicer.org                                                         #
#  [2] http://elastix.isi.uu.nl                                                       #
#  [3] Al-Dhamari et al., (2017): ACIR: automatic cochlea image registration.         #
#      In: Proceedings SPIE Medical Imaging 2017: Image Processing;. SPIE. Bd.        #
#          10133. S. 10133p1-10133p5                                                  #
#  [4] https://mtixnat.uni-koblenz.de                                                 #
#                                                                                     #
#  Updated: 8.12.2018                                                                  #    
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

import CochleaSeg
 
#TODO:
# Later:
# - Checking if all above are needed 
# - Cleaning, optimizing, commenting.  
# - Testing in both Windows and Linux. 
# - Supporting DICOM. 
# - Supporting illegal filename.  
# - Using  SlierElastix binaries.   
# - Visualizing the interimediate steps. 
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

  seglogic = CochleaSeg.CochleaSegLogic()

  #----------------------------------------   Initialization 
  def setup(self):
    print(" ")
    print("=======================================================")   
    print("   Automatic Cochlea Image Registration                ")
    print("=======================================================")    

    # to avoid conflict between slicer and elastix ITKs
    #os.environ['ITK_AUTOLOAD_PATH'] = ' '
       
    ScriptedLoadableModuleWidget.setup(self)
      
    # to access logic class functions and setup global variables
  
    # Set default VisSIm location in the user home 
    #TODO: add option user-defined path when installed first time 
    self.logic = CochleaRegLogic()
    self.logic.setGlobalVariables()
    #=================================================================
    #                     Create the GUI interface
    #=================================================================   
    # Create collapsible Button for registration, transformix and invert transform
    self.mainCollapsibleBtn = ctk.ctkCollapsibleButton()
    self.mainCollapsibleBtn.setStyleSheet("ctkCollapsibleButton { background-color: DarkSeaGreen  }")
    self.mainCollapsibleBtn.text = "ACIR: Automatic Cochlea Image Registration"
    self.layout.addWidget(self.mainCollapsibleBtn)
    self.mainFormLayout = qt.QFormLayout(self.mainCollapsibleBtn)

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
    self.mainFormLayout.addRow("Fixed Volume: ", self.fixedSelectorCoBx)

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
    self.mainFormLayout.addRow("Moving Volume: ", self.movingSelectorCoBx)       

    # Create a time label
    self.timeLbl = qt.QLabel("                 Time: 00:00")
    self.timeLbl.setFixedWidth(500)   
    self.tmLbl = self.timeLbl
    
    # Create a textbox for cochlea location
    # TODO activate input IJK values as well
    self.fixedPointEdt = qt.QLineEdit()
    self.fixedPointEdt.setFixedHeight(40)
    Pt = self.logic.fixedPoint
    self.fixedPointEdt.setText(str(Pt))

    # Create a textbox for cochlea location
    # TODO activate input IJK values as well
    self.movingPointEdt = qt.QLineEdit()
    self.movingPointEdt.setFixedHeight(40)
    Pt = self.logic.movingPoint
    self.movingPointEdt.setText(str(Pt))

    # Create a cochlea locator button
    self.fixedFiducialBtn = qt.QPushButton("Pick cochlea location in fixed image    ")
    self.fixedFiducialBtn.setFixedHeight(40)
    self.fixedFiducialBtn.setToolTip("Pick the fixed fiducial point that will be the center of the cropped image")
    self.fixedFiducialBtn.connect('clicked(bool)', lambda: self.onInputFiducialBtnClick("fixed"))
    self.mainFormLayout.addRow( self.fixedFiducialBtn, self.fixedPointEdt)    

    # Create a cochlea locator button
    self.movingFiducialBtn = qt.QPushButton("Pick cochlea location in moving image    ")
    self.movingFiducialBtn.setFixedHeight(40)
    self.movingFiducialBtn.setToolTip("Pick the moving fiducial point that will be the center of the cropped image")
    self.movingFiducialBtn.connect('clicked(bool)', lambda: self.onInputFiducialBtnClick("moving"))
    self.mainFormLayout.addRow( self.movingFiducialBtn, self.movingPointEdt)    
        
    # Add check box for disabling colors in the result of the registration
    self.colorsChkBox = qt.QCheckBox()
    self.colorsChkBox.text = "Disable colors"
    self.colorsChkBox.checked = False
    self.colorsChkBox.stateChanged.connect(self.OnColorsChkBoxChange)
    self.mainFormLayout.addRow(self.colorsChkBox)

    # Create a button to run registration
    self.applyBtn = qt.QPushButton("Run")
    self.applyBtn.setFixedHeight(50)
    self.applyBtn.setFixedWidth (250)
    self.applyBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    self.applyBtn.toolTip = ('How to use:' ' Load at least two images into Slicer. Pick cochlea locations using the buttons and the Slicer Fiducial tool ')
    self.applyBtn.connect('clicked(bool)', self.onApplyBtnClick)
    self.mainFormLayout.addRow(self.applyBtn, self.timeLbl)
    self.runBtn = self.applyBtn
    
    self.layout.addStretch(1) # Collapsible button is held in place when collapsing/expanding.

  #--------------------------------------------------------------------------------------------
  #                        Locating the cochlea: Fiducial placement
  #--------------------------------------------------------------------------------------------               
  def onInputFiducialBtnClick(self, volumeType):
      
    # Remove old Fiducial nodes
    nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for f in nodes:
        if ((f.GetName() == "FixedImageFiducial") or (f.GetName() == "MovingImageFiducial")):
            slicer.mrmlScene.RemoveNode(f)
        #endif
    #endfor
       
    # Create Fiducial Node for the cochlea location in both images
    if (volumeType=="fixed"):
        print(" ..... getting cochlea location in the fixed image")  
        self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")   
        self.logic.locateCochlea(self.fixedSelectorCoBx.currentNode(), self.fixedPointEdt, volumeType)    
        self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    elif (volumeType=="moving"):              
        print(" ..... getting cochlea location in the fixed image")  
        self.movingFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")   
        self.logic.locateCochlea(self.movingSelectorCoBx.currentNode(), self.movingPointEdt, volumeType)    
        self.movingFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    #endif    
    

  # An option to control results displaying
  def OnColorsChkBoxChange(self):
        self.logic.fuseWithOutColor(self.colorsChkBox.checked)
                        
  def onApplyBtnClick(self):
    self.runBtn.setText("...please wait")
    self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
    slicer.app.processEvents()
    self.stm=time.time()
    print("time:" + str(self.stm))
    self.timeLbl.setText("                 Time: 00:00")
    
    # create an option to use IJK point or fidicual node
    self.logic.run( self.fixedSelectorCoBx.currentNode(),self.logic.fixedMarkupNode, self.movingSelectorCoBx.currentNode(),self.logic.movingMarkupNode )
    #self.logic.register( self.fixedSelectorCoBx.currentNode(),self.movingSelectorCoBx.currentNode() )
     
    self.etm=time.time()
    tm=self.etm - self.stm
    self.timeLbl.setText("Time: "+str(tm)+"  seconds")
    self.runBtn.setText("Run")
    self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    slicer.app.processEvents()
  #enddef
  
  def cleanup(self):  
      pass
  #enddef
  
#===================================================================
#                           Logic
#===================================================================
class CochleaRegLogic(ScriptedLoadableModuleLogic):

  segLogic = CochleaSeg.CochleaSegLogic()
  
  #set global paths and parameters
  def setGlobalVariables(self):
    # most of the variables are defined in the segmentation module
    # TODO: create a common class to define common variables and functions 
    self.segLogic.setGlobalVariables()
    self.vissimPath         =  self.segLogic.vissimPath
    self.elastixBinPath     =  self.segLogic.elastixBinPath
    self.transformixBinPath =  self.segLogic.transformixBinPath
    #self.elxInvTransBinPath =  self.segLogic.elxInvTransBinPath
    #self.elastixWebLink     =  self.segLogic.elastixWebLink      
    self.noOutput           =  self.segLogic.noOutput
    self.outputPath         =  self.segLogic.outputPath
    
    self.othersWebLink  =  ("https://cloud.uni-koblenz-landau.de/s/GC82zESbzaDj4dq/download")   
    self.parsPath           = self.vissimPath +"/pars/parCochSeg.txt"
    self.downSz             = 160    
    self.winOS              =0       
    
    # initial poisition = no position
    self.fixedPoint = [0,0,0]
    self.movingPoint = [0,0,0]

    #Cropping Parameters default = 10 
    self.croppingLength = self.segLogic.croppingLength 

    #Resampling parameters, default [0.125, 0.125,0.125]
    self.RSxyz =  self.segLogic.RSxyz

    # color is displayed by default
    self.firstNodeColor  = slicer.modules.colors.logic().GetColorTableNodeID(20)
    self.secondNodeColor = slicer.modules.colors.logic().GetColorTableNodeID(16)

    # windows
    if platform.system()=='Windows':
           #self.elastixWebLink =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/TAc8toxaajSdfy7/download")   
           self.downSz= 500    
    #endif
    #check if VisSimTools folder is found 
    self.checkVisSimTools( )
  #enddef
  
  # Check if image is valid
  def hasImageData(self,inputVolumeNode):
    #check fixed image 
    if not inputVolumeNode:
      logging.debug('hasImageData failed: no input volume node')
      return False
    if inputVolumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in input volume node')
      return False
    return True
  #enddef
  
  def locateCochlea(self, inputVolumeNode,  inputPointEdt, volumeType):

        # Create Fiducial Node for the cochlea location  
        print(" ..... getting cochlea location in "+volumeType+ " image")

        #  Display Coronal during locating the cochlea
        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_cn = green_logic.GetSliceCompositeNode()
        green_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())
        lm = slicer.app.layoutManager()
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpGreenSliceView)
        # Fit slice to window
        sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
        layoutManager = slicer.app.layoutManager()
        for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()
            #endif
        #endfor
               
        # Remove old Fiducial nodes
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        for f in nodes:
            if ((f.GetName() == volumeType+"Point") ):
                    slicer.mrmlScene.RemoveNode(f)
            #endif
        #endfor
  
        # Reset global point label
        inputPoint = [0,0,0]
        inputPointEdt.setText("[0, 0, 0]")
        # Check if a volume is selected
        #
        if not inputVolumeNode:
            print >> sys.stderr, "You need to select a "+ volumeType +" image first."
            return -1
        #endif
        #TODO: reduce this code
        if volumeType=="fixed":
             self.fixedVolumeNode = inputVolumeNode  
             self.fixedPointEdt = inputPointEdt
             self.fixedMarkupNode = slicer.vtkMRMLMarkupsFiducialNode()
             self.fixedMarkupNode.SetName(volumeType+"Point")
             slicer.mrmlScene.AddNode(self.fixedMarkupNode)
             # Start Fiducial Placement Mode in Slicer
             placeModePersistance = 0
             slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)
             # Observe scene for updates
             self.fixedMarkupNode.AddObserver(self.fixedMarkupNode.MarkupAddedEvent,   self.convRAS2IJK)
             self.fixedMarkupNode.AddObserver(self.fixedMarkupNode.MarkupRemovedEvent, self.convRAS2IJK)
             self.fixedMarkupNode.AddObserver(self.fixedMarkupNode.PointModifiedEvent, self.convRAS2IJK)             
        elif volumeType=="moving":
             self.movingVolumeNode = inputVolumeNode  
             self.movingPointEdt = inputPointEdt
             self.movingMarkupNode = slicer.vtkMRMLMarkupsFiducialNode()
             self.movingMarkupNode.SetName(volumeType+"Point")
             slicer.mrmlScene.AddNode(self.movingMarkupNode)
             # Start Fiducial Placement Mode in Slicer
             placeModePersistance = 0
             slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

             # Observe scene for updates
             self.movingMarkupNode.AddObserver(self.movingMarkupNode.MarkupAddedEvent,   self.convRAS2IJK)
             self.movingMarkupNode.AddObserver(self.movingMarkupNode.MarkupRemovedEvent, self.convRAS2IJK)
             self.movingMarkupNode.AddObserver(self.movingMarkupNode.PointModifiedEvent, self.convRAS2IJK)

        #endif
    #enddef

  #--------------------------------------------------------------------------------------------
  #    RAS to  IJK Event
  #--------------------------------------------------------------------------------------------
  # The fiducial point saved in RAS, we need to convert to IJK
  #  more info in our wiki 
  def convRAS2IJK(self, caller, event):
        rasPt = [0,0,0] 
        if (caller.GetName() == "fixedPoint" ):
              vtxt= "fixed"
              ijkIntCoordinates = self.segLogic.ptRAS2IJK(self.fixedMarkupNode, self.fixedVolumeNode)
              self.fixedMarkupNode.GetNthFiducialPosition(0,rasPt)
              self.fixedPoint = ijkIntCoordinates
              self.fixedPointEdt.setText(str(ijkIntCoordinates))
        elif (caller.GetName() == "movingPoint"):
              vtxt= "moving"
              ijkIntCoordinates = self.segLogic.ptRAS2IJK(self.movingMarkupNode, self.movingVolumeNode)
              self.movingMarkupNode.GetNthFiducialPosition(0,rasPt)
              self.movingPoint = ijkIntCoordinates
              self.movingPointEdt.setText(str(ijkIntCoordinates))
        #endif        
        print(" ..... cochlea location RAS: " + str())  
        print(" ..... cochlea location in the "+vtxt+" image set to: " + str(ijkIntCoordinates))  
  #enddef
   
  #===========================================================================================
  #                       Registration Process 
  #--------------------------------------------------------------------------------------------
  # This method perform the registration steps
  def run(self, fixedVolumeNode, fixedFiducialNode, movingVolumeNode, movingFiducialNode):
      # to be used fromoutside we need to do:
      # import CochleaReg
      # logic= CochleaReg.CochleaRegLogic()
      # logic.run(with the parameters above)
      
        """
        Run the actual algorithm
        """
        # we need to run this again in case of external call
        self.setGlobalVariables()
        
        #check if the images are valid
        self.hasImageData(fixedVolumeNode)
        self.hasImageData(movingVolumeNode)
        
        self.fixedVolumeNode = fixedVolumeNode
        self.fixedFiducialNode = fixedFiducialNode
        
        self.movingVolumeNode = movingVolumeNode
        self.movingFiducialNode = movingFiducialNode

        # Create a temporary nodes as workaround for bad path or filename 
        #TODO: create a temp folder and remove temp node before display
        fixedTmpName= self.vissimPath+"/fixedImage.nrrd"
        slicer.util.saveNode( fixedVolumeNode, fixedTmpName)
        [success, self.fixedVolumeNode] = slicer.util.loadVolume(fixedTmpName, returnNode=True)    
        self.fixedVolumeNode.SetName("fixedImage")

        movingTmpName= self.vissimPath+"/movingImage.nrrd"
        slicer.util.saveNode( movingVolumeNode, movingTmpName)
        [success, self.movingVolumeNode] = slicer.util.loadVolume(movingTmpName, returnNode=True)    
        self.movingVolumeNode.SetName("movingImage")

        logging.info('Processing started')

        # Get IJK point from the fiducial to use in cropping  
        self.fixedPoint = self.segLogic.ptRAS2IJK(self.fixedFiducialNode,fixedVolumeNode)
        self.movingPoint = self.segLogic.ptRAS2IJK(self.movingFiducialNode,movingVolumeNode)

        #remove old files if exist
        if os.path.isdir(self.outputPath.strip()): 
           print("removing old output folder!")
           shutil.rmtree(self.outputPath) 
        #endif   
        os.mkdir(self.outputPath)      

        # results paths        
        resTransPath = self.outputPath  + "/TransformParameters.0.txt"
        resInvTransPath  = self.outputPath  + "/TF_Invert.txt"
        resCropPath = self.outputPath + "/result.0.nrrd"          
        resPath = self.outputPath + "/result.nrrd"          
 
        self.fixedPath = self.fixedVolumeNode.GetStorageNode().GetFileName()
        self.fixedFnm  = basename(os.path.splitext(self.fixedPath)[0])   
        
        self.movingPath = self.movingVolumeNode.GetStorageNode().GetFileName()
        self.movingFnm  = basename(os.path.splitext(self.movingPath)[0])   
        
        # Remove old nodes
        rNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in rNodes:
            if f.GetName()[0:3]=='res':
                 slicer.mrmlScene.RemoveNode(f)
            #endif
        #endfor    
        
        # TODO: add better condition
        if ((np.sum(self.fixedPoint)== 0) and (np.sum(self.movingPoint)== 0)) :
            #qt.QMessageBox.critical(slicer.util.mainWindow(),'SlicerCochleaRegistration', 'Cochlea locations are missing')
            print("Error: select cochlea points in fixed and moving images")
            return False
        #endif  

        print("=================== Cropping =====================")                           
        self.fixedCropPath  = self.segLogic.doCropping( self.fixedVolumeNode , self.fixedPoint, self.croppingLength)                     
        self.movingCropPath = self.segLogic.doCropping( self.movingVolumeNode, self.movingPoint, self.croppingLength)                     
        
        print("================= Registration =====================")
        #--------------------------------------- define results paths    -------------------------------
        # transform resulted from regitration of cropped image
        transCropPath     = self.outputPath + "/TransformParameters.0.txt"
        # invert of the above transform 
        transCropInvPath  = self.outputPath + "/TF_Invert.txt"  
        # final transform that transform the input moving image to the input fixed image
        transFinalPath      = self.outputPath + "/TransformParameters.txt"
        # final registered image    
        self.resPath = self.outputPath + "/result.nrrd"          
        #---------------------------------------   registration cropped images      -------------------------------
        # register the cropped images  
        c = self.segLogic.runElastix(self.fixedCropPath,  self.movingCropPath, self.outputPath, self.parsPath, self.noOutput, "496")                               
        #---------------------------------------  transform the moving image     -------------------------------
        # change the size in the transform file to be like the large one                  
        transFinalPath = self.modifyTransformFile(transCropPath,transFinalPath)  
        c = self.segLogic.runTransformix(self.movingPath, self.outputPath, transFinalPath , self.noOutput, "500")
        [success , resultVolumeNode ]=slicer.util.loadVolume(self.resPath, returnNode=True)  
        resultVolumeNode.SetName('result')           
        #remove the temprary cropped files and nodes  
        os.remove(self.fixedCropPath)     
        os.remove(self.movingCropPath)  
        os.remove(fixedTmpName)     
        os.remove(movingTmpName)  
        slicer.mrmlScene.RemoveNode(self.fixedVolumeNode)
        slicer.mrmlScene.RemoveNode(self.movingVolumeNode)
        slicer.mrmlScene.RemoveNode(slicer.util.getNode('fixedImage_crop'))
        slicer.mrmlScene.RemoveNode(slicer.util.getNode('movingImage_crop'))
        slicer.mrmlScene.RemoveNode(slicer.util.getNode('fixedPoint'))
        slicer.mrmlScene.RemoveNode(slicer.util.getNode('movingPoint'))
                                
        #Display the result if no error
        # Clear cochlea location labels
        print("Result saved as: " + self.resPath)
        if  c==0:
            # fixed and registered image are displayed in different colors.
            self.fuseTwoImages(fixedVolumeNode,  resultVolumeNode)
            
            print("==============================================")
            print(" All tasks are done!")
            print("==============================================")        

        else:
            print("error happened during registration, no display ")   
        #endif 
       
  
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
        #endif 
  #enddef

  #--------------------------------------------------------------------------------------------
  #                      Modify Cropped Transform File
  #--------------------------------------------------------------------------------------------
  def modifyTransformFile(self, transCropPath,transFinalPath):       
        # Get spacing, size and origin of the fixed volume
        mFid     = self.fixedVolumeNode.GetImageData()
        mDimensions = mFid.GetDimensions()
        mSpacing = self.fixedVolumeNode.GetSpacing()
        mOrigin  = self.fixedVolumeNode.GetOrigin()

        # TODO check this part 
        # try to solve AX2MR problem         
        mOrigin= (-1 * mOrigin[0] , -1 * mOrigin[1],  mOrigin[2] )

       # Get IJKToRAS direction matrix of the fixed volume and save its values into a list
        matrix = vtk.vtkMatrix4x4()
        self.fixedVolumeNode.GetIJKToRASDirectionMatrix(matrix)
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
        modifiedTransFile = open(transFinalPath, "w+")
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
        return transFinalPath

  def fuseWithOutColor(self, disableColor):
        if not disableColor:
          # Green and Magenta colors 
          self.firstNodeColor  = slicer.modules.colors.logic().GetColorTableNodeID(20)
          self.secondNodeColor = slicer.modules.colors.logic().GetColorTableNodeID(16)
        else:
          self.firstNodeColor  = slicer.modules.colors.logic().GetColorTableNodeID(1)
          self.secondNodeColor = slicer.modules.colors.logic().GetColorTableNodeID(1)
       #endif
  #--------------------------------------------------------------------------------------------
  #                       Display Results
  #--------------------------------------------------------------------------------------------
  # This method displays two images together with different colors which gives fusion effect.
  #TODO: add more options, colors and transparent values in addition to save fused image
  def fuseTwoImages(self, firstNode, secondNode):
        print("=================  Displaying Two Images =====================")

        #TODO: replace this with a loop or short code
        # The volumes are added to the Red (axial) slice of the scene.
        red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
        red_cn = red_logic.GetSliceCompositeNode()
        red_cn.SetBackgroundVolumeID(firstNode.GetID())
        red_cn.SetForegroundVolumeID(secondNode.GetID())
        red_cn.SetForegroundOpacity(0.5)

        # The volumes are added to the Yellow (sagittal) slice of the scene.
        yellow_logic = slicer.app.layoutManager().sliceWidget("Yellow").sliceLogic()
        yellow_cn = yellow_logic.GetSliceCompositeNode()
        yellow_cn.SetBackgroundVolumeID(firstNode.GetID())
        yellow_cn.SetForegroundVolumeID(secondNode.GetID())
        yellow_cn.SetForegroundOpacity(0.5)

        # The volumes are added to the Green (coronal) slice of the scene.
        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_cn = green_logic.GetSliceCompositeNode()
        green_cn.SetBackgroundVolumeID(firstNode.GetID())
        green_cn.SetForegroundVolumeID(secondNode.GetID())
        green_cn.SetForegroundOpacity(0.5)

        # The layout is set to show only the Red slice by default.
        lm = slicer.app.layoutManager()
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)

        # The window level for each image is set to be the same value.
        scan1Node = firstNode.GetScalarVolumeDisplayNode()
        scan2Node = secondNode.GetScalarVolumeDisplayNode()
        scan1Node.AutoWindowLevelOff()
        scan2Node.AutoWindowLevelOff()
        scan1Node.SetWindowLevel(1000, 400)
        scan2Node.SetWindowLevel(1000, 400)

        # Get lookup tables for colors
        self.Nd0 = slicer.modules.colors.logic().GetColorTableNodeID(1)   # original color table
        self.NdG = slicer.modules.colors.logic().GetColorTableNodeID(16)  # green color table
        self.NdM = slicer.modules.colors.logic().GetColorTableNodeID(20)  # magenta color table

        # The lookup tables for each image are applied. Green is assigned to the first Image, 
        # while magenta is assigned to the second  image.
        self.s1DisplayNode = firstNode.GetDisplayNode()
        self.s2DisplayNode = secondNode.GetDisplayNode()
        #resultDisplayNode = resultNode.GetDisplayNode()
        self.s1DisplayNode.SetAndObserveColorNodeID(self.firstNodeColor)
        self.s2DisplayNode.SetAndObserveColorNodeID(self.secondNodeColor)
        #resultDisplayNode.SetAndObserveColorNodeID(magentaNode)

        # Fit slices to window
        sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
        layoutManager = slicer.app.layoutManager()
        for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()
            #endif
        #endfor




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
                   #os.chmod(self.elxInvTransBinPath.strip(),  md)
                #endif 
                msg.setInformativeText("VisSimTools folder is downloaded and ready to use!")
                msg.exec_()                      
                                          
            except Exception as e:
                  print("Error: can not download and extract VisSimTools ...")
                  print(e)   
                  return -1
            #end try-except  
  
#===================================================================
#                           Test
#===================================================================
class CochleaRegTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  logic = CochleaRegLogic()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)
    self.logic.setGlobalVariables()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.testSlicerCochleaRegistration()

  def testSlicerCochleaRegistration(self):

    self.delayDisplay("Starting the test")

    # to get the links from datastore open http://slicer.kitware.com/midas3/community/23 then select a file and click share to get
    # the download link
    # TODO: fix datastore link download problem, the file is created before downloaded   
    #   imgLaWeb = "http://slicer.kitware.com/midas3/download/item/381221/P100001_DV_L_a"
    #   imgLbWeb=  "http://slicer.kitware.com/midas3/download/item/381255/P100001_DV_L_b" 

    fixedFnm   = self.logic.vissimPath +"/imgF.nrrd"
    movingFnm  = self.logic.vissimPath +"/imgM.nrrd"

	# remove testing nodes if found:
    try:
        os.remove(fixedFnm)     
        os.remove(movingFnm)     
    except Exception as e:
				  print("Downloading cochlea sample images ...")

    try:         
        import urllib
        imgLaWebLink = "https://cloud.uni-koblenz-landau.de/s/EwQiQidXqTcGySB/download"
        imgLbWebLink = "https://cloud.uni-koblenz-landau.de/s/qMG2WPjTXabzcbX/download"
        urllib.urlretrieve (imgLaWebLink ,fixedFnm )
        urllib.urlretrieve (imgLbWebLink ,movingFnm )
        print("Downloading complete ...")
 
    except Exception as e:
                  print("Error: can not download sample files  ...")
                  print(e)   
                  return -1
    #end try-except 
    [success, fixedVolumeNode] = slicer.util.loadVolume( fixedFnm, returnNode=True)
    [success, movingVolumeNode] = slicer.util.loadVolume( movingFnm, returnNode=True)
    
    # create a fiducial node for cochlea location for cropping    
    fixedRASpt  = [3.275 , -0.313 , 9.395]
    movingRASpt = [6.089 , 3.438  , 5.365]
    
    fixedFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    fixedFiducialNode.CreateDefaultDisplayNodes()
    fixedFiducialNode.SetName("fixedPoint")  
    fixedFiducialNode.AddFiducial(fixedRASpt[0],fixedRASpt[1],fixedRASpt[2])

    movingFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    movingFiducialNode.CreateDefaultDisplayNodes()
    movingFiducialNode.SetName("movingPoint")  
    movingFiducialNode.AddFiducial(movingRASpt[0],movingRASpt[1],movingRASpt[2])

    # run the registration
    self.logic.run(fixedVolumeNode, fixedFiducialNode, movingVolumeNode, movingFiducialNode)
    self.delayDisplay('Test passed!')
            
