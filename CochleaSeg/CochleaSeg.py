
#======================================================================================
#  3D Slicer [1] plugin that uses elastix toolbox [2] Plugin for Automatic Cochlea    # 
#  Image Segmentation [3]. More info can be found at [4].                             #
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
#  [3] Al-Dhamari et al.,(2018), Automatic Cochlear Length and Volume Size Estimation #
#       First  International Workshop on Context-Aware Operating                      #
#       Theater OR 2, MICCAI 2018, Granada Spain.                                     #
#  [4] https://mtixnat.uni-koblenz.de                                                 #
#                                                                                     #
#  Updated: 6.12.2018                                                                #    
#                                                                                     #  
#======================================================================================

import os, re , datetime, time ,shutil, unittest, logging, zipfile, urllib2, stat,  inspect
import sitkUtils, sys ,math, platform  
import numpy as np, SimpleITK as sitk
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
#   Registration download all stuff for both registration and segmentation.
#   use registration module and commong functions:  

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

class CochleaSeg(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Cochlea Segmentation"
        parent.categories = ["VisSimTools"]
        parent.dependencies = []
        parent.contributors = ["Christopher Guy",
                               "Ibraheem Al-Dhamari",
                               "Michel Peltriauxe",
                               "Anna Gessler",
                               "Jasper Grimmig",
                               "Pepe Eulzer"  
         ]
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        parent.acknowledgementText = " This work is sponsored by Cochlear as part of COMBS project "
        self.parent = parent
  #end def init
#end class CochleaSeg

    
#===================================================================
#                           Main Widget
#===================================================================
class CochleaSegWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    print(" ")
    print("=======================================================")   
    print("   Automatic Cochlea Image Segmentation               ")
    print("=======================================================")           
        
    # to avoid conflict between slicer and elastix ITKs
    os.environ['ITK_AUTOLOAD_PATH'] = ' '

    ScriptedLoadableModuleWidget.setup(self)
    
    # to access logic class functions and setup global variables
    self.logic = CochleaSegLogic()
    # Set default VisSIm location in the user home 
    #TODO: add option user-defined path when installed first time 
    self.logic.setGlobalVariables()
    
    #=================================================================
    #                     Create the GUI interface
    #=================================================================   
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

    # Create a time label
    self.timeLbl = qt.QLabel("                 Time: 00:00")
    self.timeLbl.setFixedWidth(500)   
    self.tmLbl = self.timeLbl
    
    # Create a textbox for cochlea location
    # TODO activate input IJK values as well
    self.inputPointEdt = qt.QLineEdit()
    self.inputPointEdt.setFixedHeight(40)
    Pt = self.logic.inputPoint
    self.inputPointEdt.setText(str(Pt))
    
    # Create a cochlea locator button
    self.inputFiducialBtn = qt.QPushButton("Pick cochlea location in input image    ")
    self.inputFiducialBtn.setFixedHeight(40)
    self.inputFiducialBtn.setToolTip("Pick the input fiducial point that will be the center of the cropped image")
    self.inputFiducialBtn.connect('clicked(bool)', lambda: self.onInputFiducialBtnClick("input"))
    self.mainFormLayout.addRow( self.inputFiducialBtn, self.inputPointEdt)    


    # Create a button to run segmentation
    self.applyBtn = qt.QPushButton("Run")
    self.applyBtn.setFixedHeight(50)
    self.applyBtn.setFixedWidth (250)
    self.applyBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    self.applyBtn.toolTip = ('How to use:' ' Load an images into Slicer. Pick cochlea locations using the buttons and the Slicer Fiducial tool ')
    self.applyBtn.connect('clicked(bool)', self.onApplyBtnClick)
    self.mainFormLayout.addRow(self.applyBtn, self.timeLbl)
    self.runBtn = self.applyBtn
    
    # Add check box for right ear side  
    self.sideChkBox = qt.QCheckBox()
    self.sideChkBox.text = "Right side cochlea"
    self.sideChkBox.stateChanged.connect(self.onSideChkBoxChange)
    self.mainFormLayout.addRow(self.sideChkBox) 
    
    # Create and link Btn to update measuerments    
    self.updateLengthBtn = qt.QPushButton("Update Length")
    self.updateLengthBtn.setFixedHeight(40)
    self.updateLengthBtn.setFixedWidth(250)        
    self.updateLengthBtn.toolTip = ('How to use:' ' Run segmentation first. ')
    self.updateLengthBtn.connect('clicked(bool)', self.onUpdateLengthBtnClick)
    self.mainFormLayout.addRow(self.updateLengthBtn )
  #enddef

    self.layout.addStretch(1) # Collapsible button is held in place when collapsing/expanding.

  def cleanup(self):#nothing to do 
    pass
  #enddef

  # decide which cochlea side to segment
  # TODO: automate the process
  def onSideChkBoxChange(self):
        if  (self.sideChkBox.checked == True):
            self.sideActivated = "R"
        else:
            self.sideActivated = "L"
        #endif    
  #enddef

  # Create a button for updating the cochlea length if the 
  #  the user move some fiducial points
  def onUpdateLengthBtnClick(self):
    self.logic.updateCochleaLength()
    pass
  #enddef

  def onInputFiducialBtnClick(self,volumeType):
    self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")   
    self.logic.locateCochlea(self.inputSelectorCoBx.currentNode(), self.inputPointEdt, volumeType)    
    self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
  #enddef


  def onApplyBtnClick(self):
    self.runBtn.setText("...please wait")
    self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
    slicer.app.processEvents()
    self.stm=time.time()
    print("time:" + str(self.stm))
    self.timeLbl.setText("                 Time: 00:00")
    
    # create an option to use IJK point or fidicual node
    self.logic.run( self.inputSelectorCoBx.currentNode(),self.logic.inputFiducialNode, self.logic.sideActivated )
     
    self.etm=time.time()
    tm=self.etm - self.stm
    self.timeLbl.setText("Time: "+str(tm)+"  seconds")
    self.runBtn.setText("Run")
    self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    slicer.app.processEvents()
  #enddef
  
#===================================================================
#                           Logic
#===================================================================
class CochleaSegLogic(ScriptedLoadableModuleLogic):

  #set global paths and parameters
  def setGlobalVariables(self):
    self.vissimPath    = expanduser("~/VisSimTools")
    self.elastixBinPath    =  self.vissimPath + "/sw/elastix-4.9.0/bin/elastix"
    self.transformixBinPath = self.vissimPath + "/sw/elastix-4.9.0/bin/transformix"
    self.elxInvTransBinPath = self.vissimPath + "/sw/elastix-4.9.0/bin/elxInvertTransform"
    self.elastixWebLink =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/VoxfbJ1kHw0EAQ6/download")      
    self.othersWebLink  =  ("https://mtixnat.uni-koblenz.de/owncloud/index.php/s/TCLlSzwoGK5yX0v/download")   
    self.noOutput= " >> /dev/null"
    self.outputPath = self.vissimPath+"/outputs"
    self.parsPath = self.vissimPath +"/pars/parCochSeg.txt"
    self.modelPath = self.vissimPath +"/models/modelCochlea"        
    self.downSz= 160    
    self.winOS=0       
    self.sideActivated       = "L"    # use left side by default   
    
    # initial poisition = no position
    self.inputPoint = [0,0,0]

    #Cropping Parameters
    self.croppingLength = 10      

    #Resampling parameters 
    self.RSxyz = [0.125, 0.125,0.125]
       
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
    self.checkVisSimTools( )
  #enddef

  # Check if image is valid
  def hasImageData(self,inputVolumeNode):
    #check input image 
    if not inputVolumeNode:
      logging.debug('hasImageData failed: no input volume node')
      return False
    if inputVolumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in input volume node')
      return False
    return True
  #enddef

  def locateCochlea(self, inputVolumeNode,  inputPointEdt, volumeType):

        # redefine to be used in the logic class.
        self.inputVolumeNode   = inputVolumeNode
        self.inputPointEdt     = inputPointEdt  


        # Create Fiducial Node for the cochlea location  
        if (volumeType=="input"):
            print(" ..... getting cochlea location in the input image")  
            # Reset global point label
            self.inputPoint = [0,0,0]
            self.inputPointEdt.setText("[0, 0, 0]")
            # Check if a volume is selected
            if not self.inputVolumeNode:
                print >> sys.stderr, "You need to pick a input volume first."
                return -1
            #endif
        #endif
        
        #  Display Coronal during locating the cochlea
        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_cn = green_logic.GetSliceCompositeNode()
        green_cn.SetBackgroundVolumeID(self.inputVolumeNode.GetID())
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
            if ((f.GetName() == "CochleaLocation") ):
                    slicer.mrmlScene.RemoveNode(f)
            #endif
        #endfor
        
        self.inputFiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
        self.inputFiducialNode.SetName("CochleaLocation")
        slicer.mrmlScene.AddNode(self.inputFiducialNode)

        # Start Fiducial Placement Mode in Slicer
        placeModePersistance = 0
        slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

        # Observe scene for updates
        self.inputFiducialNode.AddObserver(self.inputFiducialNode.MarkupAddedEvent,   self.convRAS2IJK)
        self.inputFiducialNode.AddObserver(self.inputFiducialNode.MarkupRemovedEvent, self.convRAS2IJK)
        self.inputFiducialNode.AddObserver(self.inputFiducialNode.PointModifiedEvent, self.convRAS2IJK)
  #enddef

  #--------------------------------------------------------------------------------------------
  #    RAS to  IJK Event
  #--------------------------------------------------------------------------------------------
  # The fiducial point saved in RAS, we need to convert to IJK
  #  more info in our wiki 
  def convRAS2IJK(self, caller, event):
        ijkIntCoordinates = self.ptRAS2IJK(self.inputFiducialNode, self.inputVolumeNode)
        rasPt = [0,0,0] 
        self.inputFiducialNode.GetNthFiducialPosition(0,rasPt)
        self.inputPoint = ijkIntCoordinates
        self.inputPointEdt.setText(str(ijkIntCoordinates))
        print(" ..... cochlea location RAS: " + str())  
        print(" ..... cochlea location in the fixed image set to: " + str(ijkIntCoordinates))  
  #enddef

 
      
  #--------------------------------------------------------------------------------------------
  #                       Cropping Process  
  #--------------------------------------------------------------------------------------------
  # Using the location as a center point, we cropp around it using the defined cropLength 
  def doCropping(self, inputVolume, point):
        print("================= Begin cropping ... =====================")
        print("Cochlea location: " + str(point))
        
        # resampling spacing  
        self.RSx= self.RSxyz[0] ; self.RSy= self.RSxyz[1];     self.RSz= self.RSxyz[2]

        #Get input image information 
        spacing = inputVolume.GetSpacing()
        imgData = inputVolume.GetImageData()
        dimensions = imgData.GetDimensions()
        
        # compute cropping bounds from image information and cropping parameters
        croppingBounds = [[0,0,0],[0,0,0]];   size = [0,0,0];    lower = [0,0,0] ;     upper = [0,0,0]
        for i in range(0,3):
            size[i] = int((self.croppingLength/spacing[i])/2)
            lower[i] = point[i] - size[i]
            upper[i] = dimensions[i] - (point[i]+size[i])
            # Check if calculated boundaries exceed image dimensions
            if lower[i] < 0:
                    lower[i] = 0
            #endif        
            if upper[i] > dimensions[i]:
                   upper[i] = dimensions[i]
            #endif
        #endfor         
        croppingBounds = [lower,upper]

        # Call SimpleITK CropImageFilter
        print("Cropping with " + str(croppingBounds[0]) + " and " + str(croppingBounds[1]) + ".")
        inputImage = sitkUtils.PullVolumeFromSlicer(inputVolume.GetID())
        cropper = sitkUtils.sitk.CropImageFilter()
        #this generates itk image
        croppedImage = cropper.Execute(inputImage, croppingBounds[0], croppingBounds[1])          
        nodeName = str(inputVolume.GetName()) + "_crop"
        self.inputCropPath = os.path.splitext(inputVolume.GetStorageNode().GetFileName())[0] + "_crop.nrrd"
        # Make a node with cropped image 
        sitkUtils.PushVolumeToSlicer(croppedImage, None, nodeName , 'vtkMRMLScalarVolumeNode' )
        crNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in crNodes:
            if nodeName in f.GetName():
                 f.SetName(nodeName) 
                 break           #endif
        #endfor    

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
        self.inputCropPath = os.path.splitext(inputVolume.GetStorageNode().GetFileName())[0] + "_crop_iso.nrrd"                                    
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
  # This method perform the atlas segementation steps
  def run(self, inputVolumeNode, inputFiducialNode, side):
      #to be used fromoutside we need to do:
      # import CochleaSeg
      # logic= CochleaSeg.CochleaSegLogic()
      # logic.run(with the parameters above)
      
        """
        Run the actual algorithm
        """
        # we need to run this again in case of external call
        self.setGlobalVariables()
        
        self.hasImageData(inputVolumeNode)
        self.inputVolumeNode = inputVolumeNode
        self.inputFiducialNode = inputFiducialNode
        
        # Create a temporary node as workaround for bad path or filename 
        #TODO: create a temp folder and remove temp node before display
        tmpName= self.vissimPath+"/inputImage.nrrd"
        slicer.util.saveNode( inputVolumeNode, tmpName)
        [success, self.inputVolumeNode] = slicer.util.loadVolume(tmpName, returnNode=True)    
        self.inputVolumeNode.SetName("inputImage")

        logging.info('Processing started')

        # Get IJK point from the fiducial to use in cropping  
        self.inputPoint = self.ptRAS2IJK(self.inputFiducialNode,inputVolumeNode)

        #remove old files if exist
        if os.path.isdir(self.outputPath.strip()): 
           print("removing old output folder!")
           shutil.rmtree(self.outputPath) 
        #endif   
        os.mkdir(self.outputPath)      

        # modality type CBCT, CT or MRI
        # it seems using CBCT atlas is enough 
        Styp="Dv"

        # segmentation atlas model paths
        self.modelCropPath       = self.modelPath + "/Mdl"+Styp+side+"c.nrrd" 
        self.modelCropSegPath    = self.modelPath + "/Mdl"+Styp+side+"cS.nrrd" 
        self.modelCropImgStPath  = self.modelPath + "/Mdl"+Styp+side+"cSt.nrrd" 
                
        # results paths        
        resTransPath = self.outputPath  + "/TransformParameters.0.txt"
        res0ImgPath  = self.outputPath  + "/result.0.nrrd"
        resImgPath   = self.outputPath  + "/result.nrrd"
             
        # check if the model is found
        if not isfile(self.modelCropPath): 
            print >> sys.stderr, "ERROR: model is not found"
            return False
        # endif

        self.inputPath = self.inputVolumeNode.GetStorageNode().GetFileName()
        self.inputFnm  = basename(os.path.splitext(self.inputPath)[0])    
        
        # Remove old nodes
        rNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in rNodes:
            if f.GetName()[0:3]=='res':
                 slicer.mrmlScene.RemoveNode(f)
            #endif
        #endfor    
        
        # TODO: add better condition
        if  np.sum(self.inputPoint)== 0 :
            print("Error: select cochlea point")
            return False
        #endif  

        print("=================== Cropping =====================")                           
        self.intputCropPath = self.doCropping(self.inputVolumeNode, self.inputPoint)                     
        print("=================== Segmentation =====================")            
        print ("************  Compute the Transform **********************")
        # register the cropped image  to the model  
        cmd =self.elastixBinPath + " -f " + self.intputCropPath +" -m "+ self.modelCropPath  + " -out " + self.outputPath  + " -p " + self.parsPath + self.noOutput
        print("Executing: " + cmd)
        cTI=os.system(cmd)
        errStr="elastix error at line 546, check the log files"
        self.chkElxER(cTI,errStr) # Check if errors happen during elastix execution
        print ("************  Transform The Segmentation **********************")
        # Apply the transformation to the segmentation image:
        Cmd = self.transformixBinPath + " -in " +self.modelCropSegPath + " -out " + self.outputPath  + " -tp " + resTransPath + self.noOutput
        print("Executing... " + str(Cmd))
        cTS=os.system(Cmd)
        errStr="Transformix error at line 553, check the log file"
        self.chkElxER(cTS,errStr) # Check if errors happen during elastix execution
       
        #rename the result file        
        self.resImgLabelPath = self.outputPath  +"/"+ self.inputVolumeNode.GetName()  +"-label.nrrd"
        self.resultFnm = basename(os.path.splitext(self.resImgLabelPath)[0])                 
        os.rename(resImgPath, self.resImgLabelPath)
                         
        print ("************  Transform The Points **********************")
        # transform the Scala Tympani Points for length Computation 
        Cmd = self.transformixBinPath + " -in " +self.modelCropImgStPath + " -out " + self.outputPath  + " -tp " + resTransPath + self.noOutput
        print("Executing... " + str(Cmd))
        cTP=os.system(Cmd)
        errStr="Transformix error, transform points at line 566, check the log file"
        self.chkElxER(cTP,errStr) # Check if errors happen during elastix execution

        #rename the result image file        
        self.resImgPtsPath = self.outputPath  +"/"+ self.inputVolumeNode.GetName()  + "-IPtsSt.nrrd"
        os.rename(resImgPath, self.resImgPtsPath )       
    
        # Load the image poins, calculat the length and crete the points          
        [success, self.resImgPtsNode] = slicer.util.loadVolume(self.resImgPtsPath, returnNode=True)    
       
        # Convert the image points to fiducials
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
                   self.resultNode = f
                   self.resultNode.SetName('result') 
                #endif
             #endfor      

             # Generate a .seg node
             self.segNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode')
             self.segNode.SetName(self.inputFnm+"S")
             slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(self.resultNode, self.segNode) 
             sg=self.segNode.GetSegmentation()
             st=sg.GetSegment("7")
             st.SetName(self.inputFnm+"_St")
             sv=sg.GetSegment("300")
             sv.SetName(self.inputFnm+"_Sv")

             #remove the temprary loaded labelmap  
             slicer.mrmlScene.RemoveNode(self.resultNode )
        
             # change the model type from vtk to stl 
             msn=slicer.vtkMRMLModelStorageNode()
             msn.SetDefaultWriteFileExtension('stl')
             slicer.mrmlScene.AddDefaultNode(msn)

             # get Cochlea measuerments
             self.getCochleaSize()
        else:
            print("error happened during segmentation ")
        #endif
        
        # TODO: find a reduced code to replace this block  
        lm = slicer.app.layoutManager()
        r_logic = lm.sliceWidget("Red").sliceLogic()
        r_cn = r_logic.GetSliceCompositeNode()
        r_cn.SetBackgroundVolumeID(self.inputVolumeNode.GetID())
        y_logic = lm.sliceWidget("Yellow").sliceLogic()
        y_cn = y_logic.GetSliceCompositeNode()
        y_cn.SetBackgroundVolumeID(self.inputVolumeNode.GetID())
        g_logic = lm.sliceWidget("Green").sliceLogic()
        g_cn = g_logic.GetSliceCompositeNode()
        g_cn.SetBackgroundVolumeID(self.inputVolumeNode.GetID())
        
        print("================= Cochlea analysis is complete  =====================")
        logging.info('Processing completed')
        return True
    #enddef
 
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
        #endif
 #enddef 
  
  #--------------------------------------------------------------------------------------------
  #                        Calculate length and volume of scalas
  #--------------------------------------------------------------------------------------------
  # This function compute the distance between all the fiducials in a markupnode       
  def getFiducilsDistance(self,markupsNode):
        markupsDistance = 0 ; rasPt0 = [0,0,0] ; rasPt1 = [0,0,0] ; 
        
        for i in range(0, markupsNode.GetNumberOfFiducials()):
            markupsNode.GetNthFiducialPosition(i,rasPt1)
            if (i>0) : 
                markupsNode.GetNthFiducialPosition(i-1,rasPt0)
                x0 =rasPt0[0]  ; y0= rasPt0[1]  ; z0=rasPt0[2]                 
                x1 =rasPt1[0]  ; y1= rasPt1[1]  ; z1=rasPt1[2]                 
                markupsDistance = markupsDistance +  math.sqrt( (x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2 )
            #endif
        #endfor    
        return markupsDistance
   
  #--------------------------------------------------------------------------------------------
  #                        Calculate length and volume of scalas
  #--------------------------------------------------------------------------------------------
  def getCochleaSize(self):
        import SegmentStatistics
        print("================= Begin Size Calculation ... =====================")
        volInput = self.inputVolumeNode 
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
            #endif
        #endfor
                    
        # compute the fiducial length
        self.StLength = 0       
        self.StLength = self.getFiducilsDistance(self.markupsNode)

        print(" length = " + str(self.StLength))

        #make a new table  with few columns
        self.resultsTableNode = slicer.vtkMRMLTableNode()
        self.resultsTableNode.SetName(self.inputFnm+"Tbl") #TODO: changet he name to the image+"_t"         
        segStatLogic.exportToTable(self.resultsTableNode)      
        for i in range (0,7):
            self.resultsTableNode.RemoveColumn(3)
        #endfor           
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
        print("Measurements are completed !!! ")
        

  #------------------------------------------------------
  #         image2points
  #------------------------------------------------------
  # This function converted a binary image to fiducials points
  # it is used for computing the length of the scala tympani
  # as skelton methods don't work. 
  # Todo: redesign to be used externally
  #        we need to input a volume node only
  def image2points(self,inputImg):          
        print("====================================================")            
        print("=                Image to Points                  =")      
        print("====================================================") 
        # clone the input image       
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
        #end for 
                
        ptsIJKtmp=  sorted(ptsIJKtmp, key = lambda t: t[-1])     

        self.markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
        self.markupsNode.CreateDefaultDisplayNodes()
        self.markupsNode.SetName("P")  
        mrk=slicer.modules.markups.logic()
        mrk.SetDefaultMarkupsDisplayNodeColor(9)
        mrk.SetDefaultMarkupsDisplayNodeGlyphScale(1.5)   
        mrk.SetDefaultMarkupsDisplayNodeTextScale(0.20)            
        ptsRAS   =np.zeros((NoPts,3))        
        for j in range (0,NoPts):
           ptsIJK[j][0:3] = ptsIJKtmp[j][0:3] 
           ptsRAS[j]=  self.ptIJK2RAS(ptsIJK[j],inputImg)           
           x = ptsRAS[j][0]; y = ptsRAS[j][1] ; z = ptsRAS[j][2]
           self.markupsNode.AddFiducial(x,y,z)
        #endfor        
                
        # get the file name at the first of the plugin
        # this keeps the fiducials with short name
        self.markupsNode.SetName(self.inputFnm + "StPts")
        # Remove image points 
        slicer.mrmlScene.RemoveNode(self.resImgPtsNode)        

#------------------------------------------------------
#                  IJK to RAS  
#------------------------------------------------------
# This function convert an IJK point to RAS point 
#  input:  a point vector and volume node
#  output: a point vector                     
  def ptIJK2RAS(self,ptIJK,inputImg):
        #TODO: add option for printing                   
        # create a IJK2RAs transformation matrix 
        ijk2rasM = vtk.vtkMatrix4x4()
        inputImg.GetIJKToRASMatrix(ijk2rasM)
        ptRAS=np.zeros((len(ptIJK),3))
        ijk= ptIJK
        # create a 4 elements array to get the converted values
        ijkv=[ijk[0],ijk[1],ijk[2],1]             
        rasPt=ijk2rasM.MultiplyDoublePoint(ijkv)
        ptRAS=rasPt[0:3]
        #print("IJK= " + str(ptIJK)+ "   RAS= " + str(ptRAS))
        return  ptRAS       

#------------------------------------------------------
#                 RAS  to IJK 
#------------------------------------------------------   
# This function convert RAS ro an IJK point 
#  input:  a point vector and volume node
#  output: a point vector                     
  def ptRAS2IJK(self,ptRAS,inputImg): 
        #TODO: add option for printing                   
        # create a RAS2IJK transformation matrix 
        ras2ijkM = vtk.vtkMatrix4x4()
        inputImg.GetRASToIJKMatrix(ras2ijkM)       
        ras=[0,0,0]
        ptRAS.GetNthFiducialPosition(0,ras)        
        # create a 4 elements array to get the converted values
        rasv=[ras[0],ras[1],ras[2],1]             
        ptIJKf=np.zeros(3);
        ijkPt=ras2ijkM.MultiplyPoint(rasv)
        ptIJKf[0]=ijkPt[0];ptIJKf [1]=ijkPt[1];ptIJKf [2]=ijkPt[2];
        ptIJK = ptIJKf.astype(np.int64)
        #print("RAS= " + str(ras)+ "   IJK= " + str(ptIJK))
        return  ptIJK       
    
  def updateCochleaLength (self): # input a markup points set                   
        # Convert thepoints from IJK to RAS
        print("Updating cochlea length")
        self.StLength = self.getFiducilsDistance(self.markupsNode)
        print(" length = " + str(self.StLength))      
        self.resultsTableNode.SetCellText(0,2,str(self.StLength))           
        print("Length is updated !!! ")

  def checkVisSimTools(self ):
        # TODO: optimise this part to download only the missing files        
        # Check if elastix exist or download it 
        print(" Defaults paths: " + self.vissimPath)
        print("      VisSimTools folder: " + self.vissimPath)
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
           print("  Parameter file: " + self.parsPath)
           print("  Output folder : " + self.outputPath)            
           print("  Cropping Length: " + str(self.croppingLength))           
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
                   os.chmod(self.elastixBinPath.strip(),  md)
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
class CochleaSegTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  logic = CochleaSegLogic()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)
    self.logic.setGlobalVariables()
    

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.testSlicerCochleaSegmentation()

  def testSlicerCochleaSegmentation(self):

    self.delayDisplay("Starting the test")

    # to get the links from datastore open http://slicer.kitware.com/midas3/community/23 then select a file and click share to get
    # the download link
    # TODO: fix datastore link download problem, the file is created before downloaded   
    #   imgLaWeb = "http://slicer.kitware.com/midas3/download/item/381221/P100001_DV_L_a"

    fnm = self.logic.outputPath +"/imgA.nrrd"
    try:         
        print("Downloading cochlea sample image ...")
        import urllib
        imgLaWebLink = "https://mtixnat.uni-koblenz.de/owncloud/index.php/s/eMvm9LHNHEHoZg3/download"
        #imgLaWebLink = "http://slicer.kitware.com/midas3/download/item/381221/P100001_DV_L_a"
        urllib.urlretrieve (imgLaWebLink ,fnm )
    except Exception as e:
                  print("Error: can not download sample file  ...")
                  print(e)   
                  return -1
    #end try-except 
    [success, inputVolumeNode] = slicer.util.loadVolume( fnm, returnNode=True)
    
    # create a fiducial node for cochlea location for cropping    
    RASpt = [3.275,-0.313,9.395]
    inputFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    inputFiducialNode.CreateDefaultDisplayNodes()
    inputFiducialNode.SetName("cochleaLocationPoint")  
    inputFiducialNode.AddFiducial(RASpt[0],RASpt[1],RASpt[2])
    # use left side cochlea
    side = "L"
    # run the segmentation
    self.logic.run(inputVolumeNode, inputFiducialNode, side)
    self.delayDisplay('Test passed!')
  #enddef

