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
#-------------------------------------------------------------------------------------#
#  Slicer 4.11.0                                                                      #
#  Updated: 24.6.2019                                                                 #
#-------------------------------------------------------------------------------------#
#  - Add branches to github to support new Slicer versions                            #                              
#  - Using VisSimCommon for shared functions.                                         #
#  - Use transformation directly to transform the points.                             #
#  - Add more support for windows and mac.                                            #   
#  - Logic functions are independent and can be called from external script.          #
#  - test function can be used in external scripts. A demo example is provided        #
#======================================================================================

from __future__ import print_function
import os, time, unittest, logging
from shutil import copyfile
import numpy as np
from __main__ import qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import SampleData

import VisSimCommon

#TODO:

# 1. solve right model problem:
#    - it seems hardening transform does not work after loading the segmentation. Slicer crashes.
#      with error : Failed to get reference image geometry 
#    - the problem above solved by exporting to label with model as reference then export to .seg.
#      but still we have bad results. More testing is needed.    
# 2. cleaning, removing temp nodes and files  
# 3. add alternative links for models and sample images.
# 4. local testing with different machines. 

# Later:
# - Checking if all above are needed 
# - Cleaning, optimizing, commenting.  
# - Testing in both Windows and Linux. 
# - Supporting DICOM. 
# - Supporting illegal filename.  
# - Add alternative to use elastix binaries directly by downloading the binary release.   
# - Visualizing the interimediate steps. 
# 
#  
#  
# Terminology
#  img         : ITK image 
#  imgNode     : Slicer Node
#  imgName     :  Filename without the path and without extension
#  imgPath     : wholePath + Filename and extension

#===================================================================
#                           Main Class
#===================================================================

class CochleaSeg(ScriptedLoadableModule):
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


  def setup(self):
    print(" ")
    print("=======================================================")   
    print("   Automatic Cochlea Image Segmentation               ")
    print("=======================================================")           

    ScriptedLoadableModuleWidget.setup(self)
    
    # to access logic class functions and setup global variables
    self.logic = CochleaSegLogic()

    # Set default VisSIm location in the user home 
    #TODO: add option user-defined path when installed first time 
    self.vsc   = VisSimCommon.VisSimCommonLogic()   
    self.vsc.setGlobalVariables(0)
    
    #----------------------------------------------------------------- 
    #                     Create the GUI interface
    #----------------------------------------------------------------- 
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
    self.inputPointEdt.setText("[0,0,0]")
    
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
  #------------------------------------------------------------------------
  #                        Define GUI Elements Functions
  #------------------------------------------------------------------------
  # decide which cochlea side to segment
  # TODO: automate the process
  def onSideChkBoxChange(self):
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      self.vsc.setItemChk("cochleaSide", self.sideChkBox.checked, "cochleaSide", nodes)
  #enddef

  # Create a button for updating the cochlea length if the 
  #  the user move some fiducial points
  def onUpdateLengthBtnClick(self):
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if ( "_StPts" in f.GetName()) :
             vtImgStNode = f  
             break   
          #endif
      #endfor      
  
      nodes = slicer.util.getNodesByClass('vtkMRMLTableNode')
      for f in nodes:
          if ( "_tbl" in f.GetName() ):
             spTblNode = f  
             break   
          #endif
      #endfor      
      self.vsc.getFiducilsDistance(vtImgStNode, spTblNode)
  #enddef

  def onInputFiducialBtnClick(self,volumeType):       
      self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")   
 
      self.inputVolumeNode=self.inputSelectorCoBx.currentNode()
      self.logic.inputFiducialNode = None    
      #remove old nodes

      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if ((f.GetName() == self.inputVolumeNode.GetName()+"_CochleaLocation") ):
             #replace  current 
             #print("inputFiducialNode exist")
             self.logic.inputFiducialNode = f  
             newNode= False
            #endif
      #endfor    
      if not hasattr(self.vsc, 'vtVars'):
         self.vsc.setGlobalVariables(0)
      #end 
      self.vsc.locateItem(self.inputSelectorCoBx.currentNode(), self.inputPointEdt, 0 , 0)    
      self.logic.inputFiducialNode= self.vsc.inputFiducialNodes[0]

      self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
  #enddef


  def onApplyBtnClick(self):
      self.runBtn.setText("...please wait")
      self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
      slicer.app.processEvents()
      self.stm=time.time()
      print("time:" + str(self.stm))
      self.timeLbl.setText("                 Time: 00:00")
    
      self.logic.run( self.inputSelectorCoBx.currentNode(),self.logic.inputFiducialNode, self.vsc.vtVars['cochleaSide'] )

      slicer.app.layoutManager().setLayout( slicer.modules.tables.logic().GetLayoutWithTable(slicer.app.layoutManager().layout))
      slicer.app.applicationLogic().GetSelectionNode().SetActiveTableID(self.logic.spTblNode.GetID())
      slicer.app.applicationLogic().PropagateTableSelection()

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
  #--------------------------------------------------------------------------------------------        
  #                       Segmentation Process 
  #--------------------------------------------------------------------------------------------        
  # This method perform the atlas segementation steps
  def run(self, inputVolumeNode, inputFiducialNode, cochleaSide):
      logging.info('Processing started')
 
      self.vsc   = VisSimCommon.VisSimCommonLogic()
      self.vsc.setGlobalVariables(0)

      self.inputVolumeNode = inputVolumeNode
      self.inputFiducialNode = inputFiducialNode
      # modality type CBCT, CT or MRI
      # it seems using CBCT atlas is enough 
      Styp="Dv"

      # segmentation atlas model paths
                     
      modelPath      =   self.vsc.vtVars['modelPath']+ ",Mdl"+Styp +cochleaSide +"c"+self.vsc.vtVars['imgType'] 
      modelPath      =   os.path.join(*modelPath.split(","))
      modelSegPath   =   self.vsc.vtVars['modelPath']+ ",Mdl"+Styp +cochleaSide +"cS.seg"+self.vsc.vtVars['imgType']        
      modelSegPath   =   os.path.join(*modelSegPath.split(","))
      modelImgStPath =   self.vsc.vtVars['modelPath']+ ",Mdl"+Styp +cochleaSide +"cSt.fcsv"       
      modelImgStPath   =   os.path.join(*modelImgStPath.split(","))       
      # set the results paths:       
      resTransPathOld  = os.path.join(self.vsc.vtVars['outputPath'] ,"TransformParameters.0.txt")
      resTransPath=resTransPathOld[0:-6]+'Pars.txt'        
      resOldDefPath = os.path.join(self.vsc.vtVars['outputPath'] , "deformationField"+self.vsc.vtVars['imgType'])
      resDefPath    = os.path.join(self.vsc.vtVars['outputPath'] , inputVolumeNode.GetName()+"_dFld"+self.vsc.vtVars['imgType'])
      inputImgName  = inputVolumeNode.GetStorageNode().GetFileName()
      inputImgName  = basename(os.path.splitext(inputImgName)[0])    
        
      segNodeName   = inputVolumeNode.GetName() + "_S.Seg"                 
      stpNodeName   = inputVolumeNode.GetName() + "_StPts"
      transNodeName = inputVolumeNode.GetName() + "_Transform"

      self.vsc.removeOtputsFolderContents()
      # check if the model is found
      if not isfile(modelPath): 
            print >> sys.stderr, "ERROR: model is not found"            
            print("modelPath: " + modelPath)
            return -1
      # endif

      # Get IJK point from the fiducial to use in cropping          
      inputPoint = self.vsc.ptRAS2IJK(inputFiducialNode,inputVolumeNode,0)
      # TODO: add better condition
      if  np.sum(inputPoint)== 0 :
            print("Error: select cochlea point")
            return -1
      #endif  
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , inputVolumeNode.GetName()+"_Cochlea_Pos.fcsv")                           
      sR = slicer.util.saveNode(inputFiducialNode, fnm )  

      #Remove old resulted nodes
      for node in slicer.util.getNodes():
          if ( segNodeName   == node): slicer.mrmlScene.RemoveNode(node) #endif
          if ( transNodeName == node): slicer.mrmlScene.RemoveNode(node) #endif
      #endfor    
        
      inputPointT = self.vsc.v2t(inputPoint) 
        
      print("=================== Cropping =====================")                           
      self.vsc.vtVars['intputCropPath'] = self.vsc.runCropping(inputVolumeNode, inputPointT,self.vsc.vtVars['croppingLength'],  self.vsc.vtVars['RSxyz'],  self.vsc.vtVars['hrChk'],0)                    
      [success, croppedNode] = slicer.util.loadVolume(self.vsc.vtVars['intputCropPath'], returnNode=True)
      croppedNode.SetName(inputVolumeNode.GetName()+"_Crop")                                                        
      print ("************  Register model to cropped input image **********************")
      cTI = self.vsc.runElastix(self.vsc.vtVars['elastixBinPath'],self.vsc.vtVars['intputCropPath'],  modelPath, self.vsc.vtVars['outputPath'], self.vsc.vtVars['parsPath'], self.vsc.vtVars['noOutput'], "292")
      copyfile(resTransPathOld, resTransPath)
      #genrates deformation field 
      cTR = self.vsc.runTransformix(self.vsc.vtVars['transformixBinPath'],modelPath, self.vsc.vtVars['outputPath'], resTransPath, self.vsc.vtVars['noOutput'], "295")
      # rename fthe file:
      os.rename(resOldDefPath,resDefPath)   
      print ("************  Load deformation field Transform  **********************")
      [success, vtTransformNode] = slicer.util.loadTransform(resDefPath, returnNode = True)
      vtTransformNode.SetName(transNodeName)   
      print ("************  Transform The Segmentation **********************")     
      [success, vtSegNode] = slicer.util.loadSegmentation(modelSegPath, returnNode = True)
      vtSegNode.SetName(segNodeName)
      vtSegNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) 
      #export seg to lbl then export back with input image as reference
      slicer.vtkSlicerTransformLogic().hardenTransform(vtSegNode)     # apply the transform
      vtSegNode.CreateClosedSurfaceRepresentation() 
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , vtSegNode.GetName()+".nrrd")                             
      sR = slicer.util.saveNode(vtSegNode, fnm )  

      print ("************  Transform The Scala Points **********************")
      #TODO:check the right side model
      #TODO:add scala vestibuli model
      # transform the Scala Tympani Points for length Computation 
      [success, vtImgStNode] = slicer.util.loadMarkupsFiducialList  (modelImgStPath, returnNode = True)
      vtImgStNode.GetDisplayNode().SetSelectedColor(1,0,0)  
      vtImgStNode.GetDisplayNode().SetTextScale(0.5)                               
      vtImgStNode.SetName(stpNodeName)
      vtImgStNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) 
      slicer.vtkSlicerTransformLogic().hardenTransform(vtImgStNode) # apply the transform
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , vtImgStNode.GetName()+".fcsv")                             
      sR = slicer.util.saveNode(vtImgStNode, fnm ) 

      # Display the result if no error
      # Clear cochlea location labels
      if  (cTI==0) and (cTR==0):
          # change the model type from vtk to stl 
          msn=slicer.vtkMRMLModelStorageNode()
          msn.SetDefaultWriteFileExtension('stl')
          slicer.mrmlScene.AddDefaultNode(msn)        
          print("get Cochlea information")
          tableName =  inputVolumeNode.GetName()+"_tbl"
          # create only if it does not exist
          try:
             spTblNode =  slicer.util.getNode(tableName)
          except Exception as e:
             print(e)   
             spTblNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
             spTblNode.SetName(tableName)
          #endtry
          spTblNode = self.vsc.getItemInfo( vtSegNode, croppedNode, spTblNode,0)
          for i in range (0,8):
              spTblNode.RemoveColumn(3)
          #endfor           
          spTblNode.GetTable().GetColumn(1).SetName("Size (mm^3)")
          spTblNode.GetTable().GetColumn(2).SetName("Length (mm)")   
          spTblNode.SetCellText(0,0,"Scala Tympani")
          spTblNode.SetCellText(1,0,"Scala Vestibuli")
          #spTblNode.resultsTableNode.SetCellText(0,2,"ST value")   
          #spTblNode.resultsTableNode.SetCellText(1,2,"0")   

          self.vsc.getFiducilsDistance(vtImgStNode,spTblNode )      
          spTblNode.RemoveRow(spTblNode.GetNumberOfRows())              
          self.spTblNode=spTblNode                     
      else:
           print("error happened during segmentation ")
      #endif

      #Remove temporary files and nodes:
      self.vsc.removeTmpsFiles()     
      print("================= Cochlea analysis is complete  =====================")
      logging.info('Processing completed')
      return vtSegNode

    #enddef
 
#===================================================================
#                           Test
#===================================================================
class CochleaSegTest(ScriptedLoadableModuleTest):
  def setUp(self):
      slicer.mrmlScene.Clear(0)
  #endef   

  def runTest(self):
      self.setUp()
      self.testSlicerCochleaSegmentation()
  #enddef

  def testSlicerCochleaSegmentation(self, imgPath=None, cochleaPoint=None, cochleaSide=None):

      self.delayDisplay("Starting testSlicerCochleaSegmentation test")
      self.stm=time.time()

      self.vsc   = VisSimCommon.VisSimCommonLogic()   
      self.vsc.vtVars = self.vsc.setGlobalVariables(0)
      self.logic = CochleaSegLogic()
      # remove contents of output folder
      self.vsc.removeOtputsFolderContents()
      #TODO: error handling to select the download link
      if cochleaSide is None:
         cochleaSide  = "L"  ;    beforORafter ="_b" # _a= before, _b=after     
         if ( cochleaSide=="L" and beforORafter=="_b" ):        
             cochleaPoint = [195,218,95]
             urisUniKo    = "https://cloud.uni-koblenz-landau.de/s/qMG2WPjTXabzcbX/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100001_DV_L_b.nrrd'
             uris = urisGitHub          
             fileNames    = 'P100001_DV_L_b.nrrd'
             nodeNames    = 'P100001_DV_L_b'
             checksums    = 'SHA256:9a5722679caa978b1a566f4a148c8759ce38158ca75813925a2d4f964fdeebf5'
         elif(cochleaSide=="L" and beforORafter=="_a"  ):
             cochleaPoint = [214,242,78]
             urisUniKo         = "https://cloud.uni-koblenz-landau.de/s/EwQiQidXqTcGySB/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100001_DV_L_a.nrrd'
             uris = urisGitHub          
             fileNames    = 'P100001_DV_L_a.nrrd'
             nodeNames    = 'P100001_DV_L_a'
             checksums    = 'SHA256:d7cda4e106294a59591f03e74fbe9ecffa322dd1a9010b4d0590b377acc05eb5'
         elif(cochleaSide=="R" and beforORafter=="_b" ):
             cochleaPoint = [194,216,93]
             urisUniKo   = "https://cloud.uni-koblenz-landau.de/s/4K5gAwisgqSHK4j/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100003_DV_R_b.nrrd'
             uris = urisGitHub          
             fileNames    = 'P100003_DV_R_b.nrrd' 
             nodeNames    = 'P100003_DV_R_b'
             checksums    = 'SHA256:4478778377982b6789ddf8f5ccd20f66757d6733853cce3f89faf75df2fa4faa'
         elif(cochleaSide=="R" and beforORafter=="_a" ):
             cochleaPoint = [294,250,60]
             urisUniKo    = "https://cloud.uni-koblenz-landau.de/s/WAxHyqLC3JsKY2x/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100003_DV_R_a.nrrd'
             uris = urisGitHub          
             fileNames    = 'P100003_DV_R_a.nrrd'
             nodeNames    = 'P100003_DV_R_a'
             checksums    = 'SHA256:c62d37e13596eafc8550f488006995d811c8d6503445d5324810248a3c3b6f89'
         else:
             print("error in cochlea side or before after type")
             return -1
      #endif  
      #sampledata loads the volume as well but didn't provide storage node.
      if imgPath is None: 
         try:
            tmpVolumeNode =  SampleData.downloadFromURL(uris, fileNames, nodeNames, checksums )[0]
            imgPath       =  os.path.join(slicer.mrmlScene.GetCacheManager().GetRemoteCacheDirectory(),fileNames)
            slicer.mrmlScene.RemoveNode(tmpVolumeNode)
         except Exception as e:
            print("Error: can not download sample data")
            print (e)
            return -1 
         #endtry
      else:
         nodeNames = os.path.splitext(os.path.basename(imgPath))[0]
      #endif 
      [success, inputVolumeNode]  = slicer.util.loadVolume(imgPath, returnNode=True)
      inputVolumeNode.SetName(nodeNames)

      # create a fiducial node for cochlea location for cropping    
      cochleaPointRAS = self.vsc.ptIJK2RAS(cochleaPoint,inputVolumeNode) 
      inputFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      inputFiducialNode.CreateDefaultDisplayNodes()
      inputFiducialNode.SetName("cochleaLocationPoint")  
      inputFiducialNode.AddFiducialFromArray(cochleaPointRAS)

      # run the segmentation
      segNode = self.logic.run(inputVolumeNode, inputFiducialNode, cochleaSide)
      #display:
      try:
         self.vsc.dispSeg(inputVolumeNode,segNode,34) # 34: 4up table layout
      except Exception as e:
             print("Can not display results! probably an external call ...")
             print(e)   
      #endtry  
      self.etm=time.time()
      tm=self.etm - self.stm
      print("Time: "+str(tm)+"  seconds")
      self.delayDisplay('Test testSlicerCochleaSegmentation passed!')
  #enddef
#endclass
