#======================================================================================
#  3D Slicer [1] plugin that uses elastix toolbox [2] Plugin for Automatic Cochlea    #
#  Image Segmentation [3,4]. More info can be found at [5].                           #
#  Sample cochlea datasets can be downloaded using Slicer Datastore module            #
#                                                                                     #
#  Contributers:                                                                      #
#      - Ibraheem Al-Dhamari, ia@idhamari.com                                         #
#  [1] https://www.slicer.org                                                         #
#  [2] http://elastix.isi.uu.nl                                                       #
#  [3] Al-Dhamari et al.,(2018), Automatic Cochlear Length and Volume Size Estimation #
#       First  International Workshop on Context-Aware Operating                      #
#       Theater OR 2, MICCAI 2018, Granada Spain.                                     #
#  [4] Al-Dhamari et al.,(2023), Automatic Cochlea Multimodal 3D Image Segmentation   #
#      And Analysis Using Atlas-model-based Method. Cochlea Implant International     #
#      https://doi.org/10.1080/14670100.2023.2274199                                  #
#  [5] https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea                 #
#                                                                                     #
#-------------------------------------------------------------------------------------#
#  Slicer 5.4.0                                                                       #
#  Updated: 18.11.2023                                                                #
#======================================================================================
# Non Slicer libs
from __future__ import print_function
import os, sys, time, re, shutil,  math, unittest, logging, zipfile, platform, subprocess, hashlib
from shutil import copyfile

from six.moves.urllib.request import urlretrieve
import numpy as np
import SimpleITK as sitk

# Slicer related
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import sitkUtils
import SampleData
import SegmentStatistics
import Elastix
import VisSimCommon

# TODOS:
# Update the models 
 
# Terminology
#  img         : ITK image
#  imgNode     : Slicer Node
#  imgName     : Filename without the path and without extension
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
#===================================================================
#                           Main Widget
#===================================================================
class CochleaSegWidget(ScriptedLoadableModuleWidget):


  def setup(self):
    print(" ")
    print("=======================================================")
    print("   Automatic Cochlea Image Segmentation and Analysis   ")
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
    self.inputSelectorCoBx.addEnabled         = False
    self.inputSelectorCoBx.removeEnabled      = False
    self.inputSelectorCoBx.noneEnabled        = False
    self.inputSelectorCoBx.showHidden         = False
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

    self.layout.addStretch(1) # Collapsible button is held in place when collapsing/expanding.
  #------------------------------------------------------------------------
  #                        Define GUI Elements Functions
  #------------------------------------------------------------------------
  # decide which cochlea side to segment
  # TODO: automate the process
  def onSideChkBoxChange(self):
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      self.vsc.setItemChk("cochleaSide", self.sideChkBox.checked, "cochleaSide", nodes)
  
  # Updating the cochlea length if the
  #  the user correct some fiducial points
  def onUpdateLengthBtnClick(self):
      # get table node:
      spTblNode = None
      nodes = slicer.util.getNodesByClass('vtkMRMLTableNode')
      for f in nodes:
          if ( "_tbl" in f.GetName() ):
             spTblNode = f
             break
      # fidNames=["_StPts","_SvPts","_StLtPts","_StOcPts","_avPts"]
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
           if ( "_StPts" in f.GetName()) :
                  # update length
                  stLength  = str(self.vsc.getFiducilsDistance(f))
                  spTblNode.SetCellText(0,2, stLength)
           elif ( "_SvPts" in f.GetName()) :
                  # update length
                  svLength  = str(self.vsc.getFiducilsDistance(f))
                  spTblNode.SetCellText(1,2, svLength)
           elif ( "_StLtPts" in f.GetName()) :
                  # update length
                  stLtLength  = str(self.vsc.getFiducilsDistance(f))
                  spTblNode.SetCellText(2,2, stLtLength)
           elif ( "_StOcPts" in f.GetName()) :
                  # update length
                  stOcLength  = str(self.vsc.getFiducilsDistance(f))
                  spTblNode.SetCellText(3,2, stOcLength)
           elif ( "_avPts" in f.GetName()) :
                  # update length
                  aVal  = self.vsc.getFiducilsDistance(f)
                  aValLengths  =self.logic.getAvalueLengths(aVal)
                  spTblNode.SetCellText(4,2, str(aValLengths[0]))
                  spTblNode.SetCellText(5,2, str(aValLengths[1]))
                  spTblNode.SetCellText(6,2, str(aValLengths[2]))
           spTblNode.Modified()

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

      if not hasattr(self.vsc, 'vtVars'):
         self.vsc.setGlobalVariables(0)

      self.vsc.locateItem(self.inputSelectorCoBx.currentNode(), self.inputPointEdt, 0 , 0)
      self.logic.inputFiducialNode= self.vsc.inputFiducialNodes[0]

      self.inputFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")



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
  
#===================================================================
#                           Logic
#===================================================================
class CochleaSegLogic(ScriptedLoadableModuleLogic):
  #--------------------------------------------------------------------------------------------
  #                       Segmentation Process
  #--------------------------------------------------------------------------------------------
  # This method perform the atlas segementation steps
  def run(self, inputVolumeNode, inputFiducialNode, cochleaSide, customisedOutputPath=None,customisedParPath=None):
    logging.info('Processing started')
 
    self.vsc   = VisSimCommon.VisSimCommonLogic()
    self.vsc.setGlobalVariables(0)
    
    if customisedOutputPath is not None: 
       self.vsc.vtVars['outputPath'] = customisedOutputPath
    if customisedParPath is not None: 
       self.vsc.vtVars['parsPath'] = customisedParPath
    
    print("inputVolumeNode       = ",inputVolumeNode.GetName())
    print("inputFiducialNode     = ", inputFiducialNode.GetName())
    print("outputPath            = ", self.vsc.vtVars['outputPath'])
    print("parsPath              = ", self.vsc.vtVars['parsPath'])
    
    self.inputVolumeNode = inputVolumeNode
    self.inputFiducialNode = inputFiducialNode
    
    # modality type CBCT, CT or MRI
    # it seems using CBCT atlas is enough
    Styp="Dv"
 
    # segmentation atlas model paths
    modelPath      =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"c"+self.vsc.vtVars['imgType'])
    modelSegPath   =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"cS.seg"+self.vsc.vtVars['imgType'])
    #modelSegPath   =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"cHS.seg"+self.vsc.vtVars['imgType']) #High resolution
 
    # points model path
    modelImgStPtPath   =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"c_StPt.fcsv") # Scala Tympani
    modelImgSvPtPath   =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"c_SvPt.fcsv") # Scala Vestibuli
    modelImgStLtPtPath =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"c_StLtPt.fcsv") # Scala Tympani A-value Lateral
    modelImgStOcPtPath =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"c_StOcPt.fcsv") # Scala Tympani A-value Organ of corti
    modelImgAvPtPath   =   os.path.join(self.vsc.vtVars['modelPath'] , "Mdl"+Styp +cochleaSide +"c_AvPt.fcsv") # A-value two points 
 
    # set the results paths:
    resImgPathOld  = os.path.join(self.vsc.vtVars['outputPath'] ,"result.0.nrrd")
    resImgRgPath   = os.path.join(self.vsc.vtVars['outputPath'] ,"result.Rg.nrrd") 
    resImgNRgPath   = os.path.join(self.vsc.vtVars['outputPath'] ,"result.NRg.nrrd") 
    
    resTransPathOld  = os.path.join(self.vsc.vtVars['outputPath'] ,"TransformParameters.0.txt")
    resTransRgPath=resTransPathOld[0:-6]+'_Rg_Pars.txt'
    resTransNRgPath=resTransPathOld[0:-6]+'_NRg_Pars.txt'

    node_name = inputVolumeNode.GetName()

    resOldDefPath = os.path.join(self.vsc.vtVars['outputPath'] , "deformationField"+self.vsc.vtVars['imgType'])
    resDefRgPath    = os.path.join(self.vsc.vtVars['outputPath'] , node_name+"_Rg_dFld"+self.vsc.vtVars['imgType'])
    resDefNRgPath    = os.path.join(self.vsc.vtVars['outputPath'] , node_name+"_NRg_dFld"+self.vsc.vtVars['imgType'])

    segNodeName   = node_name + "_S.Seg"

    stPtNodeName   = node_name + "_StPts"
    svPtNodeName   = node_name + "_SvPts"
    stLtPtNodeName = node_name + "_StLtPts"
    stOcPtNodeName = node_name + "_StOcPts"
    avPtNodeName   = node_name + "_avPts"

    transRgNodeName = node_name  + "_Rg_Transform"
    transNRgNodeName = node_name + "_NRg_Transform"

    self.vsc.removeOtputsFolderContents()
      # check if the model is found
    if not os.path.isfile(modelPath):
        print("ERROR: model is not found", file=sys.stderr)
        print("modelPath: " + modelPath)
        return -1
    # endif
 
    # Get IJK point from the fiducial to use in cropping
    inputPoint = self.vsc.ptRAS2IJK(inputFiducialNode,inputVolumeNode,0)
    # TODO: add better condition
    if  np.sum(inputPoint)== 0 :
           print("Error: select cochlea point")
           return -1

    fnm = os.path.join(self.vsc.vtVars['outputPath'] , inputVolumeNode.GetName()+"_Cochlea_Pos.fcsv")
    sR = slicer.util.saveNode(inputFiducialNode, fnm )
    
    #Remove old resulted nodes
    for node in slicer.util.getNodes():
         if ( segNodeName   == node): slicer.mrmlScene.RemoveNode(node)  
         if ( transRgNodeName == node): slicer.mrmlScene.RemoveNode(node)  
         if ( transNRgNodeName == node): slicer.mrmlScene.RemoveNode(node)  
 
    inputPointT = self.vsc.v2t(inputPoint)
    
    print("=================== Cropping =====================")
    self.vsc.vtVars['intputCropPath'] = self.vsc.runCropping(inputVolumeNode, inputPointT,self.vsc.vtVars['croppingLength'],  self.vsc.vtVars['RSxyz'],  self.vsc.vtVars['hrChk'],0)
    croppedNode = slicer.util.loadVolume(self.vsc.vtVars['intputCropPath'])
    croppedNode.SetName(inputVolumeNode.GetName()+"_Crop")
    
    print("=================== Registration =====================")
    
    print ("************  Rigid Registeration: model to cropped input image **********************")
     
    cTIr = self.vsc.runElastix(self.vsc.vtVars['elastixBinPath'],self.vsc.vtVars['intputCropPath'],  modelPath, self.vsc.vtVars['outputPath'], self.vsc.vtVars['parsPath'], self.vsc.vtVars['noOutput'], "292")
    
    os.rename(resImgPathOld,resImgRgPath)
    os.rename(resTransPathOld,resTransRgPath)
    
    #genrates deformation field
    cTRr = self.vsc.runTransformix(self.vsc.vtVars['transformixBinPath'],modelPath, self.vsc.vtVars['outputPath'], resTransRgPath, self.vsc.vtVars['noOutput'], "295")
    # rename fthe file:
    os.rename(resOldDefPath,resDefRgPath)
     
    print ("************  Non-Rigid Registeration: registered model to cropped input image **********************")
     
    cTInr = self.vsc.runElastix(self.vsc.vtVars['elastixBinPath'],self.vsc.vtVars['intputCropPath'],  resImgRgPath, self.vsc.vtVars['outputPath'], self.vsc.vtVars['parsNRPath'], self.vsc.vtVars['noOutput'], "292")
     
    os.rename(resImgPathOld,resImgNRgPath)
    os.rename(resTransPathOld,resTransNRgPath)
    
    #genrates deformation field
    cTRnr = self.vsc.runTransformix(self.vsc.vtVars['transformixBinPath'],resImgNRgPath, self.vsc.vtVars['outputPath'], resTransNRgPath, self.vsc.vtVars['noOutput'], "295")
    # rename fthe file:
    os.rename(resOldDefPath,resDefNRgPath)
         
    print ("************  Load deformation field Transforms  **********************")
    vtRgTransformNode = slicer.util.loadTransform(resDefRgPath)
    vtRgTransformNode.SetName(transRgNodeName)
  
    vtNRgTransformNode = slicer.util.loadTransform(resDefNRgPath)
    vtNRgTransformNode.SetName(transNRgNodeName)
  
    #combine the transforms      
    vtNRgTransformNode.SetAndObserveTransformNodeID(vtRgTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(vtNRgTransformNode)    
    chTransformNode = vtNRgTransformNode
       
    print ("************  Transform The Segmentation **********************")
    chSegNode = slicer.util.loadSegmentation(modelSegPath)
    chSegNode.SetName(segNodeName)
      
    chSegNode.SetAndObserveTransformNodeID(chTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(chSegNode)     # apply the transform
    #export seg to lbl then export back with input image as reference
    chSegNode.CreateClosedSurfaceRepresentation()
    fnm = os.path.join(self.vsc.vtVars['outputPath'] , chSegNode.GetName()+".nrrd")
    sR = slicer.util.saveNode(chSegNode, fnm )

    print ("************  Transform The Scala Tympani (St) Points **********************")
    # transform the Scala Tympani Points for length Computation
    chImgStPtNode = slicer.util.loadMarkupsFiducialList  (modelImgStPtPath, returnNode = True)
    chImgStPtNode.GetDisplayNode().SetSelectedColor(1,1,0)
    chImgStPtNode.GetDisplayNode().SetTextScale(0)
    chImgStPtNode.GetDisplayNode().SetGlyphScale(0.1)
    chImgStPtNode.SetName(stPtNodeName)
# 
# #     chTransformNode.Inverse()
    chImgStPtNode.SetAndObserveTransformNodeID(chTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(chImgStPtNode)     # apply the transform
    fnm = os.path.join(self.vsc.vtVars['outputPath'] , chImgStPtNode.GetName()+".fcsv")
    sP1 = slicer.util.saveNode(chImgStPtNode, fnm )
    self.vsc.vtVars['StLength']  = str(self.vsc.getFiducilsDistance(chImgStPtNode ))
 
    print ("************  Transform The Scala Vestibuli (Sv) Points **********************")
    chImgSvPtNode = slicer.util.loadMarkupsFiducialList  (modelImgSvPtPath, returnNode = True)
    chImgSvPtNode.GetDisplayNode().SetSelectedColor(0,0,0)
    chImgSvPtNode.GetDisplayNode().SetTextScale(0)
    chImgSvPtNode.GetDisplayNode().SetGlyphScale(0.1)
    chImgSvPtNode.SetName(svPtNodeName)
 
    chImgSvPtNode.SetAndObserveTransformNodeID(chTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(chImgSvPtNode) # apply the transform
 
    fnm = os.path.join(self.vsc.vtVars['outputPath'] , chImgSvPtNode.GetName()+".fcsv")
    sP2 = slicer.util.saveNode(chImgSvPtNode, fnm )
    self.vsc.vtVars['SvLength']  = str(self.vsc.getFiducilsDistance(chImgSvPtNode ))
 
    print ("************  Transform The A-value Scala Tempany Lateral (StL) Points *******")
    chImgStLtPtNode = slicer.util.loadMarkupsFiducialList  (modelImgStLtPtPath, returnNode = True)
    chImgStLtPtNode.GetDisplayNode().SetSelectedColor(1,0,1)
    chImgStLtPtNode.GetDisplayNode().SetTextScale(0)
    chImgStLtPtNode.GetDisplayNode().SetGlyphScale(0.1)
    chImgStLtPtNode.SetName(stLtPtNodeName)
#       
    chImgStLtPtNode.SetAndObserveTransformNodeID(chTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(chImgStLtPtNode) # apply the transform
       
    fnm = os.path.join(self.vsc.vtVars['outputPath'] , chImgStLtPtNode.GetName()+".fcsv")
    sP3 = slicer.util.saveNode(chImgStLtPtNode, fnm )
    self.vsc.vtVars['StLtLength']  = str(self.vsc.getFiducilsDistance(chImgStLtPtNode ))

    print ("************  Transform The A-value Organ od Corti (OC) Points ***************")
    chImgStOcPtNode = slicer.util.loadMarkupsFiducialList  (modelImgStOcPtPath, returnNode = True)
    chImgStOcPtNode.GetDisplayNode().SetSelectedColor(0,1,1)
    chImgStOcPtNode.GetDisplayNode().SetTextScale(0)
    chImgStOcPtNode.GetDisplayNode().SetGlyphScale(0.1)
    chImgStOcPtNode.SetName(stOcPtNodeName)
       
    chImgStOcPtNode.SetAndObserveTransformNodeID(chTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(chImgStOcPtNode) # apply the transform
       
    fnm = os.path.join(self.vsc.vtVars['outputPath'] , chImgStOcPtNode.GetName()+".fcsv")
    sP4 = slicer.util.saveNode(chImgStOcPtNode, fnm )
    self.vsc.vtVars['StOcLength']  = str(self.vsc.getFiducilsDistance(chImgStOcPtNode ))

    print ("************  Transform The A-value (Av)Points **********************")
    chImgAvPtNode = slicer.util.loadMarkupsFiducialList  (modelImgAvPtPath, returnNode = True)
    chImgAvPtNode.GetDisplayNode().SetSelectedColor(1,0,0)
    chImgAvPtNode.GetDisplayNode().SetTextScale(0)
    chImgAvPtNode.GetDisplayNode().SetGlyphScale(0.1)
    chImgAvPtNode.SetName(avPtNodeName)
       
    chImgAvPtNode.SetAndObserveTransformNodeID(chTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(chImgAvPtNode) # apply the transform
 
    print ("************  Computing information using A-value **********************")     
    fnm = os.path.join(self.vsc.vtVars['outputPath'] , chImgAvPtNode.GetName()+".fcsv")
    sP5 = slicer.util.saveNode(chImgAvPtNode, fnm )
    aVal =  self.vsc.getFiducilsDistance(chImgAvPtNode)

    print("A-value = " , aVal)
    aValLengths = self.getAvalueLengths(aVal )
    print("A-value Distance"  , aValLengths[0])     
    print("A-value Length Lt" , aValLengths[1])
    print("A-value Length Oc" , aValLengths[2])
    self.vsc.vtVars['AvalueDistance']   = aValLengths[0]
    self.vsc.vtVars['AvalueStLtLength'] = aValLengths[1]
    self.vsc.vtVars['AvalueStOcLength'] = aValLengths[2]      
    print("A-value Length Lt" , self.vsc.vtVars['AvalueStLtLength'])
    print("A-value Length Oc" , self.vsc.vtVars['AvalueStOcLength'])
    print("A-value Length Lt" , type(self.vsc.vtVars['AvalueStLtLength']))
    print("A-value Length Oc" , type(self.vsc.vtVars['AvalueStOcLength']))
 
    # Display the result if no error
    # Clear cochlea location labels
    if  (cTInr==0) and (cTRnr==0):
        # change the model type from vtk to stl
        msn=slicer.vtkMRMLModelStorageNode()
        msn.SetDefaultWriteFileExtension('stl')
        slicer.mrmlScene.AddDefaultNode(msn)
        print("get Cochlea information")
        tableName =  inputVolumeNode.GetName()+"_tbl"
        # create only if it does not exist
        try:
           spTblNode =  slicer.util.getNode(tableName)
           print("found ", tableName)    
        except Exception as e:
           print(e)
           print("creating  ", tableName)    
           spTblNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
           spTblNode.SetName(tableName)

        spTblNode = self.vsc.getItemInfo( chSegNode, croppedNode, spTblNode,0)
        for i in range (0,8):
            spTblNode.RemoveColumn(3)

        if spTblNode.GetNumberOfRows()>3:
           spTblNode.RemoveRow(spTblNode.GetNumberOfRows()-2)
           spTblNode.RemoveRow(spTblNode.GetNumberOfRows()-2)

        stVol = spTblNode.GetCellText(0,1)
        svVol = spTblNode.GetCellText(1,1)
        spTblNode.AddEmptyRow();          spTblNode.AddEmptyRow()     
        spTblNode.AddEmptyRow();          spTblNode.AddEmptyRow()
        spTblNode.AddEmptyRow()     
                                    
        spTblNode.GetTable().GetColumn(1).SetName("Size (mm^3)")
        spTblNode.GetTable().GetColumn(2).SetName("Length (mm)")
        spTblNode.SetCellText(0,0,"Scala Tympani")
        spTblNode.SetCellText(1,0,"Scala Vestibuli")
        
        # Lengths
        spTblNode.SetCellText(2,0,"Length StLt")
        spTblNode.SetCellText(3,0,"Length StOc")
        spTblNode.SetCellText(4,0,"A-value Distance")
        spTblNode.SetCellText(5,0,"A-value StLt")
        spTblNode.SetCellText(6,0,"A-value StOc")
        
        # volume 
        spTblNode.SetCellText(0,1,stVol)
        spTblNode.SetCellText(1,1,svVol)
        spTblNode.SetCellText(2,1,"")
        spTblNode.SetCellText(3,1,"")
 
        #TODO: change column width
        spTblNode.SetColumnProperty('A','Width','250')
           
        spTblNode.SetCellText(0,2,self.vsc.vtVars['StLength'])
        spTblNode.SetCellText(1,2,self.vsc.vtVars['SvLength'])
        spTblNode.SetCellText(2,2,self.vsc.vtVars['StLtLength'])          
        spTblNode.SetCellText(3,2,self.vsc.vtVars['StOcLength'])
        spTblNode.SetCellText(4,2, str( self.vsc.vtVars['AvalueDistance'] ) )
        spTblNode.SetCellText(5,2, str( self.vsc.vtVars['AvalueStLtLength'] ) )
        spTblNode.SetCellText(6,2, str( self.vsc.vtVars['AvalueStOcLength'] ) )          
        
        #spTblNode.RemoveRow(spTblNode.GetNumberOfRows())            
        self.spTblNode=spTblNode
        fnm = os.path.join(self.vsc.vtVars['outputPath'] , spTblNode.GetName()+".tsv")
        sR = slicer.util.saveNode(spTblNode, fnm )
    else:
         print("error happened during segmentation ")
 
    #Remove temporary files and nodes:
    self.vsc.removeTmpsFiles()
    print("================= Cochlea analysis is complete  =====================")
    logging.info('Processing completed')
    return chSegNode
      
 
  def getAvalueLengths(self,Aval):
      #  L= 8.58; cl1=L*3.86+4.99; cl2=L*4.16-5.05; print("CL1 = :", cl1, "      CL2 = :", cl2); 
      l1 =  Aval * 3.86 + 4.99 # lateral length
      l2 =  Aval * 4.16 - 5.05 # organ of corti length
      return Aval, l1, l2
 
#===================================================================
#                           Test
#===================================================================
class CochleaSegTest(ScriptedLoadableModuleTest):
  def setUp(self):
      slicer.mrmlScene.Clear(0)
 
  def runTest(self):
      self.setUp()
      #TODO: error handling to select the download link
      cochleaSide  = "L"  ;    beforORafter ="_a" # _a= before, _b=after
      fileName    = 'P100001_DV'
      nodeName    = 'P100001_DV'

      if ( cochleaSide=="L" and beforORafter=="_b" ):
             cochleaPoint = [195,218,93]
             urisUniKo    = "https://cloud.uni-koblenz-landau.de/s/qMG2WPjTXabzcbX/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100001_DV_L_b.nrrd'
             uris = urisGitHub
             fileName    = fileName + '_L_b.nrrd'
             nodeName    = nodeName + '_L_b'
             checksums    = 'SHA256:9a5722679caa978b1a566f4a148c8759ce38158ca75813925a2d4f964fdeebf5'
      elif(cochleaSide=="L" and beforORafter=="_a"  ):
             cochleaPoint = [214,242,78]
             urisUniKo         = "https://cloud.uni-koblenz-landau.de/s/EwQiQidXqTcGySB/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100001_DV_L_a.nrrd'
             uris = urisGitHub
             fileName    = fileName + '_L_a.nrrd'
             nodeName    = nodeName + '_L_a'
             checksums    = 'SHA256:d7cda4e106294a59591f03e74fbe9ecffa322dd1a9010b4d0590b377acc05eb5'
      elif(cochleaSide=="R" and beforORafter=="_b" ):
             cochleaPoint = [194,216,93]
             urisUniKo   = "https://cloud.uni-koblenz-landau.de/s/4K5gAwisgqSHK4j/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100003_DV_R_b.nrrd'
             uris = urisGitHub
             fileName    = fileName + '_R_b.nrrd'
             nodeName    = nodeName + '_R_b'
             checksums    = 'SHA256:4478778377982b6789ddf8f5ccd20f66757d6733853cce3f89faf75df2fa4faa'
      elif(cochleaSide=="R" and beforORafter=="_a" ):
             cochleaPoint = [294,250,60]
             urisUniKo    = "https://cloud.uni-koblenz-landau.de/s/WAxHyqLC3JsKY2x/download"
             urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100003_DV_R_a.nrrd'
             uris = urisGitHub
             fileName    = fileName + '_R_a.nrrd'
             nodeName    = nodeName + '_R_a'
             checksums    = 'SHA256:c62d37e13596eafc8550f488006995d811c8d6503445d5324810248a3c3b6f89'
      else:
             print("error in cochlea side or before after type")
             return -1
      #sampledata loads the volume as well but didn't provide storage node.

      try:
            #check if a local fils is available:
            imgPath       =  os.path.join(slicer.mrmlScene.GetCacheManager().GetRemoteCacheDirectory(),fileName)

            print(imgPath)
            if not os.path.isfile(imgPath):
               tmpVolumeNode =  SampleData.downloadFromURL(uris, fileName, nodeName, checksums )[0]
               slicer.mrmlScene.RemoveNode(tmpVolumeNode)
            else:
               print("Sample file exist, no need to download")
      except Exception as e:
            print("Error: can not download sample data")
            print (e)
            return -1

      self.testSlicerCochleaSegmentation(imgPath,cochleaPoint,cochleaSide)


  def testSlicerCochleaSegmentation(self, imgPath, cochleaPoint, cochleaSide, customisedOutputPath=None,customisedParPath=None ):

      self.delayDisplay("Starting testSlicerCochleaSegmentation test")
      self.stm=time.time()

      self.vsc   = VisSimCommon.VisSimCommonLogic()
      self.vsc.vtVars = self.vsc.setGlobalVariables(0)
      self.logic = CochleaSegLogic()
      # remove contents of output folder
      self.vsc.removeOtputsFolderContents()

      inputVolumeNode  = slicer.util.loadVolume(imgPath)
      print(imgPath)

      inputVolumeNode.SetName(os.path.splitext(os.path.basename(imgPath))[0])

      # create a fiducial node for cochlea location for cropping
      cochleaPointRAS = self.vsc.ptIJK2RAS(cochleaPoint,inputVolumeNode)
      inputFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      inputFiducialNode.CreateDefaultDisplayNodes()
      inputFiducialNode.SetName("cochleaLocationPoint")
      inputFiducialNode.AddFiducialFromArray(cochleaPointRAS)

      # run the segmentation
      segNode = self.logic.run(inputVolumeNode, inputFiducialNode, cochleaSide,customisedOutputPath,customisedParPath)
      #display:
      try:
         self.vsc.dispSeg(inputVolumeNode,segNode,34) # 34: 4up table layout
      except Exception as e:
             print("Can not display results! probably an external call ...")
             print(e)
      self.etm=time.time()
      tm=self.etm - self.stm
      print("Time: "+str(tm)+"  seconds")
      self.delayDisplay('Test testSlicerCochleaSegmentation passed!')
