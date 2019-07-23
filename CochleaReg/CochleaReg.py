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
#-------------------------------------------------------------------------------------#
#  Slicer 4.10.0                                                                      #
#  Updated: 24.6.2019                                                                 # 
#======================================================================================

import os, time, logging, unittest
import numpy as np
from __main__ import qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import SampleData

import VisSimCommon

# TODO:
# - Visualizing the interimediate steps.


# Terminology
#  img         : ITK image
#  imgNode     : Slicer Node
#  imgName     :  Filename without the path and without extension
#  imgPath     : wholePath + Filename and extension


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
  #end def init
#end class CochleaReg

#===================================================================
#                           Main Widget
#===================================================================
class CochleaRegWidget(ScriptedLoadableModuleWidget):
  def setup(self):
    print(" ")
    print("=======================================================")
    print("   Automatic Cochlea Image Registration                ")
    print("=======================================================")

    ScriptedLoadableModuleWidget.setup(self)

    # to access logic class functions and setup global variables
    self.logic = CochleaRegLogic()

    # Set default VisSIm location in the user home
    #TODO: add option user-defined path when installed first time
    self.vsc   = VisSimCommon.VisSimCommonLogic()
    self.vsc.setGlobalVariables(0)

    #-----------------------------------------------------------------
    #                     Create the GUI interface
    #-----------------------------------------------------------------
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
    self.fixedPointEdt.setText("[0,0,0]")

    # Create a textbox for cochlea location
    # TODO activate input IJK values as well
    self.movingPointEdt = qt.QLineEdit()
    self.movingPointEdt.setFixedHeight(40)
    self.movingPointEdt.setText("[0,0,0]")

    # Create a cochlea locator button
    self.fixedFiducialBtn = qt.QPushButton("Pick cochlea location in fixed image    ")
    self.fixedFiducialBtn.setFixedHeight(40)
    self.fixedFiducialBtn.setToolTip("Pick the fixed fiducial point that will be the center of the cropped image")
    self.fixedFiducialBtn.connect('clicked(bool)', lambda: self.onInputFiducialBtnClick("F"))
    self.mainFormLayout.addRow( self.fixedFiducialBtn, self.fixedPointEdt)

    # Create a cochlea locator button
    self.movingFiducialBtn = qt.QPushButton("Pick cochlea location in moving image    ")
    self.movingFiducialBtn.setFixedHeight(40)
    self.movingFiducialBtn.setToolTip("Pick the moving fiducial point that will be the center of the cropped image")
    self.movingFiducialBtn.connect('clicked(bool)', lambda: self.onInputFiducialBtnClick("M"))
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

  #------------------------------------------------------------------------
  #                        Define GUI Elements Functions
  #------------------------------------------------------------------------
  def onInputFiducialBtnClick(self, volumeType):
      if not hasattr(self.vsc, 'vtVars'):
         self.vsc.setGlobalVariables(0)
       #end
      self.fixedVolumeNode=self.fixedSelectorCoBx.currentNode()
      self.movingVolumeNode=self.movingSelectorCoBx.currentNode()
      self.logic.fixedFiducialNode = None
      self.logic.movingFiducialNode = None
      #remove old nodes
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if (f.GetName() == "_CochleaLocation") :
               slicer.mrmlScene.RemoveNode(f)
          #endif
      #endfor
      # Create Fiducial Node for the cochlea location in both images
      if (volumeType=="F"):
         print(" ..... getting cochlea location in the fixed image")
         self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: White  }")
         self.vsc.locateItem(self.fixedSelectorCoBx.currentNode(), self.fixedPointEdt,1, 0)
         self.fixedFiducialNode= self.vsc.inputFiducialNodes[1]
         self.fixedFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
      elif (volumeType=="M"):
         print(" ..... getting cochlea location in the fixed image")
         self.vsc.locateItem(self.movingSelectorCoBx.currentNode(), self.movingPointEdt,2, 0)
         self.movingFiducialNode= self.vsc.inputFiducialNodes[2]
         self.movingFiducialBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    #endif
  #enddef


  # An option to control results displaying
  def OnColorsChkBoxChange(self):
        print("color is changed")
        self.vsc.fuseWithOutColor(self.colorsChkBox.checked)

  def onApplyBtnClick(self):
      self.runBtn.setText("...please wait")
      self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
      slicer.app.processEvents()
      self.stm=time.time()
      print("time:" + str(self.stm))
      self.timeLbl.setText("                 Time: 00:00")

      print(type(self.fixedFiducialNode))
      # create an option to use IJK point or fidicual node
      registeredMovingVolumeNode =self.logic.run( self.fixedSelectorCoBx.currentNode(),self.fixedFiducialNode, self.movingSelectorCoBx.currentNode(),self.movingFiducialNode )
      self.vsc.fuseTwoImages(self.fixedSelectorCoBx.currentNode(), registeredMovingVolumeNode, True)
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
  #--------------------------------------------------------------------------------------------
  #                       Registration Process
  #--------------------------------------------------------------------------------------------
  # This method perform the registration steps
  def run(self, fixedVolumeNode, fixedFiducialNode, movingVolumeNode, movingFiducialNode):
      logging.info('Processing started')
      print(fixedVolumeNode.GetName())
      print(movingVolumeNode.GetName())
      self.vsc   = VisSimCommon.VisSimCommonLogic()
      self.vsc.setGlobalVariables(0)

      self.vsc.removeOtputsFolderContents()

      # results paths
      resTransPath  = os.path.join(self.vsc.vtVars['outputPath'] ,"TransformParameters.0.txt")
      resOldDefPath = os.path.join(self.vsc.vtVars['outputPath'] , "deformationField"+self.vsc.vtVars['imgType'])
      resDefPath    = os.path.join(self.vsc.vtVars['outputPath'] , movingVolumeNode.GetName()+"_dFld"+self.vsc.vtVars['imgType'])
      transNodeName = movingVolumeNode.GetName() + "_Transform"

      # Save original fixed and moving images
      if fixedVolumeNode.GetStorageNode() is None:
          fixedImgPath = os.path.join(self.vsc.vtVars['vissimPath'], fixedVolumeNode.GetName()+".nrrd")
          slicer.util.saveNode(fixedVolumeNode, fixedImgPath)
      fixedPath = fixedVolumeNode.GetStorageNode().GetFileName()

      if movingVolumeNode.GetStorageNode() is None:
          movingImgPath = os.path.join(self.vsc.vtVars['vissimPath'], movingVolumeNode.GetName()+".nrrd")
          slicer.util.saveNode(movingVolumeNode, movingImgPath)
      movingPath = movingVolumeNode.GetStorageNode().GetFileName()

      # Get IJK point from the fiducial to use in cropping
      fixedPoint = self.vsc.ptRAS2IJK(fixedFiducialNode,fixedVolumeNode,0)
      print("run fixed point: ============================")
      print(fixedVolumeNode.GetName())
      print(fixedPoint)
      # TODO: add better condition
      if  np.sum(fixedPoint)== 0 :
            print("Error: select cochlea fixed point")
            return -1
      #endif
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , fixedVolumeNode.GetName()+"_F_Cochlea_Pos.fcsv")
      sR = slicer.util.saveNode(fixedFiducialNode, fnm )

      movingPoint = self.vsc.ptRAS2IJK(movingFiducialNode,movingVolumeNode,0)
      print("run moving point: ============================")
      print(movingVolumeNode.GetName())
      print(movingPoint)

      # TODO: add better condition
      if  np.sum(fixedPoint)== 0 :
            print("Error: select cochlea moving point")
            return -1
      #endif
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , movingVolumeNode.GetName()+"_M_Cochlea_Pos.fcsv")
      sR = slicer.util.saveNode(movingFiducialNode, fnm )

      #Remove old resulted nodes
      #for node in slicer.util.getNodes():
      #    if ( "result"   in [node].GetName() ): slicer.mrmlScene.RemoveNode(node) #endif
      #endfor

      # TODO: add better condition
      if  (np.sum(fixedPoint)== 0) and (np.sum(movingPoint)== 0) :
            #qt.QMessageBox.critical(slicer.util.mainWindow(),'SlicerCochleaRegistration', 'Cochlea locations are missing')
            print("Error: select cochlea points in fixed and moving images")
            return False
      #endif

      fixedPointT = self.vsc.v2t(fixedPoint)
      movingPointT = self.vsc.v2t(movingPoint)

      print("=================== Cropping =====================")
      self.vsc.vtVars['fixedCropPath'] = self.vsc.runCropping(fixedVolumeNode, fixedPointT,self.vsc.vtVars['croppingLength'],  self.vsc.vtVars['RSxyz'],  self.vsc.vtVars['hrChk'],0)
      [success, croppedFixedNode] = slicer.util.loadVolume(self.vsc.vtVars['fixedCropPath'], returnNode=True)
      croppedFixedNode.SetName(fixedVolumeNode.GetName()+"_F_Crop")

      self.vsc.vtVars['movingCropPath'] = self.vsc.runCropping(movingVolumeNode, movingPointT,self.vsc.vtVars['croppingLength'],  self.vsc.vtVars['RSxyz'],  self.vsc.vtVars['hrChk'],0)
      [success, croppedMovingNode] = slicer.util.loadVolume(self.vsc.vtVars['movingCropPath'], returnNode=True)
      croppedMovingNode.SetName(movingVolumeNode.GetName()+"_M_Crop")
      print ("************  Register cropped moving image to cropped fixed image **********************")
      cTI = self.vsc.runElastix(self.vsc.vtVars['elastixBinPath'],self.vsc.vtVars['fixedCropPath'],  self.vsc.vtVars['movingCropPath'], self.vsc.vtVars['outputPath'], self.vsc.vtVars['parsPath'], self.vsc.vtVars['noOutput'], "336")
      #copyfile(resTransPathOld, resTransPath)
      #genrates deformation field
      cTR = self.vsc.runTransformix(self.vsc.vtVars['transformixBinPath'],self.vsc.vtVars['movingCropPath'], self.vsc.vtVars['outputPath'], resTransPath, self.vsc.vtVars['noOutput'], "339")
      # rename fthe file:
      os.rename(resOldDefPath,resDefPath)

      print ("************  Load deformation field Transform  **********************")
      [success, vtTransformNode] = slicer.util.loadTransform(resDefPath, returnNode = True)
      vtTransformNode.SetName(transNodeName)
      print ("************  Transform The Original Moving image **********************")
      movingVolumeNode.SetAndObserveTransformNodeID(vtTransformNode.GetID())
      #export seg to lbl then export back with input image as reference
      slicer.vtkSlicerTransformLogic().hardenTransform(movingVolumeNode)     # apply the transform
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , movingVolumeNode.GetName()+"_Registered.nrrd")
      sR = slicer.util.saveNode(movingVolumeNode, fnm )
      [success, registeredMovingVolumeNode] = slicer.util.loadVolume(fnm, returnNode = True)
      registeredMovingVolumeNode.SetName(movingVolumeNode.GetName()+"_Registered")
      #remove the tempnode and load the original
      slicer.mrmlScene.RemoveNode(movingVolumeNode)
      [success, movingVolumeNode] = slicer.util.loadVolume(movingPath, returnNode = True)
      movingVolumeNode.SetName(os.path.splitext(os.path.basename(movingVolumeNode.GetStorageNode().GetFileName()))[0])
      if  (cTI==0) and (cTR==0):
          print("No error is reported during registeration ...")
      else:
           print("error happened during registration ")
      #endif

      #Remove temporary files and nodes:
      self.vsc.removeTmpsFiles()
      print("================= Cochlea registration is complete  =====================")
      logging.info('Processing completed')

      return registeredMovingVolumeNode
    #enddef


  def runTest(self):
      self.setUp()
      self.testSlicerCochleaRegistration()
  #enddef
  def testSlicerCochleaRegistration(self, fixedImgPath=None, fixedPoint=None, movingImgPath=None, movingPoint=None):

      self.delayDisplay("Starting testSlicerCochleaRegistration test")
      self.stm=time.time()

      if fixedPoint is None:
          fixedPoint = [220,242,78]
      #endif
      if movingPoint is None:
          movingPoint = [196,217,93]
      #endif
      nodeNames='P100001_DV_L_a'
      fileNames='P100001_DV_L_a.nrrd'
      urisUniKo         = "https://cloud.uni-koblenz-landau.de/s/EwQiQidXqTcGySB/download"
      urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100001_DV_L_a.nrrd'
      uris = urisGitHub
      checksums='SHA256:d7cda4e106294a59591f03e74fbe9ecffa322dd1a9010b4d0590b377acc05eb5'
      if fixedImgPath is None:
         tmpVolumeNode =  SampleData.downloadFromURL(uris, fileNames, nodeNames, checksums )[0]
         fixedImgPath  =  os.path.join(slicer.mrmlScene.GetCacheManager().GetRemoteCacheDirectory(),fileNames)
         slicer.mrmlScene.RemoveNode(tmpVolumeNode)
      else:
         nodeNames = os.path.splitext(os.path.basename(fixedImgPath))[0]
      #endif
      [success, fixedVolumeNode]  = slicer.util.loadVolume(fixedImgPath, returnNode=True)
      fixedVolumeNode.SetName(nodeNames)
      #endifelse
      nodeNames='P100001_DV_L_b'
      fileNames='P100001_DV_L_b.nrrd'
      urisUniKo    = "https://cloud.uni-koblenz-landau.de/s/qMG2WPjTXabzcbX/download"
      urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/P100001_DV_L_b.nrrd'
      uris = urisGitHub
      checksums='SHA256:9a5722679caa978b1a566f4a148c8759ce38158ca75813925a2d4f964fdeebf5'
      if movingImgPath is None:
         tmpVolumeNode =  SampleData.downloadFromURL(uris, fileNames, nodeNames, checksums )[0]
         movingImgPath  =  os.path.join(slicer.mrmlScene.GetCacheManager().GetRemoteCacheDirectory(),fileNames)
         slicer.mrmlScene.RemoveNode(tmpVolumeNode)
      else:
         nodeNames = os.path.splitext(os.path.basename(movingImgPath))[0]
      #endif
      [success, movingVolumeNode] = slicer.util.loadVolume(movingImgPath, returnNode=True)
      movingVolumeNode.SetName(nodeNames)
      #endifelse
      self.logic = CochleaRegLogic()
      self.vsc   = VisSimCommon.VisSimCommonLogic()
      #setGlobal variables.
      self.vsc.vtVars = self.vsc.setGlobalVariables(0)

      # remove contents of output folder
      self.vsc.removeOtputsFolderContents()

      # record duration of the test

      # create a fiducial node for cochlea locations for cropping
      fixedPointRAS = self.vsc.ptIJK2RAS(fixedPoint , fixedVolumeNode)
      fixedFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      fixedFiducialNode.CreateDefaultDisplayNodes()
      fixedFiducialNode.SetName("F_cochleaLocationPoint")
      fixedFiducialNode.AddFiducialFromArray(fixedPointRAS)
      fixedFiducialNode.SetNthFiducialLabel(0, "F_CochleaLocation")

      movingPointRAS = self.vsc.ptIJK2RAS(movingPoint , movingVolumeNode)
      movingFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      movingFiducialNode.CreateDefaultDisplayNodes()
      movingFiducialNode.SetName("M_cochleaLocationPoint")
      movingFiducialNode.AddFiducialFromArray(movingPointRAS)
      movingFiducialNode.SetNthFiducialLabel(0, "M_CochleaLocation")

      # run the segmentation
      registeredMovingVolumeNode = self.logic.run(fixedVolumeNode, fixedFiducialNode, movingVolumeNode, movingFiducialNode)

      #display:
      try:
         self.vsc.fuseTwoImages(fixedVolumeNode, registeredMovingVolumeNode , True)
      except Exception as e:
             print("Can not display results! probably an external call ...")
             print(e)
      #endtry

      self.etm=time.time()
      tm=self.etm - self.stm
      print("Time: "+str(tm)+"  seconds")
      self.delayDisplay('Test testSlicerCochleaRegistration passed!')
  #enddef
#endclass
