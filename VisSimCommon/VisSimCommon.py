
#======================================================================================
#  3D Slicer [1] common functions used by VisSim extensions                           # 
#                                                                                     #
#  Contributers:                                                                      #   
#      - Ibraheem Al-Dhamari,  idhamari@uni-koblenz.de    : Plugin design.            #
#  [1] https://www.slicer.org                                                         #
#                                                                                     # 
#-------------------------------------------------------------------------------------#
#  Slicer 4.10
#  Updated: 19.6.2019                                                                 # #-------------------------------------------------------------------------------------#
# this file can be updated andreload automatically when call dependant module by      # 
# modifing bin/python/slicer/ScriptedLoadableModule.py                                #
#    def onReload(self):
#       .
#       . 
#       VisSimModules = ["CochleaSeg","CochleaReg","CervicalSpineTools","CervicalVertebraTools"]
#       if self.moduleName in VisSimModules:
#          print("Reloading VisSimCommon ............")
#          slicer.util.reloadScriptedModule("VisSimCommon")
#       slicer.util.reloadScriptedModule(self.moduleName)    
#======================================================================================
import os, re , datetime, time ,shutil, unittest, logging, zipfile,  stat,  inspect
import sitkUtils, sys ,math, platform  , glob,subprocess, urllib, urllib2
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
import SampleData

import SegmentStatistics
import Elastix

#TODO:
# - test on vsScritps
# - update cochlea GUI
# - Update private spine plugins
#  
# - remove temporary nodes and files in a cleaning procedure after segmentation
#       - can be activated or disables within the call 
#       - cropped images
# - fix cochlea right model segmentation
# - Update python 3 branches
 
#===================================================================
#                           Main Class
#===================================================================
class VisSimCommon(ScriptedLoadableModule):
  def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "VisSim Common"
        parent.categories = ["VisSimTools"]
        self.parent = parent        
  #enddef
#endclass

class VisSimCommonWidget(ScriptedLoadableModuleWidget):
  def setup(self):
      ScriptedLoadableModuleWidget.setup(self)
      self.mainCollapsibleBtn = ctk.ctkCollapsibleButton()
      self.mainCollapsibleBtn.text = "VisSim Common"
      self.layout.addWidget(self.mainCollapsibleBtn)

class VisSimCommonLogic(ScriptedLoadableModuleLogic):
  ElastixLogic = Elastix.ElastixLogic()
  ElastixBinFolder = ElastixLogic.getElastixBinDir()+"/"

  # these should be removed later
  vsID = "testing VisSimCommonLogic"
  def tstSum(self, x,y):
      print("testing")
      return x+y
  #enddef
 

  # vsExtension = 0: Cochlea, vsExtension = 1: Spine
  def setGlobalVariables(self,vsExtension):
      # define global variables as a dictonary
      self.vtVars = {}

      #shared stuff
      self.elastixEnv                     = self.ElastixLogic.getElastixEnv()  # to load elastix libs
      self.elastixStartupInfo             = self.ElastixLogic.getStartupInfo() # to hide the console
      self.vtVars['vissimPath']           = os.path.join(expanduser("~"),"VisSimTools")
      self.vtVars['elastixBinPath']       = os.path.join(self.ElastixBinFolder, "elastix")
      self.vtVars['transformixBinPath']   =  os.path.join(self.ElastixBinFolder, "transformix")  
      self.vtVars['winOS']                = "False"
      if (sys.platform == 'win32') or (platform.system()=='Windows'):
         self.vtVars['elastixBinPath']       = os.path.join(self.ElastixBinFolder, "elastix.exe")
         self.vtVars['transformixBinPath']   = os.path.join(self.ElastixBinFolder, "transformix.exe")
         #self.vtVars['noOutput']             = " >> /dev/null"
         self.vtVars['winOS']                     = "True"
      #endif   
      self.vtVars['noOutput']             = " >> /dev/null"
      self.vtVars['outputPath']           = os.path.join(self.vtVars['vissimPath'],"outputs")
      self.vtVars['imgType']              = ".nrrd"
      self.vtVars['hrChk']                = "True"
      self.vtVars['fixedPoint']           = "[0,0,0]" # initial poisition = no position               
      self.vtVars['movingPoint']          = "[0,0,0]" # initial poisition = no position               
      # change the model type from vtk to stl 
      msn=slicer.vtkMRMLModelStorageNode()
      msn.SetDefaultWriteFileExtension('stl')
      slicer.mrmlScene.AddDefaultNode(msn)

      if vsExtension == 0: #0=cochlea
         self.vtVars['othersUniKoWebLink']  =  ("https://cloud.uni-koblenz-landau.de/s/XYXPb4Fepms2JeC/download")  
         self.vtVars['othersWebLink']       =  ("https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/VisSimToolsCochlea.zip")  
         parsPath                            = self.vtVars['vissimPath']  + ",pars,parCochSeg.txt" 
         self.vtVars['parsPath']             = os.path.join(*parsPath.split(","))
         modelPath                           = self.vtVars['vissimPath']  + ",models,modelCochlea" 
         self.vtVars['modelPath']            = os.path.join(*modelPath.split(","))            
         self.vtVars['downSz']               = "500"
         self.vtVars['inputPoint']           = "[0,0,0]" # initial poisition = no position               
         self.vtVars['croppingLength']       = "[ 10 , 10 , 10 ]"   #Cropping Parameters
         self.vtVars['RSxyz']                = "[ 0.125, 0.125 , 0.125 ]"  #Resampling parameters 
         self.vtVars['dispViewTxt']          = "Green"
         self.vtVars['cochleaSide']          = "L" # default cochlea side is left
         self.vtVars['StLength']             = "0" # initial scala tympani length
         self.vtVars['dispViewTxt']          = "Green"
         self.vtVars['dispViewID']           = "8" #Green, coronal view 
      #Only for Cervical Spine      
      elif vsExtension == 1: # Cervical Spine       
         print("VisSimCommonLogic: initializing global variables:")  
         self.vtVars['othersUniKoWebLink']   = "https://cloud.uni-koblenz-landau.de/s/yfwcdymS9QfqKc9/download"
         self.vtVars['othersWebLink']        = "https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/VisSimToolsCervicalSpine.zip"
         parsPath                            = self.vtVars['vissimPath']  + ",pars,parSpiSeg.txt" 
         self.vtVars['parsPath']             = os.path.join(*parsPath.split(","))
         modelPath                           = self.vtVars['vissimPath']  + ",models,modelCervicalSpine" 
         self.vtVars['modelPath']            = os.path.join(*modelPath.split(","))
         self.vtVars['vtID']                 = "7"
         vtMethodsegT= [",Default"]
         vtMethodsgT = ["S.seg" ]
         self.vtVars['vtMethodID']           = "0"
         self.vtVars['segT']                 = vtMethodsegT[int(self.vtVars['vtMethodID'])] 
         self.vtVars['sgT']                  = vtMethodsgT[int(self.vtVars['vtMethodID'])] 
         self.vtVars['Styp']                 = "Ct"  
         self.vtVars['vtPtsLigDir']          = ",PtsLig"
         self.vtVars['vtPtsLigSuff']         = "Lp"
         modelCropImgLigPtsTmpPath           = self.vtVars['modelPath'] +"," +self.vtVars['vtPtsLigDir']+","+self.vtVars['Styp']
         self.vtVars['modelLigPtsPath']      = os.path.join(*modelCropImgLigPtsTmpPath.split(","))  
         subVarsTemplateFnm                  = self.vtVars['modelPath'] +","+self.vtVars['vtPtsLigDir']+",simPackSubVars.txt"  
         self.vtVars['subVarsTemplateFnm']   =  os.path.join(*subVarsTemplateFnm.split(","))  
         self.vtVars['dispViewTxt']          = "Yellow"
         self.vtVars['dispViewID']           = "7" #Yellow, sagittal view 
         self.vtVars['downSz']               = "160"
         self.vtVars['winOS']                = "False"
         self.vtVars['ligChk']               = "True"
         self.vtVars['segNodeCoM']           = "[ 0 , 0 , 0 ]"
         self.vtVars['croppingLength']       = "[ 80 , 80 , 50 ]"
         self.vtVars['RSxyz']                = "[ 0.5, 0.5 , 0.5 ]"
         #endif             
      #check if VisSimTools folder is found                  
      self.checkVisSimTools(self.vtVars,vsExtension)
      
      return self.vtVars
  #enddef


  def checkVisSimTools(self,vtVars,vsExtension ):

      print(" Defaults paths: " )
      print("      VisSimTools folder: " + vtVars['vissimPath'])
      print("      Output folder     : " + vtVars['outputPath'])            
      if isfile(vtVars['elastixBinPath'].strip()): 
          print("      elastix binaries are found in " + vtVars['elastixBinPath'] )
      else:           
          print("      elastix binaries are missing, please install SlicerElastix extension ... ")
          #TODO: download elastix binaries as additional option
          return -1
      #endif
      othersWebLink =  ""      
      if vsExtension ==0: # cochlea
         # TODO: optimise this part to download only the missing files        
         # Check if elastix exist or download it 
         print("      Cochlea Extension is selected")            
         # check if other files exist
         if  os.path.exists(vtVars['modelPath']): 
            print("      Other files are found !" )
            print("      Parameter file: "  + vtVars['parsPath'])
            print("      Cropping Length: " + vtVars['croppingLength'] )
         else:
            othersWebLink = vtVars['othersWebLink']      
                           
         #endif
      elif vsExtension ==1: # Spine
         # TODO: optimise this part to download only the missing files        
         # Check if elastix exist or download it 
         print("      Spine Extension is selected")            
         if  os.path.exists(vtVars['modelPath']): 
            print("      Other files are found !" )
            print("      Parameter file: "  + vtVars['parsPath'])
            print("      Cropping Length: " + vtVars['croppingLength'] )
         else:
            othersWebLink = vtVars['othersWebLink']      
         #endif
      #endif  
      if not othersWebLink=="":
         print("      Other files are  missing, trying to download ... ")
         try:                               
                print("      Downloading VisSimTools others ...")
                vissimZip = expanduser("~/VisSimToolsTmp.zip")      
                uFile = urllib.urlretrieve(othersWebLink,vissimZip)                       
                print ("     Extracting to user home ")
                zip_ref = zipfile.ZipFile(vissimZip, 'r')
                zip_ref.extractall(expanduser("~/"))
                zip_ref.close()  
                #remove the downloaded zip file     
                os.remove(vissimZip)   
                print ("    done! ")
         except Exception as e:
                print("      Error: can not download and extract VisSimTools ...")
                print(e)   
                return -1
       #end try-except  
      
  #enddef
    
  # string to boolean converter
  def s2b(self,s):
        return s.lower() in ("yes", "true", "t", "1")
  #enddef
  
        inputPointT = self.vsc.v2t(inputPoint) # "["+str(inputPoint[0])+","+str(inputPoint[1])+","+str(inputPoint[2])+"]"

  #convert a vector to text
  def v2t(self,v):
      return  "["+str(v[0])+","+str(v[1])+","+str(v[2])+"]"
  #enddef

  #convert dictonary text to vector
  def t2v(self,txt):
      vector = [0,0,0]
      print(txt)
      t = txt.strip("]").strip("[").strip("(").strip(")")
      t = t.split(",")
      for i in range(3):
          vector[i] =float(t[i])
      return vector 
  #enddef

#------------------------------------------------------
#                  IJK to RAS  
#------------------------------------------------------
# This function convert an IJK point to RAS point 
#  input:  a point vector and volume node
#  output: a point vector                     
  def ptIJK2RAS(self,ptIJK,inputImg): # imgPath or imgNode are supported
        #TODO: add option for printing                   
        # create a IJK2RAs transformation matrix 
        inputImgNode =inputImg 
        if (isinstance(inputImg, str)):
            [success, inputImgNode] = slicer.util.loadVolume( inputImg, returnNode=True)
        #endif
        ijk2rasM = vtk.vtkMatrix4x4()
        inputImgNode.GetIJKToRASMatrix(ijk2rasM)
        ptRAS=np.zeros((len(ptIJK),3))
        ijk= ptIJK
        # create a 4 elements array to get the converted values
        ijkv=[ijk[0],ijk[1],ijk[2],1]             
        rasPt=ijk2rasM.MultiplyDoublePoint(ijkv)
        ptRAS=rasPt[0:3]
        if (isinstance(inputImg, str)):
           slicer.mrmlScene.RemoveNode(inputImgNode )
        #endif

        return  ptRAS       
   


#------------------------------------------------------
#                 RAS  to IJK 
#------------------------------------------------------   
# This function convert RAS ro an IJK point 
#  input:  a point vector and volume node
#  output: a point vector                     
  def ptRAS2IJK(self,ptRAS,inputImg,i):      
      inputImgNode =inputImg 
      if (isinstance(inputImg, str)):
            [success, inputImgNode] = slicer.util.loadVolume( inputImg, returnNode=True)
      #endif

      # create a RAS2IJK transformation matrix 
      ras2ijkM = vtk.vtkMatrix4x4()
      inputImgNode.GetRASToIJKMatrix(ras2ijkM)       
      ras=[0,0,0]
      #'vtkSlicerMarkupsModuleMRMLPython.vtkMRMLMarkupsFiducialNode'
      if not type(ptRAS) =='str':
         #TODO: add option for printing       
         if not i is None:
            print(i)
            print(type(ptRAS)) 
            ptRAS.GetNthFiducialPosition(i,ras)
         else:
            ptRAS.GetNthFiducialPosition(0,ras)             
         #endif       
      else:
            ras=self.t2v(ptRAS) 
      #endif              
      # create a 4 elements array to get the converted values
      rasv=[ras[0],ras[1],ras[2],1]             
      ptIJKf=np.zeros(3);
      ijkPt=ras2ijkM.MultiplyPoint(rasv)
      ptIJKf[0]=ijkPt[0];ptIJKf [1]=ijkPt[1];ptIJKf [2]=ijkPt[2];
      ptIJK = ptIJKf.astype(np.int64)
      #print("RAS= " + str(ras)+ "   IJK= " + str(ptIJK))
      return  ptIJK       
    
    
  #------------------------------------------------------
  #         image2points
  #------------------------------------------------------
  # This function converted a binary image to fiducials points
  # it is used for computing the length of the scala tympani
  # as skelton methods don't work. 
  # Todo: redesign to be used externally
  #        we need to input a volume node only
  def image2points(self,inputImgNode):          
        print("====================================================")            
        print("=                Image to Points                  =")      
        print("====================================================") 
        # clone the input image       
        nimg = sitkUtils.PullVolumeFromSlicer(inputImgNode.GetID())
        sz= inputImgNode.GetImageData().GetDimensions()
        tmpImgArray = slicer.util.array(inputImgNode.GetID())
        nimgMax = tmpImgArray.max()
        nimgMin = tmpImgArray.min()
          
        b= list(zip( *np.where(tmpImgArray > (nimgMax/2 ) ) ))
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
        ptsRAS   =np.zeros((NoPts,3))        
        for j in range (0,NoPts):
           ptsIJK[j][0:3] = ptsIJKtmp[j][0:3] 
           ptsRAS[j]=  self.ptIJK2RAS(ptsIJK[j],inputImgNode)           
           x = ptsRAS[j][0]; y = ptsRAS[j][1] ; z = ptsRAS[j][2]
           self.markupsNode.AddFiducial(x,y,z)
        #endfor        
                
        # get the file name at the first of the plugin
        # this keeps the fiducials with short name
        #self.markupsNode.SetName(self.inputFnm + "StPts")
        self.markupsNode.SetName(inputImgNode.GetName()+ "_StPts")
        return self.markupsNode

#enddef


  #-----------------------------------------------------------------------------------
  #                       Cropping Process  
  #--------------------------------------------------------------------------------------------
  # Using the location as a center point, we cropp around it using the defined cropLength 
  # point must be a string in IJK format e.g. "[190,214,92]"
  # this is useful to call the function from console with some arguments 
  def runCropping(self, inputVolume, pointT,croppingLengthT, samplingLengthT, hrChkT,  vtIDt):

        print("================= Begin cropping  ... =====================")
        # Create a temporary node as workaround for bad path or filename 
        #TODO: create a temp folder and remove temp node before display
        print(" location: " + pointT + "   cropping length: " + str(croppingLengthT) )
        nodeName    = inputVolume.GetName() +"_Crop" 
        nodeNameIso = inputVolume.GetName() +"_CropIso"        
        inputCropPath = self.vtVars['vissimPath']+","+nodeName +".nrrd"                                 
        inputCropPath = os.path.join(*inputCropPath.split(","))     
        inputCropIsoPath = self.vtVars['vissimPath']+","+nodeNameIso +".nrrd"                               
        inputCropIsoPath = os.path.join(*inputCropIsoPath.split(","))     

        #for Spine
        if not vtIDt == 0:
           print("Vertebra "+vtIDt+ " location: " + pointT + "   cropping length: " + str(croppingLengthT) )
           nodeName    = inputVolume.GetName() +"_C" + vtIDt
           nodeNameIso = inputVolume.GetName() +"_C" + vtIDt +"_iso"        
           inputCropPath = self.vtVars['vissimPath']+","+nodeName  +".nrrd"                                        
           inputCropPath = os.path.join(*inputCropPath.split(","))     
           inputCropIsoPath = self.vtVars['vissimPath']+","+nodeNameIso  +".nrrd"
           inputCropIsoPath = os.path.join(*inputCropIsoPath.split(","))     
        #endif
        
        croppingLength =   self.t2v(croppingLengthT)
        samplingLength =   self.t2v(samplingLengthT)
        hrChk          =   self.s2b(hrChkT)
        point          =   self.t2v(pointT)
        print("location: " + str(point) + "   cropping length: " + str(croppingLength) )
        #Remove old cropping node
        nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in nodes:
            if ("_Crop" in f.GetName()):
                slicer.mrmlScene.RemoveNode(f )
            #endif
        #endfor
        # resampling spacing  
        self.RSx= samplingLength[0] ; self.RSy=samplingLength[1];     self.RSz= samplingLength[2]

        #Get input image information 
        spacing = inputVolume.GetSpacing()
        imgData = inputVolume.GetImageData()
        dimensions = imgData.GetDimensions()
        
        # compute cropping bounds from image information and cropping parameters
        croppingBounds = [[0,0,0],[0,0,0]];   size = [0,0,0];    lower = [0,0,0] ;     upper = [0,0,0]
        for i in range(0,3):
            size[i] = int((croppingLength[i]/spacing[i])/2)
            lower[i] = int(point[i]) - int(size[i])
            upper[i] = dimensions[i] - int(point[i]+size[i])
            # Check if calculated boundaries exceed image dimensions
            if lower[i] < 0:
                    lower[i] = 0
            #endif        
            if upper[i] > dimensions[i]:
                   upper[i] = int(dimensions[i])
            #endif
        #endfor   
        croppingBounds = [lower,upper]
        # Call SimpleITK CropImageFilter
        print("Cropping with " + str(croppingBounds[0]) + " and " + str(croppingBounds[1]) + ".")
        #self.inputCropPath = "bla blaa blaaa"

        inputImage = sitkUtils.PullVolumeFromSlicer(inputVolume.GetID())
        cropper = sitkUtils.sitk.CropImageFilter()
        croppedImage = cropper.Execute(inputImage, croppingBounds[0], croppingBounds[1])                  
        # Make a node with cropped image 
        sitkUtils.PushVolumeToSlicer(croppedImage, None, nodeName , 'vtkMRMLScalarVolumeNode' )
        crNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in crNodes:
            if nodeName in f.GetName():
                 f.SetName(nodeName) 
                 break         
            #endif
        #endfor  
        croppedNode = slicer.util.getNode(nodeName)
        print("cropped:     "+inputCropPath)
        slicer.util.saveNode( croppedNode, inputCropPath)

        #-------------------------------------------------------
        # Resampling: this produces better looking models  
        #-------------------------------------------------------
        #TODO: separate this in  a new function
        if hrChk:
           # remove the old cropped node
           nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
           for f in nodes:
               if ("_Crop" in f.GetName()):
                  slicer.mrmlScene.RemoveNode(f )
               #endif
           #endfor
           #Run slicer cli module: resample scalar volume
           #inputCropIsoPath = os.path.splitext(inputVolume.GetStorageNode().GetFileName())[0] +"_C"+str(vtID) +"_crop_iso.nrrd"  
           print("iso cropped: "+inputCropIsoPath)
           resampleSpacing = " ["+ str(self.RSx) + "," + str(self.RSy) + "," + str(self.RSz) + "] "
           try:
               resamplingCommand = slicer.modules.resamplescalarvolume.path
           except AttributeError:
               #TODO: Get Slicer PATH
               SlicerPath      =os.path.abspath(os.path.join(os.path.abspath(os.path.join(os.sys.executable, os.pardir)), os.pardir))
               SlicerBinPath   = os.path.join(SlicerPath,"Slicer")
               ResampleBinPath =  os.path.join( (glob.glob(os.path.join(SlicerPath,"lib","Slicer") + '*'))[0]    , "cli-modules","ResampleScalarVolume" )
               if sys.platform == 'win32':
                   ResampleBinPath + ".exe"
                   resamplingCommand = SlicerBinPath + " --launch " + ResampleBinPath
               else:
                   #note: in windows, no need to use --launch
                   resamplingCommand = ResampleBinPath + ".exe"
           #endtry
           print(resamplingCommand)
           si = None 
           currentOS = sys.platform           
           cmdPars = " -i linear -s "+ resampleSpacing + inputCropPath +" "+inputCropIsoPath  
           Cmd = resamplingCommand  + cmdPars
           if sys.platform == 'win32':
              #note: in windows, no need to use --launch
              SlicerBinPath = SlicerBinPath +".exe"
              resamplingCommand = ResampleBinPath + ".exe"
              print(os.path.getsize(resamplingCommand))
              si = subprocess.STARTUPINFO()
              si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
              Cmd = resamplingCommand  + cmdPars
              print("Executing ... "+Cmd)
              cRs = subprocess.call(Cmd , shell = (sys.platform == currentOS) , startupinfo=si )
           else:
              print(currentOS)
              print("Executing ... "+Cmd)
              cRs = subprocess.call(Cmd , shell = (sys.platform == currentOS) , startupinfo=si )
           #endif              

           
           #inputCropPath = inputCropIsoPath
           print(" Cropping and resampling are done !!! ")
        #endif

        #inputCropPath    = inputCropPath.strip("'")
        print(" Cropped image is saved in : [%s]" % inputCropPath)
        print(" Cropping is done !!! ") 
        # so we can remove these files later
        return inputCropPath
        
  #--------------------------------------------------------------------------------------------
  #                        run elastix
  #--------------------------------------------------------------------------------------------      
  def runElastix(self, elastixBinPath, fixed, moving, output, parameters, verbose, line):
      print ("************  Compute the Transform **********************")
      currentOS = sys.platform
      print(currentOS)
      #python2
      #Cmd = elastixBinPath + ' -f ' + unicode(fixed, "utf-8") + ' -m ' +  unicode(moving, "utf-8")  + ' -out ' +  unicode(output, "utf-8") + ' -p ' + parameters 
      #python3
      Cmd = elastixBinPath + ' -f ' + fixed + ' -m ' +  moving  + ' -out ' +  output + ' -p ' + parameters 
      errStr="No error!"
#      if hasattr(subprocess, 'mswindows'):
      if currentOS in ["win32","msys","cygwin"]:
          print(" elastix is running in Windows :( !!!") 
          print(Cmd)
          si = subprocess.STARTUPINFO()
          si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
          cTI = subprocess.call(Cmd , shell = (sys.platform == currentOS) , startupinfo=si )
      elif currentOS in ["linux","linux2"]:
          print(" elastix is running in Linux :) !!!") 
          CmdList = Cmd.split()
          print(CmdList)
          cTI=  subprocess.Popen(CmdList, env=self.elastixEnv,  stdout=subprocess.PIPE, universal_newlines=True)
          cTI.wait()
          out, err = cTI.communicate()
          cTI=int(not err==None)
      elif currentOS in ["darwin","os2","os2emx"]:
          print(" elastix is running in Mac :( !!!") 
          CmdList = Cmd.split()
          print(CmdList)
          cTI=  subprocess.Popen(CmdList, env=self.elastixEnv,  stdout=subprocess.PIPE, universal_newlines=True)
          cTI.wait()
          out, err = cTI.communicate()
          cTI=int(not err==None)          
      else:
            print(" elastix is running in Unknown system :( !!!")           
            cTI=1
            errStr="elastix error at line"+ line +", check the log files"
      #endif
      print(cTI)
      #time.sleep(5)      
      self.chkElxER(cTI,errStr) # Check if errors happen during elastix execution
      
      return cTI
  #enddef

  #--------------------------------------------------------------------------------------------
  #                        run transformix
  #--------------------------------------------------------------------------------------------      
  def runTransformix(self,transformixBinPath, img, output, parameters, verbose, line):
      print ("************  Apply transform **********************")
      currentOS = sys.platform
      Cmd = transformixBinPath + ' -tp ' + parameters + ' -in ' +  img  +' -out ' +   output +' -def '+ ' all '              
      #if subprocess.mswindows:
      errStr="No error!"
#      if hasattr(subprocess, 'mswindows'):          
      if currentOS in ["win32","msys","cygwin"]:
         print(" transformix is running in Windows :( !!!") 
         print(Cmd)         
         si = subprocess.STARTUPINFO()
         si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
         cTS = subprocess.call(Cmd , shell = (sys.platform == currentOS) , startupinfo=si )
      elif currentOS in ["linux","linux2"]:
          print(" transformix is running in Linux :) !!!") 
          CmdList = Cmd.split()
          print(CmdList)
          cTS=  subprocess.Popen(CmdList, env=self.elastixEnv,  stdout=subprocess.PIPE, universal_newlines=True)
          cTS.wait()
          out, err = cTS.communicate()
          cTS=int(not err==None)
      elif currentOS in ["darwin","os2","os2emx"]:
          print(" transformix is running in Mac :( !!!") 
          CmdList = Cmd.split()
          print(CmdList)
          cTS=  subprocess.Popen(CmdList, env=self.elastixEnv,  stdout=subprocess.PIPE, universal_newlines=True)
          cTS.wait()
          out, err = cTS.communicate()
          cTS=int(not err==None)
      else:
            print(" elastix is running in Unknown system :( !!!")           
            cTS=1
            errStr="transformix error at line"+ line +", check the log files"
      #endif
      print(cTS)
      self.chkElxER(cTS,errStr) # Check if errors happen during elastix execution
      return cTS
  #enddef
        
  #--------------------------------------------------------------------------------------------
  #                       Check Elastix error
  #--------------------------------------------------------------------------------------------
  # This method checks if errors happen during elastix execution
  def chkElxER(self,c, s):
        if c>0:
           #qt.QMessageBox.critical(slicer.util.mainWindow(),'segmentation', s)
           print(s)  
           return False
        else: 
            print("done !!!")
        #endif
 #enddef             
            
    
  def openResultsFolder(self):
      if not hasattr(self, 'vtVars'):
         self.setGlobalVariables(1)
      #endif 		 
      currentOS = sys.platform
      print(currentOS)
      if currentOS in ["win32","msys","cygwin"]:
         cmd=  subprocess.Popen("explorer " + self.vtVars['outputPath'])
      elif currentOS in ["linux","linux2"]:
         cmd = os.system('xdg-open ' + self.vtVars['outputPath'])
      elif currentOS in ["darwin","os2","os2emx"]:
         cmd = os.system('open ' + self.vtVars['outputPath'])
      else:
         print("uknown system")   
      #endif
  #enddef  

  # segmenteditor effect on the resulted segmentations 
  # this function is called by functions like doSmoothing and doMargining
  def getSegmentationEditor(self,segNode,masterNode):
            # Create segment editor to get access to effects
           segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
           segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
           segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
           segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
           segmentEditorWidget.setSegmentationNode(segNode)
           segmentEditorWidget.setMasterVolumeNode(masterNode)
           return segmentEditorWidget, segmentEditorNode

  # smoothing segmentation effect    
  #TODO: add more options   
  def runSmoothing(self,segNode,masterNode,KernelSizeMm):       
           # Smoothing
           [segEditorW,segEditorN]= self.getSegmentationEditor(segNode,masterNode)               
           for i in range (0,segNode.GetSegmentation().GetNumberOfSegments()):
               segmentID = segNode.GetSegmentation().GetNthSegmentID(i)
               segEditorW.setActiveEffectByName("Smoothing")
               segEditorW.setCurrentSegmentID(segmentID)
               effect = segEditorW.activeEffect()
               effect.setParameter("SmoothingMethod", "MEDIAN")
               effect.setParameter("KernelSizeMm", KernelSizeMm)
               effect.self().onApply()
           #endfor
           # Clean up
           segEditorW = None
           slicer.mrmlScene.RemoveNode(segEditorN)
  #enddef

  # Margin segmentation effect
  # MarginSizeMm>0 Grow, else Shrink          
  def runMargining(self,segNode,masterNode,MarginSizeMm):
           #Dilation and Eroding
           [segEditorW,segEditorN]= self.getSegmentationEditor(segNode,masterNode)               
           for i in range (0,segNode.GetSegmentation().GetNumberOfSegments()):
               segmentID = segNode.GetSegmentation().GetNthSegmentID(i)
               segEditorW.setActiveEffectByName("Margin")
               segEditorW.setCurrentSegmentID(segmentID)
               effect = segEditorW.activeEffect()
               effect.setParameter("MarginSizeMm", MarginSizeMm) 
               effect.self().onApply()
           #endfor
           # Clean up
           segEditorW = None
           slicer.mrmlScene.RemoveNode(segEditorN)
  #enddef  

  def removeOtputsFolderContents(self):
      try:
          for file in os.listdir(self.vtVars['outputPath']):
              filePath = os.path.join(self.vtVars['outputPath'], file)
              if os.path.isfile(filePath):
                 os.remove(filePath)
             #endif
          #endfor                    
      except Exception as e:
            print("nothing to delete ...")
            print(e)
       #endtry 
  #enddefr 
                 
  def removeTmpsFiles(self):
      #remove old files 
      outputPath = self.vtVars['outputPath']
      print("removing temp output files!")
      fds=[]
      outoutputFolders = os.listdir(outputPath)
      for fd in outoutputFolders:        
          if os.path.isdir(os.path.join(outputPath,fd)):
             fds.append(fd)   
          #endif
      #endfor 
      fds.append(".")   
      for fd in fds:         
          print(os.path.join(outputPath,fd) )
          resfiles = os.listdir(os.path.join(outputPath,fd) )
          for fnm in resfiles:
              if "IterationInfo" in fnm:
                 os.remove(os.path.join(outputPath,fd,fnm))
              elif  "result" in fnm:
                 os.remove(os.path.join(outputPath,fd,fnm))
              elif  ".log" in fnm:
                 os.remove(os.path.join(outputPath,fd,fnm))
              elif  "TransformParameters" in fnm:
                 os.remove(os.path.join(outputPath,fd,fnm))
              #endif
          #endfor fnm
      #endfor fd
      vissimPath = self.vtVars['vissimPath']
      cropfiles = os.listdir(vissimPath) 
      try:  
         for fnm in cropfiles:
             if "Crop" in fnm:
                 os.remove(os.path.join(vissimPath, fnm))
             #endif
             if re.search('[C][1-7]', fnm):
                 os.remove(os.path.join(vissimPath, fnm))
             #endif
         #endfor
      except Exception as e:
             print(" Error: can not remove " + fnm)
             print(e)   
      #endtry  
      print("removing temp nodes ...!") 
      nodes = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')
      for f in nodes:
          if "_Crop"  in f.GetName(): slicer.mrmlScene.RemoveNode(f)
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:         
          if "Location" in f.GetName(): f.GetDisplayNode().SetVisibility(False) 
      #endfor      
      """
      #using dictionary doesn ot work: 
      #remove temp nodes
      nodes = slicer.util.getNodes().keys()
      for f in nodes:         
          if ("_Crop"    in f): slicer.mrmlScene.RemoveNode(slicer.util.getNodes()[f]); 
          #if ("Location" in f): slicer.mrmlScene.RemoveNode(slicer.util.getNodes()[f]);
          if ("Location" in f): slicer.util.getNodes()[f].GetDisplayNode().SetVisibility(False);
          #endif
       #endfor
       """
                      
  #enddef  

  def rmvSlicerNode(self,node):
    slicer.mrmlScene.RemoveNode(node)
    slicer.mrmlScene.RemoveNode(node.GetDisplayNode())
    slicer.mrmlScene.RemoveNode(node.GetStorageNode())
  #enddef
  # this can be used only if extension is running from GUI  
  def msgBox(self,txt):
      #msg = qt.QMessageBox()
      #msg.setIcon(qt.QMessageBox.Information)
      #msg.setText("information:")
      #msg.setInformativeText(txt)
      #msg.setWindowTitle("VisSimTools")
      #msg.exec_()
      print(txt)
  #enddef

  def setItemChk(self,itemT, itemChk, itemName, nodes,):
      self.vtVars[itemT] = str(itemChk)
      # Set unvisible 
      for f in nodes:
          if (itemName in f.GetName() ):
             f.GetDisplayNode().SetVisibility(self.s2b(self.vtVars[itemT]))
             print(itemT +" is " + str(itemChk))
            #break
          #endif
      #endfor
      if (itemName =="cochleaSide" ):  
          if itemChk:
             self.vtVars[itemT] = "R"  
          else:
             self.vtVars[itemT] = "L"
          #endif
          #print("side selected is " + self.vtVars[itemT]) 
      #endif

   #enddef
  
  def setVtID(self,idx,inputVolumeNode , inputFiducialNode):
      self.vtVars['vtID']=str(idx)
      print(self.vtVars['vtID']+ " is selected")
      self.inputVolumeNode = inputVolumeNode
      self.inputFiducialNode = inputFiducialNode 
      #remove old points 
        
      # Check if a markup node exists
      newNode = True
      print(self.inputVolumeNode.GetName()+"_vtLocations")
      self.inputFiducialNode = None
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if ((f.GetName() == self.inputVolumeNode.GetName()+"_vtLocations") ):
             #replace  current 
             print("inputFiducialNode exist")
             self.inputFiducialNode = f  
             newNode= False
            #endif
      #endfor      
      if not (self.inputFiducialNode is None):
         ls = slicer.modules.markups.logic()
         ls.SetActiveListID(self.inputFiducialNode)
         print(ls.GetActiveListID())
         noPts = self.inputFiducialNode.GetNumberOfFiducials() 
         newFid= True
         for j in range (0, noPts):
             if self.inputFiducialNode.GetNthFiducialLabel(j)==("C"+self.vtVars['vtID']) :
                newFid= False 
                print("C"+self.vtVars['vtID'] +" exist, removing old point at: " +str(j))
                #get the new location
                self.inputFiducialNode.RemoveMarkup(j)      
             #endif
         #endfor
      else:
         print("inputFiducialNode does not exist")
      #endif        
  #enddef
      
  # check if vertebra location is available
  def setVtIDfromEdt(self,point,vtID):
        # no external call
        #print("no external call, point= " + point)
        self.vtVars['vtID'] = str(vtID)
        isExternalCall = False
        print("point changed,  " + str(vtID) + " is selected")      
        #TODO: add option to use point from text for cropping   
        return isExternalCall     
  #enddef
      
  #this part need to be optimized
  # reg =0: no registration, 1: fixed image, 2: moving image
  def locateItem(self, inputVolumeNode, inputPointEdt, reg, vtID):
      self.inputFiducialNodes = []
      for i in range (3):
          self.inputFiducialNodes.append(slicer.vtkMRMLMarkupsFiducialNode())
          # 0 input: 1: fixed, 2:moving
      #endfor
      if reg ==1:
         regType="F"
      elif reg ==2:
         regType="M"
      else:
         regType=""
      #endif 
      if not hasattr(self, 'vtVars'):
         self.setGlobalVariables(int(not (vtID==0)) )
      #end             
      self.inputVolumeNode   = inputVolumeNode
      self.inputPointEdt =inputPointEdt
      self.vtVars['vtID']= str(vtID)
      # Reset global point label
      self.inputPoint = [0,0,0]
      inputPointEdt.setText(str(self.inputPoint))
      # Check if a volume is selected
      if not self.inputVolumeNode:
           print >> sys.stderr, "You need to pick a input volume first before locating vertebra."
           return -1
      #endif
      #  Display suitable during locating the vertebra
      disp_logic = slicer.app.layoutManager().sliceWidget(self.vtVars['dispViewTxt']).sliceLogic()
      disp_cn = disp_logic.GetSliceCompositeNode()
      disp_cn.SetBackgroundVolumeID(self.inputVolumeNode.GetID())
      lm = slicer.app.layoutManager()
      lm.setLayout(int(self.vtVars['dispViewID']))
      # Fit slice to window
      sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
      layoutManager = slicer.app.layoutManager()
      for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()
            #endif
      #endfor

      if vtID==0: #Cochlea:
         print(" ..... getting Cochlea location in the input image")  
         self.FidLabel =  regType+"_CochleaLocation"       
      else: #cervical spine
         # redefine to be used in the logic class.
         print(" ..... getting vertebra location in the input image")  
         self.FidLabel= regType+"_vtLocations"
      #endif        
      # Check if a markup node exists
      newNode = True
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
             if ((f.GetName() == inputVolumeNode.GetName()+self.FidLabel ) ):
                  #replace  current 
                  self.inputFiducialNodes[reg] = f  
                  newNode= False
              #endif
      #endfor
      #this is implemented in setVtID
      if newNode:
            self.inputFiducialNodes[reg] = slicer.vtkMRMLMarkupsFiducialNode()
            self.inputFiducialNodes[reg].SetName(inputVolumeNode.GetName()+self.FidLabel)
            slicer.mrmlScene.AddNode(self.inputFiducialNodes[reg])
      #endif     
      self.inputFiducialNodes[reg].GetDisplayNode().SetVisibility(True);        
      self.inputFiducialNodes[reg].GetDisplayNode().SetTextScale(2)
      self.inputFiducialNodes[reg].GetDisplayNode().SetSelectedColor(1,0,0)           
      # Start Fiducial Placement Mode in Slicer for one node only
      placeModePersistance = 0 # one node only
      slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

      # Observe scene for updates
      self.addObs = self.inputFiducialNodes[reg].AddObserver(self.inputFiducialNodes[reg].MarkupAddedEvent,   self.onInputFiducialNodeMarkupAddedEvent)
      self.modObs = self.inputFiducialNodes[reg].AddObserver(self.inputFiducialNodes[reg].PointModifiedEvent, self.onInputFiducialNodePointModifiedEvent)
      self.rmvObs = self.inputFiducialNodes[reg].AddObserver(self.inputFiducialNodes[reg].MarkupRemovedEvent, self.onInputFiducialNodeMarkupRemovedEvent)
      return  self.inputFiducialNodes[reg]   
  #enddef

  #--------------------------------------------------------------------------------------------
  #    InputFiducialNode MarkupAddedEvent
  #--------------------------------------------------------------------------------------------
  def onInputFiducialNodeMarkupAddedEvent(self, caller, event):
      # it seems this action happened after adding new fiducial
      print("Fiducial adding event!")
      #remove previous observer
      caller.RemoveObserver(self.addObs)
      noPts = caller.GetNumberOfFiducials() 
      rasPt = [0,0,0] 
      if self.vtVars['vtID']== "0":  #Cochlea
         caller.SetNthFiducialLabel(noPts-1, self.FidLabel)
         caller.GetNthFiducialPosition(noPts-1,rasPt)
      else: #Spine
         caller.SetNthFiducialLabel(noPts-1, "C"+self.vtVars['vtID'])
         caller.GetNthFiducialPosition(noPts-1,rasPt)
      #self.inputPoint = self.vsc.ptRAS2IJK(caller, rasPt,self.inputVolumeNode,noPts-1)
      self.inputPoint = self.ptRAS2IJK( caller,self.inputVolumeNode,noPts-1)       
      self.inputPointEdt.setText(str(self.inputPoint))
      print(" ..... location RAS: " + str(rasPt))  
      print(" ..... location in the input image set to: " + str(self.inputPoint))
  #enddef  
  
  #--------------------------------------------------------------------------------------------
  #    InputFiducialNode PointModifiedEvent
  #--------------------------------------------------------------------------------------------
  # The fiducial point saved in RAS, we need to convert to IJK
  #  more info in our wiki 
  def onInputFiducialNodePointModifiedEvent(self, caller, event):
      #caller.RemoveObserver(self.modObs)
      # get the new IJK position and display it
      rasPt = [0,0,0] 
      i = caller.GetAttribute('Markups.MovingMarkupIndex')
      if not (i is None):
         i=int(i)
         caller.GetNthFiducialPosition(i,rasPt)
         self.inputPoint = self.ptRAS2IJK(caller, self.inputVolumeNode, i)        
         self.inputPointEdt.setText(str(self.inputPoint))
      #endif   
  #enddef

  def onInputFiducialNodeMarkupRemovedEvent(self, caller, event):
      #print("Fiducial removed event!")
      caller.RemoveObserver(self.rmvObs)
      #i = caller.GetNumberOfFiducials()-1
      #print("number of rmaining fiducials: " + str(i))
  #endif

  #--------------------------------------------------------------------------------------------
  #                        Calculate Segmentation Information    
  #--------------------------------------------------------------------------------------------
  def getItemInfo(self, segNode, masterNode, tblNode, vtID):
        segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
        segStatLogic.getParameterNode().SetParameter("Segmentation", segNode.GetID())
        segStatLogic.getParameterNode().SetParameter("ScalarVolume", masterNode.GetID())
        segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.enabled","False")
        segStatLogic.getParameterNode().SetParameter("ScalarVolumeSegmentStatisticsPlugin.voxel_count.enabled","False")
        segStatLogic.computeStatistics()
        if vtID == 0:  # Cochlea
           # compute the fiducial length
           segStatLogic.exportToTable(tblNode)           
           tblNode.SetCellText(0,2,self.vtVars['StLength'])           
           tblNode.SetCellText(1,2,"0")   
        else           :  #Spine
           #  egt volume size of a vertebra
           if tblNode is None:
              print("create new table  .........................................")
              tblNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
              #tblNode.SetName(masterNode.GetName()[0:-7]+"_tbl")
              #this work for SpineTools
              #TODO: check if it works with vertebra tools
              tblNode.SetName(masterNode.GetName()[0:-3]+"_tbl")

              for i in range(5):
                 tblNode.AddColumn()
              #endfor 
              tblNode.GetTable().GetColumn(0).SetName("Vertebra")
              tblNode.GetTable().GetColumn(1).SetName("Volume mm3")
              tblNode.GetTable().GetColumn(2).SetName("CoM X")
              tblNode.GetTable().GetColumn(3).SetName("CoM Y")
              tblNode.GetTable().GetColumn(4).SetName("CoM Z")
              tblNode.Modified()
           #endif

           #remove old row of this vertebra if exists
           for i in range (tblNode.GetNumberOfRows()):
               if "C"+str(vtID) == tblNode.GetCellText(i,0):
                  print( "C"+str(vtID) + " table row exists, old values will be removed in row." + str(i)) 
                  tblNode.RemoveRow(i)   
               #endif
           #endfor
           tblNode.AddEmptyRow() # empty row for current vertebra info    

           spTmpTblNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
           spTmpTblNode.SetName("tmpTable")
           segStatLogic.exportToTable(spTmpTblNode)        
           idx = tblNode.GetNumberOfRows()-1 # this should be the empty row we added 
           print(" last row index: "+ str(idx))
           tblNode.GetTable().SetRow(  idx  ,  spTmpTblNode.GetTable().GetRow(0)  ) 
           tblNode.SetCellText(idx,2," ") 
           tblNode.SetCellText(idx,3," ") 
           tblNode.SetCellText(idx,4," ") 
           slicer.mrmlScene.RemoveNode(spTmpTblNode )

           # if this is C7 compute center of mass        
           if (vtID ==7) and (self.vtVars['vtMethodID']== "0"): # for testing
              print("updating COM in table ..............")     
              segID = segNode.GetSegmentation().GetSegmentIdBySegmentName("C"+str(vtID))
              modelNode = segNode.GetClosedSurfaceRepresentation(segID)
              com = vtk.vtkCenterOfMass(); com.SetInputData(modelNode);   com.Update()
              segNodeCoM = com.GetCenter()
              tblNode.SetCellText(idx,2,str(segNodeCoM[0])) 
              tblNode.SetCellText(idx,3,str(segNodeCoM[1])) 
              tblNode.SetCellText(idx,4,str(segNodeCoM[2]))
              self.vtVars['segNodeCoM']=str(segNodeCoM)
              #if  tblNode.GetCellText(tblNode.GetNumberOfRows(),0)=="":
              #    tblNode.RemoveRow(tblNode.GetNumberOfRows()-1)
              #endif
           #endif
        #endifelse
        # update table and set it the active table in slicer
        print("updating table ..............") 
        tblNode.Modified()
        slicer.app.applicationLogic().GetSelectionNode().SetActiveTableID(tblNode.GetID())
        slicer.app.applicationLogic().PropagateTableSelection()
        return tblNode
  #enddef
  
  #--------------------------------------------------------------------------------------------
  #                        Calculate length and volume of scalas
  #--------------------------------------------------------------------------------------------
  # This function compute the distance between all the fiducials in a markupnode       
  def getFiducilsDistance(self,markupsNode,tblNode):
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
        self.vtVars['StLength'] = str(markupsDistance)
        if not tblNode is None: 
           tblNode.SetCellText(0,2,self.vtVars['StLength'])
        #endif
           
        return markupsDistance
  #enddef  
  def fitAllSlicesViews(self):
      sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
      layoutManager = slicer.app.layoutManager()
      for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()
            #endif
      #endfor
  #enddef

  # An option to control results displaying

  def fuseWithOutColor(self, disableColor):
      print ("Disable enabled is " + str(disableColor))
      firstNode = slicer.util.getNode(slicer.app.layoutManager().sliceWidget(slicer.app.layoutManager().sliceViewNames()[0]).sliceLogic().GetSliceCompositeNode().GetForegroundVolumeID())
      secondNode = slicer.util.getNode(slicer.app.layoutManager().sliceWidget(slicer.app.layoutManager().sliceViewNames()[0]).sliceLogic().GetSliceCompositeNode().GetBackgroundVolumeID())
      self.fuseTwoImages( firstNode, secondNode , not disableColor)            
  #enddef

  #fuse two images with two different colors
  def fuseTwoImages(self, firstNode, secondNode, colorful):
      if colorful:
         self.vtVars['nodeColorFG']          = slicer.modules.colors.logic().GetColorTableNodeID(20)  # green color
         self.vtVars['nodeColorBG']          = slicer.modules.colors.logic().GetColorTableNodeID(16)  # magnta color
      else:
         self.vtVars['nodeColorFG']          = slicer.modules.colors.logic().GetColorTableNodeID(1)  # gray color
         self.vtVars['nodeColorBG']          = slicer.modules.colors.logic().GetColorTableNodeID(1)  # gray color
      #endif
      #TODO: replace this with a loop or short code
      for i in range(3):     
          print(i) 
          viewSliceCompositeNode = slicer.app.layoutManager().sliceWidget(slicer.app.layoutManager().sliceViewNames()[i]).sliceLogic().GetSliceCompositeNode()
          viewSliceCompositeNode.SetBackgroundVolumeID(firstNode.GetID())
          viewSliceCompositeNode.SetForegroundVolumeID(secondNode.GetID())
          viewSliceCompositeNode.SetForegroundOpacity(0.5)
      #endfor

      # The layout is set to show only the default view.
      slicer.app.layoutManager().setLayout(int(self.vtVars['dispViewID']))                             
      # The window level for each image is set to be the same value.
      firstNode.GetScalarVolumeDisplayNode().AutoWindowLevelOff()
      firstNode.GetScalarVolumeDisplayNode().SetWindowLevel(1000, 400)
      secondNode.GetScalarVolumeDisplayNode().AutoWindowLevelOff()
      secondNode.GetScalarVolumeDisplayNode().SetWindowLevel(1000, 400)

      # Assign colors, note that original color value is 1.
      firstNode.GetDisplayNode ().SetAndObserveColorNodeID(self.vtVars['nodeColorFG']) # green color table
      secondNode.GetDisplayNode().SetAndObserveColorNodeID(self.vtVars['nodeColorBG']) # magenta color table

      # Fit slices to window
      self.fitAllSlicesViews() 
  #enddef

  def dispSeg(self,inputVolumeNode, vtSegNode, view):
        lm = slicer.app.layoutManager();   
        lm.setLayout(view)
        r_logic = lm.sliceWidget("Red").sliceLogic()
        r_cn = r_logic.GetSliceCompositeNode()
        r_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())
        y_logic = lm.sliceWidget("Yellow").sliceLogic()
        y_cn = y_logic.GetSliceCompositeNode()
        y_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())
        g_logic = lm.sliceWidget("Green").sliceLogic()
        g_cn = g_logic.GetSliceCompositeNode()
        g_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())

        #center 3D view and zoom in 3 times
        v3DDWidget = lm.threeDWidget(0)
        v3DDWidgetV = v3DDWidget.threeDView()
        v3DDWidgetV.resetFocalPoint() 
        v3DDWidgetV.zoomFactor =3
        v3DDWidgetV.zoomIn()
        v3DDWidgetV.zoomFactor =0.05 # back to default value
  #enddef
    
    
    
class VisSimCommonTest(ScriptedLoadableModuleLogic):

  def setUp(self):
    slicer.mrmlScene.Clear(0)   
  #enddef

  def runTest(self):
    self.setUp()
    print(VisSimCommonLogic().tstSum(10,20))
  #enddef
#enclass 
