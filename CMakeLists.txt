#-----------------------------------------
#          SlicerCochlea master
#          updated: 16.6.2019
#-----------------------------------------
cmake_minimum_required(VERSION 3.5)

project(SlicerCochlea)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://medicalimageanalysistutorials.github.io/SlicerCochlea")
set(EXTENSION_CATEGORY "VisSimTools")
set(EXTENSION_CONTRIBUTORS "Ibraheem Al-Dhamari (VisSim, Koblenz-Landau University)")
set(EXTENSION_DESCRIPTION [=[The cochlea is a very important part of the inner ear. It is responsible for the transfer of audio signals to the brain.

There are two modules in this extension:
<ul>
  <li><b>Cochlea Registration</b>: Register and fuse cochlea images from different modalities.</li>
  <li><b>Cochlea Segmentation</b>: Segment out the cochlea structure, scala tympani and other scalae (media and vestibuli) and measure the length and the size.</li>
</ul>

Both modules required a few seconds to complete the tasks. They use a customized set of parameters from <a href="https://github.com/SuperElastix/elastix">elastix toolbox</a> (itself using <a href="https://itk.org/">Insight Toolkit (ITK)</a>). The extension will download elastix binaries and other necessary files so it will take some minutes when use first time. For testing, cochlea sample datasets can be downloaded using Slicer Data Store module.</p>

You are welcome to contribute by optimising the code or sharing your cochlea dataset. For questions or problems using this extension please post in the <a href="https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/issues">GitHub issue tracker</a> or on the <a href="https://discourse.slicer.org/">Slicer Discussion Forum</a>.

For more details please read the related publications:
<ol>
  <li>Al-Dhamari I et al. (2022) Automatic Cochlea Multimodal 3D Image Segmentation And Analysis Using Atlas-model-based Method (to be published!).</li>
  <li>Al-Dhamari I, Helal R, Morozova O, Abdelaziz T, Jacob R, et al. (2022) Automatic intra-subject registration and fusion of multimodal cochlea 3D clinical images. PLOS ONE 17(3): e0264449. https://doi.org/10.1371/journal.pone.0264449. Link</li>
  <li>Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Rania Helal, Friedrich Lisseck and Roland Jacob, (2018), Automatic Cochlear Length and Volume Size Estimation, Accepted in: First International Workshop on Context-Aware Operating Theater OR 2 MICCAI 2018, Granada Spain.</li>
  <li>Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Roland Jacob, (2018), Automatic Cochlea Multi-modal Images Segmentation Using Adaptive Stochastic Gradient Descent. In: CI2018 DC Emerging Issues in Cochlear Implantation, Washington DC, USA. <a href="https://cdn.ymaws.com/www.acialliance.org/resource/resmgr/CI2018/CI2018_PPT/CI2018_DC_Podium_Abstracts.pdf">(link).</a></li>
  <li>Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus and Roland Jacob, (2017): Automatic Cochlea Segmentation Using Diffusion Snakes. In: 15th Symposium on Cochlear Implants in Children, CI 2017, San Francisco, USA <a href="https://cdn.ymaws.com/www.acialliance.org/resource/resmgr/ci2017/CI2017_PPT/CI2017_Podium_Abstracts.pdf">(link).</a></li>
  <li>Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Friedrich Lisseck and Roland Jacob, (2017): ACIR: automatic cochlea image registration. In: Proceedings SPIE Medical Imaging 2017: Image Processing;. SPIE. Bd. 10133. S. 10133p1-10133p5 <a href="http://spie.org/Publications/Proceedings/Paper/10.1117/12.2254396">(link).</a></li>
  <li>Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Friedrich Lissek, Roland Jacob(2016): Automatic Multimodal Registration and Fusion of 3D Human Cochlea Images. In: 14th International Conference on Cochlear Implants, Toronto, Canada. S. 511-514.</li>
</ol>]=])
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/SlicerCochlea.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/c.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/r1.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/r2.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s1.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s2.png")

set(EXTENSION_DEPENDS SlicerElastix) # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})



#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CochleaReg)
add_subdirectory(CochleaSeg)
add_subdirectory(VisSimCommon)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
