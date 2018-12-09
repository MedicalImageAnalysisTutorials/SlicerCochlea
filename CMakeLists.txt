cmake_minimum_required(VERSION 3.5)

project(SlicerCochlea)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SlicerCochlea")
set(EXTENSION_CATEGORY "VisSimTools")
set(EXTENSION_CONTRIBUTORS "Ibraheem Al-Dhamari (VisSim, Koblenz-Landau University)")
set(EXTENSION_DESCRIPTION "The cochlea is a very important part of the inner ear. It is responsible for the transfer of audio signals to the brain. There are two modules in this extension. The first module, Cochlea Registration, registers and fuses cochlea images from different modalities. The second plugin, Cochlea Segmenttaion, segment out the cochlea structure, scala tympani and other scalae (media and vestibuli) in addition to measure the length and the size of the scala tympani. Both modules required a few seconds to complete the tasks. They use a cusomised set of parameters for elastix toolbox. The extension will download elastix binaries and other necessary files so it will take some minutes when use first time. For testing, cochlea sample datasets can be downloaded using Slicer Data Store module. You are welcome to contribute by optimising the code or sharing your cochlea dataset. For questions or problems using this extension please post in the gitHub issue or in Slicer forum. For more details please read the related publications:
Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Rania Helal, Friedrich Lisseck and Roland Jacob, (2018), Automatic Cochlear Length and Volume Size Estimation, Accepted in: First International Workshop on Context-Aware Operating Theater OR 2, MICCAI 2018, Granada Spain.
Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Roland Jacob, (2018), Automatic Cochlea Multi-modal Images Segmentation Using Adaptive Stochastic Gradient Descent. In: CI2018 DC Emerging Issues in Cochlear Implantation, Washington DC, USA. (link).
Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus and Roland Jacob, (2017): Automatic Cochlea Segmentation Using Diffusion Snakes. In: 15th Symposium on Cochlear Implants in Children, CI 2017, San Francisco, USA (link).
Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Friedrich Lisseck and Roland Jacob, (2017): ACIR: automatic cochlea image registration. In: Proceedings SPIE Medical Imaging 2017: Image Processing;. SPIE. Bd. 10133. S. 10133p1-10133p5 (link).
Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Friedrich Lissek, Roland Jacob(2016): Automatic Multimodal Registration and Fusion of 3D Human Cochlea Images. In: 14th International Conference on Cochlear Implants, Toronto, Canada. S. 511-514. ")
#set(EXTENSION_ICONURL "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/SlicerCochlea.png")
#set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/c.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/r1.png  https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/r2.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s1.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s2.png")
set(EXTENSION_ICONURL "https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/raw/master/SlicerCochlea.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/c.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s2.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/r1.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s1.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCochlea/master/Screenshots/s2.png")

set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

set(EXTENSION_DEPENDS SlicerElastix)

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CochleaReg)
add_subdirectory(CochleaSeg)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
