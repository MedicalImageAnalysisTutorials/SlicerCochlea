# Cochlea Image Analysis

<img src="https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/blob/master/Cochlea.png" width="400" height="400">

This is a [3D Slicer](https://gaithub.com/Slicer/Slicer) plugin that uses [elastix toolbox](https://github.com/SuperElastix/elastix) for Multi-modal cochlea Images registration, segmentation and analysis. The elastix parameters file can be found [here](https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/tree/master/docs/elastixPars)

# Tested on

* Slicer 5.0.3 [for windows](https://slicer-packages.kitware.com/api/v1/file/62d5d2ebe911182f1dc285b2/download) and [for linux](https://slicer-packages.kitware.com/api/v1/file/62cc52d2aa08d161a31c1af2/download). Tested on Windows 10 and Ubuntu 20.04

* Slicer 4.10.2 [for windows](https://slicer-packages.kitware.com/api/v1/file/60add732ae4540bf6a89c029/download) and [for linux](https://slicer-packages.kitware.com/api/v1/file/60add73aae4540bf6a89c03b/download). Tested on Windows 10 and Ubuntu 18.04


This project contains two modules:

  1. Cochlea image registration.
  2. Cochlea image segmentation and measuerments.

Notice: one can modify the optimiser in elastix parameters to get new results. 

# How to use:

* In case you are new to 3D Slicer check my [Slicer tutorials](https://www.youtube.com/playlist?list=PLW9iOMxMvikpyCUMmuqiloNp7rpaUl2M1) (notice that there maybe small difference between different versions):
  - [3D Slicer installation](https://www.youtube.com/watch?v=7XHhgpk0m78&list=PLW9iOMxMvikpyCUMmuqiloNp7rpaUl2M1&index=3).
  - [3D Slicer introduction](https://www.youtube.com/watch?v=mmf5eb0WrR8&list=PLW9iOMxMvikpyCUMmuqiloNp7rpaUl2M1&index=4).
* See the [documentation](https://medicalimageanalysistutorials.github.io/SlicerCochlea/)
* See the video tutorials for [cochlea registration and fusion](https://www.youtube.com/watch?v=JfEaPO3N47U&t=4s) and [cochlea segmentation and analysis](https://www.youtube.com/watch?v=A_mTcT3eT_c&t=1s)

# Publications:

Please cite our papers:

*  Al-Dhamari I, Helal R, Morozova O, Abdelaziz T, Jacob R, et al. (2022) Automatic intra-subject registration and fusion of multimodal cochlea 3D clinical images. PLOS ONE 17(3): e0264449. https://doi.org/10.1371/journal.pone.0264449. [Link](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0264449&type=printable)
*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Roland Jacob, (2018), Automatic Cochlea Multi-modal Images Segmentation Using Adaptive Stochastic Gradient Descent. In: CI2018 DC Emerging Issues in Cochlear Implantation, Washington DC, USA.
*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Friedrich Lisseck and Roland Jacob, (2017): ACIR: automatic cochlea image registration. In: Proceedings SPIE Medical Imaging 2017: Image Processing;. SPIE. Bd. 10133. S. 10133p1-10133p5. [link](http://spie.org/Publications/Proceedings/Paper/10.1117/12.2254396)



# Updates

Updtaed on 6 August 2022

* Links are updated 
* Master branch is now compitable with Slicer 5 
* Tested on Slicer 5.0.3 on Ubuntu 20.04 and Windows 10

# Notes:  

* For general 3D Slicer questions, please use Slicer [forum](https://discourse.slicer.org), many experts will be able to help you there. 
* For Slicer cochlea related questions, comemnts, or feedback, please use github [discussion](https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/discussions/categories/q-a) section. 
* For  Slicer cochlea related bugs, please open a new [issue](https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/issues) if needed. Please mention your operating system, slicer version, and the error message you see in python interactor.
* For sharing private information, dataset, a future project proposal or cooperation, please use the [email](ibr_ex@yahoo.com), use SlicerCochlea in the subject. 
* It would be nice to share your anonymised cochlea dataset.
* It would be nice to contribute to this repository!