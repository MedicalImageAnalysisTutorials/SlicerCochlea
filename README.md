**Cochlea Image Analysis**

<img src="https://github.com/MedicalImageAnalysisTutorials/SlicerCochlea/blob/master/Cochlea.png" width="400" height="400">

This is a [3D Slicer](https://gaithub.com/Slicer/Slicer) plugin that uses [elastix toolbox](https://github.com/SuperElastix/elastix) for Multi-modal cochlea Images registration, segmentation and analysis. More information can be found [here](https://mtixnat.uni-koblenz.de). The elastix parameters file can be found [here](http://elastix.bigr.nl/wiki/index.php/Par0053)

**Tested on:**
Slicer 4.10, Windows 10 and Ubuntu 18.04

This project contains two modules:

  1. Cochlea image registration.
  2. Cochlea image segmentation and measuerments.

Please cite our papers:
*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Rania Helal, Friedrich Lisseck and Roland Jacob, (2018), Automatic Cochlear Length and Volume Size Estimation, Accepted in: First  International Workshop on Context-Aware Operating Theater OR 2, MICCAI 2018, Granada Spain.[link](https://or20.univ-rennes1.fr/sites/or20.univ-rennes1.fr/files/asset/document/aldhamarietal2018_0.pdf)

*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Roland Jacob, (2018), Automatic Cochlea Multi-modal Images Segmentation Using Adaptive Stochastic Gradient Descent. In: CI2018 DC Emerging Issues in Cochlear Implantation, Washington DC, USA.

*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, Friedrich Lisseck and Roland Jacob, (2017): ACIR: automatic cochlea image registration. In: Proceedings SPIE Medical Imaging 2017: Image Processing;. SPIE. Bd. 10133. S. 10133p1-10133p5. [link](http://spie.org/Publications/Proceedings/Paper/10.1117/12.2254396)

Please share your cochlea dataset with us.

Your contribution is welcome!


Updtae : 6/12/2018

- Logic and testing classes are added. Some useful functions can be called from external modules. User can test the module using Reload and test button.
- Runtime libraries are included.
- Suppport DICOM and other formats.
- Support filenames with spaces.

Tested on Slicer 4.10, Ubuntu 18.04 and Windows 10
