# AutoTube

AutoTube quantifies vascular parameters such as the vessel area, width, skeleton length or the number of branching points 
of vascular networks in tissues or in in vitro assays. AutoTube is freely available, comprises an intuitive Graphical 
User Interface (GUI) and helps to perform otherwise highly time-consuming image analyses in a rapid, automated objective 
and reproducible manner.

This repository contains the source files needed to run AutoTube and also a manual to understand AutoTube's functionality.

If you use AutoTube, we appreciate it if you cite an appropriate subset of the following papers:

@article{einstein,
    author =       "Javier A. Montoya-Zegarra and Erica Russo and Peter Runge and Maria Jadhav and Ann-Helen Willrodt and Szymon Stoma and Simon F. Nørrelykke and Michael Detmar and Cornelia Halin",
    title =        "AutoTube: a novel software for the automated morphometric analysis of vascular networks in tissues",
    journal =      "Journal of Angiogenesis",
    volume =       "X",
    number =       "X",
    pages =        "X--X",
    year =         "2018",
    DOI =          "http://dx.doi.org/X/X"
}
 

as well as, the following papers that are used as part of this library:

@inproceedings{DollarICCV13edges,
  author    = {Piotr Doll\'ar and C. Lawrence Zitnick},
  title     = {Structured Forests for Fast Edge Detection},
  booktitle = {ICCV},
  year      = {2013},
}

@inproceedings{ZitnickECCV14edgeBoxes,
  author    = {C. Lawrence Zitnick and Piotr Doll\'ar},
  title     = {Edge Boxes: Locating Object Proposals from Edges},
  booktitle = {ECCV},
  year      = {2014},
}

###################################################################

1. Installation.

 - This code is written for the Matlab interpreter (tested with versions R2017a) and requires the Matlab Image Processing Toolbox. 

 - To use AutoTube follow the instructions described in https://github.com/autotubularity/autotube/wiki

 - Finally, optionally download the Lymphatic Vessel dataset (used in the paper):
 https://doi.org/10.3929/ethz-b-000262426
 After downloading the dataset, you can run AutoTube with the same parameters described in the paper. You should obtain the same statistics.

###################################################################

2. Getting Started.

 - Make sure to carefully follow the installation instructions above.
 - For any potential issues with the software, you can report an issue: https://github.com/autotubularity/autotube/issues
 
 ###################################################################

3. Wiki.

 - We have created a wiki page, where you can find more information on how to use [AutoTube](https://github.com/autotubularity/autotube/wiki).

###################################################################

4. History.

Version 1.0 (05/08/2018)
 - initial version corresponding to the Journal of Angiogenesis
 - now hosting on github (https://github.com/autotubularity/autotube)

###################################################################
