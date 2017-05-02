This directory has all the codes needed to implement my project for CS766
All images are saved and read as .nrrd files, that can be loaded using Slicer 3D software

A list of the codes and their functions are as follows:

main.m
Script that is used for locating, segmenting, and placing liver sphere of test patients. This is the only script that is needed to implement the entire process. 

It uses the following codes:
nrrdWriter (to save images as nrrd files)
nrrdread (to load nrrd files)
loadAmOrNrrd (this calls the nrrdread function and loads in a more digestible format)
FindLiver (this generates points of interest and runs them through the SVM model)
and
liverSegment (this segments the liver using a region growing algorithm)


THE DATA
Two example patients, as well as the training data for the SVM model, can be downloaded at the following link:
https://uwmadison.box.com/s/1zpsux7we41zjfvhmffmcv02lrfl1sy5
