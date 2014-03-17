Pluto Analyzer
==============

Simple tool to analyze ROOT trees created with the Pluto event generator. 

Build process
-------------

To compile the code, [ROOT](http://root.cern.ch/ "ROOT") has to be installed. For building use the following command:
``g++ -std=gnu++11 -o main main.cpp `root-config --cflags --glibs` -lSpectrum``
