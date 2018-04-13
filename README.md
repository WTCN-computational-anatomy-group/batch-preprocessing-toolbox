# Batch pre-processing of neuroimaging data using SPM

This code performs batch pre-processing of neuroimaging data using the SPM software.

## Folder structure
This code requires a well-defined folder structure of the input data. What follows is an example of this structure. 

Let's say we have eight subjects who each have T1-, T2- and PD-weighted MR data; and an image with categorical labels of tissue types. If the root directory is called *data1* then the folder structure should be:

|-- *data1*  
|-- |-- *S1*  
|-- |-- |-- *scans*  
|-- |-- |-- | -- *T1*  
|-- |-- |-- | -- *T2*  
|-- |-- |-- | -- *PD*  
|-- |-- |-- *labels*    
.  
.  
.  
|-- |-- *S8*  
|-- |-- |-- *scans*  
|-- |-- |-- | -- *T1*  
|-- |-- |-- | -- *T2*  
|-- |-- |-- | -- *PD*  
|-- |-- |-- *labels*   

The folder that contains the imaging data for each subjects must be named *scans* and the folder that contains the labels *labels*. The categorical labels are assumed to be stored in one image.

## Basic use cases

Four basic use cases are in the main function *preprocess_images.m*. It uses data available in the shared folder *Ashburner_group*.

### Example 1 (T1-weighted data)
* Rigidly realign to MNI space
* Remove data outside of head (+neck)
* Bias-field correct
* Skull strip
* Make ML-labels
* Normalise intensities
* Write 2D versions

### Example 2 (T1, T2-, PD-weighted data)
* Rigidly realign to MNI space
* Remove data outside of head
* Co-register
* Reslice to image with largest FOV
* Make 1 mm isotropic voxels
* Write 2D versions

### Example 3 (T1-weighted data)
* Rigidly realign to MNI space
* Remove data outside of head (+neck)
* Segment
* Write 2D versions

### Example 4 (CT data)
* Rigidly realign to MNI space
* Remove data outside of head (+neck)
* Skull strip
* Write 2D versions

## Dependencies

This project has strong dependencies to SPM12 and its `Shoot` toolbox. Both of them should be added to Matlab's path. SPM can be downloaded at [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/).

Core functions also depend on our [`auxiliary-functions` toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions), which gathers lots of low-level functions.

Furthermore, executable scripts depend on our [`distributed-computing` toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing), which helps parallelising parts of the code either on the local workstation (using Matlab's parallel processing toolbox) or on a computing cluster (see the toolbox help file for use cases and limitations).

## TODO

* Improve functionality for normalising MRI intensities
* Compare IXI co-reg vs non-co-reg
