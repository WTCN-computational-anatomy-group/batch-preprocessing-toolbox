# Batch pre-processing of neuroimaging data using SPM

This code performs batch pre-processing of neuroimaging data using the SPM software.

## Basic use cases

### Example one

## Dependencies

This project has strong dependencies to SPM12 and its `Shoot` toolbox. Both of them should be added to Matlab's path. SPM can be downloaded at [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/).

Core functions also depend on our [`auxiliary-functions` toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions), which gathers lots of low-level functions.

Furthermore, executable scripts depend on our [`distributed-computing` toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing), which helps parallelising parts of the code either on the local workstation (using Matlab's parallel processing toolbox) or on a computing cluster (see the toolbox help file for use cases and limitations).

## TODO

* Add functionality to DICOM convert (with new DICOM convert code)
* Include labels in pre-processing (when applicable)
* Add functionality to normalise MRI intensities
* Compare IXI co-reg vs non-co-reg
