# SeisElastic2D
SeisElastic2D is an open-source package for multiparameter FWI in isotropic-, anisotropic- and visco-elastic media built upon several exiting packages including SPECFEM2D, seisDD, etc. Note that the SPECFEM2D package in SeisElastic2D has been modifed to calculate the sensitivity kernels in isotropic-, VTI-, TTI- and visco-elastic media with various model parameterizations.

This package is developed on parallel computing cluster with OpenMPI/MPI. It provides different misfit functions, model parameterizations, advanced optimization methods, flexible inversion strategies and workflows, etc.  These features make it promising to overcome the difficulties in multiparameter elastic FWI and bridge the gap between adademic studies and industrial application. With this package, we have applied isotropic-, VTI- and visco-elastic FWI to practical seismic data successfully.





For downloading and using this package for your own studies, please cite the following reference publications:

Pan, W. & Innanen, K., 2019. Parameterization analysis and field validation of VTI-elastic full-waveform inversion in a walk-away vertical seismic profile configuration, Geophysics. submitted.

Pan, W. & Innanen, K., 2019. Amplitude-based misfit functions in viscoelastic full-waveform inversion applied to walk-away vertical seismic profile data, Geophysics. submitted.

Pan, W., Innanen, K., Geng Y. & Li, J., 2019. Interparameter tradeoff quantification for isotropic-elastic full-waveform inversion with various model parameterizations, Geophysics. 84, R185-R206.

Pan, W., Geng, Y. & Innanen, A. K., 2018. Interparameter tradeoff quantification and reduction in isotropic-elastic full-waveform inversion: synthetic experiments and Hussar land dataset application, Geophysical Journal International. 213, 1305-1333.

Pan, W., Innanen, K. & Geng, Y., 2018. Elastic full-waveform inversion and parameterization analysis applied to walk-away vertical seismic profile data for unconventional (heavy oil) reservoir characterization, Geophysical Journal International. 213, 1934-1968.

Yuan, O. Y., Simons, F. J. & Tromp, J., 2016. Double-difference adjoint seismic tomography, Geophysical Journal International. 206, 1599-1618.

Yuan, O. Y., Simons, F. J. & Bozdag, E., 2015. Multiscale adjoint waveform tomography for surface and body waves, Geophysics. 80, R281-R302.

Yuan, O. Y. & Simons, F. J., 2014. Multiscale adjoint waveform-difference tomography using wavelets, Geophysics. 79, WA79-WA95.
