## Release summary

This is a new version that is compatible with updates to 'rstantools'.

At the same time, we have made minor updates to avoid R CMD check NOTEs with R devel. 

## Test environments

* local macOS install, R 4.2.2
* Windows (on github-actions), R 4.2.2
* Ubuntu 20.04.4 (on github-actions), R-devel
* Ubuntu 20.04.4 (on github-actions), R-release
* Ubuntu 20.04.4 (on github-actions), R-old-release
* Windows (winbuilder), R devel
* Windows (winbuilder), R release

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking C++ specification ... NOTE
  Specified C++14: please drop specification unless essential
  
This is created by 'rstantools' on build and is presumably essential.
'rstantools' itself has a system requirement of C++14.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
  
This is correct.
