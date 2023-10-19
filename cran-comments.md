This version has minor updates on the request of the 'rstan' team. It removes
deprecated syntax for future 'rstan' compatibility.

At the same time, we have made minor updates to avoid R CMD check NOTEs with
R devel. 

## Test environments

* local macOS M2 install, R 4.3.1
* Windows (on github-actions), R 4.3.1
* Ubuntu 22.04.3 (on github-actions), R-devel
* Windows (winbuilder), R-devel

With memory management checks:
 
* Ubuntu 22.04.3 (on github-actions), R-devel with valgrind

## R CMD check results

0 errors | 0 warnings | 1 notes

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
  
This is correct.
