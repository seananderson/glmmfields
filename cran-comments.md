## Release summary

This is a resubmission that addresses three requested changes:

- Added copyright holder 
  `person("Trustees of", "Columbia University", role = "cph"))` 
  to `Authors@R` field.
  
- Added DOI for method to 'Description' field. This DOI was just assigned by the
  publisher today and will be live when the paper is "posted as an Accepted
  Article online shortly".
  Format: `Anderson and Ward (2018) <doi:10.1002/ecy.2403>`.
  
- Replaced `\dontrun{}` by `\donttest{}` in examples.

## Test environments

* local OS X install, R 3.5.0
* Ubuntu 14.04 (on travis-ci), devel
* win-builder (devel)
* Windows (on AppVeyor), 3.5.0 Patched

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
