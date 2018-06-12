## Release summary

This is a resubmission that addresses three requested changes:

- Added copyright holder 
  `person("Trustees of", "Columbia University", role = "cph"))`
  to `Authors@R` field.
  
- Added DOI for method to 'Description' field. This DOI was just assigned by the
  publisher today and will be live when the paper is posted as an Accepted
  Article online shortly.
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

* Possibly mis-spelled words in DESCRIPTION:
  Spatiotemporal (4:9)
  spatiotemporal (11:46)
  
These are correct.

* Found the following (possibly) invalid DOIs:
  DOI: 10.1002/ecy.2403
    From: DESCRIPTION
          inst/CITATION
    Status: Not Found
    Message: 404
    
These will become valid very soon and I am hoping this package can 
be released on CRAN slightly before the paper is published.

Note from the journal editorial office:

"I've received the DOI information from the publisher, and your
upcoming publication in Ecology will be published under the 
following DOI: 10.1002/ecy.2403

The publisher should be posting your Accepted Article online shortly;
please provide the article DOI to CRAN to release your related R package."
