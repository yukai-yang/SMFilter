## Summary

This is a maintenance update to the CRAN package "SMFilter".

Version 1.0.5 fixes Rd markup issues reported by CRAN, updates package citation metadata, and modernises the package-level roxygen documentation while keeping the statistical methodology and the user-facing workflow unchanged.

## Main changes

* Fixed Rd markup for lagged dependent-variable subscripts in `SimModel1()` and `SimModel2()`, resolving the "Lost braces" notes in the Rd files check.
* Updated the package citation metadata to refer to version 1.0.5 while retaining the original software year.
* Updated package-level roxygen documentation to use the recommended `"_PACKAGE"` sentinel instead of deprecated `@docType package`.
* Updated README and NEWS entries for the maintenance release.
* Excluded repository-only files from the source package build.

Existing user code following documented usage should continue to work as before.

## Reverse dependencies

There are no reverse depends, imports, or linking-to dependencies on CRAN.

There is one reverse suggests dependency on CRAN: `tip`.

## R CMD check results

Platform: aarch64-apple-darwin23  
R version: 4.6.1  

0 errors | 0 warnings | 0 notes
