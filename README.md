# nim
### An R package for fitting non-stationary index-flood models (NIMs)

**The package is still in pre-mature state and has to be used with caution!**

Currently implemented features:
- fitting stationary index-flood models
- fitting non-stationary index-flood models with linear and smooth trends
- bootstrapping
- calculation of the Anderson-Darling statistics
- various methods for working with fitted `nim` and `nims` objects

To be implemented soon:
- evaluation of critical values for the Anderson-Darling test
- multicore/parallel support
- documentation

To be implemented later:
- methods for bandwidth selection

To install the package in R use 

```
devtools::install_github('hanel/nim')
```

For help, examples etc. see ?`nim-package`
