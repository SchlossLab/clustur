# clustur 0.1.4
* Updated the `cluster()` function to change the current clustering metric it is using if it is not suitable for your data.

* Fixed memory leaks that were found using Valgrind.

# clustur 0.1.3

* Fixed `cluster()` returning duplicate cluster names when using method = "opticlust."
