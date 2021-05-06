HLIBcore
========

*HLIBcore* contains basic type definitions and function implementations from
[HLIBpro](https://hlibpro.com) needed by [libHLR](http://libhlr.org), e.g., matrix types, clustering
algorithms and input/output functions. 

The API is a fully compatible subset of the API provided by *HLIBpro*, i.e., you can directly
replace *HLIBcore* with *HLIBpro*. However, all the arithmetic functions for H-matrices and even
basic matrix construction is missing.  For this, please use the functions provided by *libHLR*.

