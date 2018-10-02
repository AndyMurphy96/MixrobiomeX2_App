# MixrobiomeX2_App

The app currently allows simple microbiome analysis starting from a phyloseq object.


## To Improve

- general note for coding on the app: each tab can be minimized both in the server function and in the ui. Do this first, it makes working on it easier.


- give more user options, for example on the test applied for example in Phylum analysis
- allow do download plots
- add a Phewas tab
- if you want "real" relative abundances, this can easily changed in tab 3 by uncommenting the line: SFs_RA <- sample_sums(rv$ps)
    - however, both the alpha diversity code and DESeq have problems with non integer numbers therefore I used the count way currently
