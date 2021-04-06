Informative Spatial Selection
=============================

`InformativeSpatialSelection` is an R package that contains the source
code to reproduce the graphs and simulations of the article
"InformativeSpatialSelection", by Daniel Bonnéry, Francesco Pantalone
and Giovanna Ranalli.

To install the package:
-----------------------

You will need to have Rtools installed, as well as the devtools package.

    devtools::install_github("DanielBonnery/InformativeSpatialSelection",dependencies=TRUE)

Produce the paper graphs and tables.
------------------------------------

### Load the package

    library(InformativeSpatialSelection)

### Figure 1.

    demo(figure1)

![Test Image 1](latex/fig/figure1.png)

### Figure 2.

    library(InformativeSpatialSelection)
    demo(figure2)

![Paper image 2](latex/fig/figure2.png)

### Figure 3.

    demo(figure3)

![Paper image 3](latex/fig/figure3.png)

### Figure 4.

    library(InformativeSpatialSelection)
    demo(figure4)

![Paper figure 4](latex/fig/figure4.png)

### Figure 5.

    demo(figure5)

![Paper figure 5](latex/fig/figure4.png)

Details - demo of functions - selection.
========================================
