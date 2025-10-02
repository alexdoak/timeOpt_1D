A partial Python conversion of the timeOpt algorithm from the astrochron R package

I'm working on a project related to cyclostratigraphy. There is a very elegent method developed by Stephen Meyers called timeOpt, included in the R package astrochron https://cran.r-project.org/web/packages/astrochron/index.html

While exploring new methods, for myself it is useful to be able to reconsctruct timeOpt's output in python, rather than running the code in R. Hence, I have coded the method into python.

Note that the astrochron R package, built by Stephen Meyers, Alberto Malinverno, Linda Hinnov, Christian Zeeden, Huaran Liu, Vincent Moron, and Michel Crucifix, is significantly more developed than this minimalist python version. Please use the python code with caution, and if you are working on a project involving the timeOpt code, you are better off using the R package. At the very least, check outputs agree for some parameters.

For the paper introducing the timeOpt algorithm, please see Meyers, S. R. (2015), The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization, Paleoceanography, 30, 1625–1640, doi:10.1002/2015PA002850

The data used here is taken from the supplementary details of the paper Josso P., et al (2021), Geochemical evidence of Milankovitch cycles in Atlantic Ocean ferromanganese crusts, Earth and Planetary Science Letters, 553, https://doi.org/10.1016/j.epsl.2020.116651

Running "run_timeOpt.py" will load the data used to create figure 5 of Josso et al (2021), perform the timeOpt algorithm, and plot figures 5E-5F from the paper.
