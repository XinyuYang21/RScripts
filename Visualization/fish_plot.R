#install devtools if you don't have it already for easy installation
install.packages("devtools")
library(devtools)
install_github("chrisamiller/fishplot",force=TRUE)

library(fishplot)

#provide a list of timepoints to plot
#You may need to add interpolated points to end up with the desired
#visualization. Example here was actually sampled at days 0 and 150
timepoints=c(0,10,20)      

#provide a matrix with the fraction of each population
#present at each timepoint
frac.table = matrix(
  c(50,20,0,
    98, 40, 55,
    100, 90, 70),
  ncol=length(timepoints))
frac.table = matrix(
  c(100, 45, 00, 00,
    02, 00, 00, 00,
    02, 00, 02, 01,
    98, 00, 95, 40),
  ncol=4)
#provide a vector listing each clone's parent
#(0 indicates no parent)
parents = c(0,1,2)

#create a fish object
fish = createFishObject(frac.table,parents,timepoints=timepoints)

#calculate the layout of the drawing
fish = layoutClones(fish)

#draw the plot, using the splining method (recommended)
#and providing both timepoints to label and a plot title
fishPlot(fish,shape="spline",title.btm="COADREAD",
         cex.title=0.5, vlines=c(10,20), 
         vlab=c("ESig3","ESig4"))

