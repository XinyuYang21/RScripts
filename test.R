install.packages("httpgd")
install.packages('vscDebugger')

devtools::install_github("ManuelHentschel/VSCode-R-Debugger")
library(vscDebugger)
devtools::install_github("ManuelHentschel/vscDebugger")

library(vscDebugger)
library(languageserver)
.vsc.listen()

{
    "type": "R-Debugger",
    "request": "attach",
    "name": "Attach to R process"
}

httpgd

hgd_browse()
hgd()
library(httpgd)
x = seq(0, 3 * pi, by = 0.1)
plot(x, sin(x), type = "l")

library(ggplot2)
ggplot(mpg, aes(displ, hwy, colour = class)) +
  geom_point()


gene = maf@data$Hugo_Symbol  %>% unique(.)
om = createOncoMatrix(m = maf, g = gene)
numMat = om$numericMatrix
om$oncoMatrix
