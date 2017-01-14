## ---- eval=FALSE, message=FALSE------------------------------------------
#  ## read simulated example data
#  library("DEtime")
#  data(SimulatedData)
#  res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, PerturbedTimes = PerturbedTimes, PerturbedData=PerturbedData)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  ## read simulated example data
#  library("DEtime")
#  data(SimulatedData)
#  res <- DEtime_rank(ControlTimes = ControlTimes, ControlData = ControlData, PerturbedTimes = PerturbedTimes, PerturbedData=PerturbedData, savefile=TRUE)

## ---- message=FALSE------------------------------------------------------
library("DEtime")
## read simulated example data
data(SimulatedData)
res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, PerturbedTimes = PerturbedTimes, PerturbedData=PerturbedData)
print_DEtime(res)

## ---- message=FALSE------------------------------------------------------
library("DEtime")
## read simulated example data
data(SimulatedData)
res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, PerturbedTimes=PerturbedTimes, PerturbedData=PerturbedData)
###plot_DEtime(res,BestPerturbPos="mean", plot_gene_ID='1')

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

