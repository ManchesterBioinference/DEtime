library(DEtime)

test_that("missing arguments handling", {
    data(SimulatedData)
    expect_that(DEtime_infer(times), throws_error())
    expect_that(DEtime_infer(times, ControlData), throws_error())
    expect_that(DEtime_infer(times=times,ControlData=NULL, PerturbedData=PerturbedData), throws_error())
    
    })
