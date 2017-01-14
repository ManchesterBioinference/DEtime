library(DEtime)

test_that("missing arguments handling", {
    data(SimulatedData)
    expect_that(DEtime_infer(ControlTimes), throws_error())
    expect_that(DEtime_infer(ControlTimes, ControlData), throws_error())
    expect_that(DEtime_infer(ControlTimes, ControlData, PerturbedTimes), throws_error())
    })
