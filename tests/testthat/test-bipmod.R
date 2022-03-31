test_that("Algorithm Correctness Simple", {
  data(example_data,package = 'BipartiteModularityMaximization')
  Q_Part=BipartiteModularityMaximization::bipmod(example_data)
  expect_equal(Q_Part$MODULARITY,0.2622772)
})
