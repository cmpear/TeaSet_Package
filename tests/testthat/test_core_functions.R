context("Test TeaSet's Core Functions")

test_that("Basic Functions Work",{
  testData <- matrix(c(1,1,1,
                       2,0,0,
                       0,3,0,
                       0,0,4,
                       5,5,0,
                       0,6,6,
                       7,0,7,
                       8,8,8,
                       -9,-9,-9,
                       0,-10,-10,
                       -11,0,-11,
                       -12,-12,0),ncol=3,nrow=12,byrow=TRUE)
  teaSet <- brew_tea(testData)

  expect_equal(teaSet$myIJK,teaSet$get_ijk(raw=TRUE))
  expect_equal(sum(teaSet$myIJK[1:3,1]),0)
  expect_equal(sum(teaSet$myIJK[1:3,2]),0)
  teaSet <- brew_tea(testData)
  xy<-teaSet$get_xy()
  expect_equal(c(0,0),xy[1,])
  expect_equal(teaSet$myIJK[1,],xy[2,])
  expect_equal(teaSet$myIJK[2,],xy[3,])
  expect_equal(teaSet$myIJK[3,],xy[4,])
  expect_equal((teaSet$myIJK[1,]+teaSet$myIJK[2,])/c(2,2),xy[5,])
  expect_equal((teaSet$myIJK[2,]+teaSet$myIJK[3,])/c(2,2),xy[6,])
  expect_equal((teaSet$myIJK[1,]+teaSet$myIJK[3,])/c(2,2),xy[7,])
  expect_equal(c(0,0),xy[8,])
  expect_equal(c(0,0),xy[9,])
  expect_equal(teaSet$myIJK[1,],xy[10,])
  expect_equal(teaSet$myIJK[2,],xy[11,])
  expect_equal(teaSet$myIJK[3,],xy[12,])

  teaSet$set_frame(stretch=c(-2,-2),inplace=TRUE)
  flip <- teaSet$get_xy()
  expect_equal(xy,flip * c(-1/2,-1/2))
})
