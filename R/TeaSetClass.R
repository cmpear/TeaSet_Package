#' @title TeaSet class
#' @description This class and accompanying functions are designed to help with the creation of ternary graphs and the manipulatoin of ternary data.  The class uses snake_names for functions and methods, camel-case for variables, and Pascal Case for everything else.  For graphical methods and functions, I have tried to emulate similar functions' parameter names.
#' @exportClass TeaSet
#' @importFrom methods new
#' @importFrom graphics polygon
#' @importFrom stats sd
#' @field myCoord a numeric matrix
#' @field myCenter a numeric list of xy-coordinates
#' @field myStretch a numeric list of length two
#' @field myIJK a numeric matrix, dimensions 4X2, i is first and last
#' @note allows for blank class to be built
#' @examples
#'   \dontrun{testData <- matrix(c(rnorm(300)),ncol=3)
#'            teaSet   <- TeaSet$new(testData)}
#'

TeaSet <- setRefClass("TeaSet",
                     fields = list(myCoord   = "matrix",
                                   myCenter  = "numeric",    # x_val, y_val
                                   myStretch = "numeric",    # x_val, y_val
                                   myIJK     = "matrix"      # has a second i-value at the end to make plotting easier
                     ),   # Close fields, comma to methods
                     methods = list(
                       ## INITIALIZE ##
                       initialize = function(data_ = NA, columns_ = c(1,2,3), center_ = NA, xFrame_ =NA, yFrame_ = NA)
                       {
                         "Initializes a TeaSet object, used when TeaSet$new() is called.  No arguments necessary."
                         myStretch <<- p_size_frame(xFrame_, yFrame_)
                         myCenter  <<- p_center_frame(xFrame_, yFrame_)
                         myIJK     <<- get_ijk()
                         if ((!is.data.frame(data_) && !is.matrix(data_)))
                         {
                           myCoord <<- get_ijk() + matrix(rep(myCenter,4), ncol=2, byrow=TRUE)
                         }
                         else{
                           myCoord <<- as.matrix(data_[,columns_])
                         }
                       }, # close initialize
                       set_frame = function(xFrame = NA, yFrame = NA, center = c(NA,NA), stretch = c(NA,NA),inplace = FALSE){
                         rValue <- matrix(NA, nrow = 6, ncol = 2)
                         if (!any(is.na(stretch))){
                           if(inplace){
                             myStretch <<- stretch
                           }
                           rValue[2,] <- stretch
                         }
                         else{
                           rValue[  2,]<-p_size_frame(xFrame,yFrame,inplace = inplace)
                         }
                         if (!any(is.na(center)))
                         {
                           if(inplace){
                             myCenter<<-center
                           }
                           rValue[1,]<-center
                         }
                         else{
                           rValue[1,]<-p_center_frame(xFrame, yFrame, inplace = inplace)
                         }
                         rValue[3:6,]<- get_ijk(inplace = inplace)
                         return(rValue)
                       }, # close set_frame
                       get_ijk = function(raw = FALSE, inplace = FALSE)
                       {
                         "Returns myIJK. If raw=TRUE, returns default ijk-vectors. If inplace=TRUE, recalculates myIJK or sets myIJK to default values."
                         if (raw){
                           rValue <- matrix(data = c( 0  , sqrt(3)/3,
                                                      -0.5,-sqrt(3)/6,
                                                      0.5,-sqrt(3)/6,
                                                      0  , sqrt(3)/3), ncol = 2, nrow = 4, byrow = TRUE)
                         }
                         else{
                           rValue <- matrix(data = c( 0   * myStretch[1], sqrt(3)/3 * myStretch[2],
                                                      -0.5 * myStretch[1],-sqrt(3)/6 * myStretch[2],
                                                      0.5 * myStretch[1],-sqrt(3)/6 * myStretch[2],
                                                      0   * myStretch[1], sqrt(3)/3 * myStretch[2]), ncol = 2, nrow = 4, byrow = TRUE)
                         }
                         if(inplace){
                           myIJK <<- rValue
                         }
                         return(rValue)
                       }, # close get_ijk
                       get_xy = function(data=myCoord){
                         "Given t-coordinates, returns xy-coordinates. Defaults to myCoord."
                         # scale is used to ensure that the TeaPoint fits a given plot--by default a 1*1 centered on the origin

                         if (is.matrix(data) || is.data.frame(data)){
                           return(t(apply(data,MARGIN=1, FUN = get_xy)))
                         }
                         # in ternary space, negative values instead mean moving in the other two directions simultaneously
                         # at this point, data should just be a list of 3 numbers

                         data <- p_normalize_ternary_pt(data)  # apply transposes data, so must un-transpose it
                         rValue <- c(0,0)
                         #                names(rValue) <- c('Xs','Ys')
                         rValue[1] <- (data[1] * myIJK[1,1] + data[2] * myIJK[2,1] + data[3] * myIJK[3,1]) + myCenter[1]
                         rValue[2] <- (data[1] * myIJK[1,2] + data[2] * myIJK[2,2] + data[3] * myIJK[3,2]) + myCenter[2]
                         return(rValue)
                       }, # close get_xy
                       get_leaf = function(data){
                         "Given t-coordinates, returns what section of the t-plot they are in."
                         # takes a ternary point or raw list of ternary points
                         if (length(data)>3){
                           return(matrix_apply(data,nrow=3,fun=get_leaf))
                         }
                         x <- data[1]; y <- data[2]; z <- data[3]

                         # exact ratio does not matter if any are negative
                         # but negatives will still screw with function--setting to zero
                         if (x <0) x<-0
                         if (y<0) y<-0
                         if (z<0) z<-0
                         rValue <- 0L  # will be tested at ### rValue01
                         if      (x == y && x == z) return( 0L)  # point at center
                         if      (x == y && z<x)    {rValue <- 4L} # on i-leaf/j-leaf seperator
                         else if (z == x && y<z)    {rValue <- 5L} # on i-leaf/k-leaf seperator
                         else if (y == z && x<y)    {rValue <- 6L} # on j-leaf/k-leaf seperator
                         if (rValue %in% c(4,5,6)){       ### rValue01
                           u <- max(data); v <- min(data)
                         }
                         else {
                           if (max(data)==x){
                             rValue <- 1L  # i-leaf
                             u <- x
                             v <- y + z
                           }
                           if (max(data) ==y){
                             rValue <- 2L  # j-leaf
                             u <- y
                             v <- x + z
                           }
                           if (max(data) ==z){
                             rValue <- 3L  # k-leaf
                             u <- z
                             v <- x + y
                           }
                         }
                         if      (any(data==0))  {rValue <- rValue + 40L}  # outer edge
                         else if (u <  v)       {rValue <- rValue + 10L}  # inner triangle
                         else if (u == v)       {rValue <- rValue + 20L}  # inner-outer border
                         else if (u >  v)       {rValue <- rValue + 30L}  # outer triangle
                         else{
                           stop("ERROR: u should be either greater than, equal to, or lesser than v")
                         }
                         return(rValue)
                         # 0 Origin
                         # 1s for triangle
                         #  1-Top, 2-Left, 3-Right
                         #  4-Top-left Border, 5-Top-right Border, 6-Left-Right Border
                         # 10s for inner, inner-outer border, outer, edge
                         # 10 - Inner Triangle
                         # 20 - Inner-Outer Border
                         # 30 - Outer Triangle
                         # 40 - Edge
                       }, # close get_leaf
                       tea_segments = function(p0, p1=NA, draw = TRUE,...){
                         # each value must be from 0 to 1, or 1 to 100
                         if (is.matrix(p0) | is.data.frame(p0)){
                           if (ncol(p0) == 3 && !any(is.na(p1)) && ncol(p1) == 3){
                             p0 <- cbind(p0, p1) # easier to work with one list of 6 elements
                           }
                           if (ncol(p0)==6){
                             return(apply(p0, MARGIN = 1, FUN = tea_segments,draw=draw,...))
                           }

                           else{
                             message("WARNING: tea_segments() requires p0 and p1 to be t-points, 3-length numeric lists or 3-column matrices; p0 can also be length-6, representing two tea-points")
                           }
                         }
                         else{
                           if (length(p0) ==3 && !any(is.na(p1)) && length(p1)==3){
                             xy <- get_xy(matrix(c(p0,p1), ncol = 3, byrow=TRUE))
                           }
                           if (length(p0) ==6){
                             xy <- get_xy(matrix(p0, ncol = 3, byrow=TRUE))
                           }
                           while(max(p0)>1) p0<-p0/10
                           xy <- unlist(list(xy))  # plain English aside, this actually makes sense
                           xy <- xy[c(1,3,2,4)]  # we're filling by rows, but unlist treats matricies as though they were filled by columns
                           if(draw)
                           {
                             p_segments(x0 = xy[1], y0 = xy[2], x1 = xy[3], y1 = xy[4],...)
                           }
                           return(xy)
                         }
                         stop("ERROR: Invalid arguments for tea_segments--must be a vector of length 3 or 6, or a matrix of 3 or 6 columns")
                       }, # close tea_segments
                       tea_lines = function(x = NULL, y = NULL, z = NULL, draw = TRUE,overDraw = c(FALSE,FALSE),col = 'grey',...){
                         "Draws a line given a single value along one of the ternary axes.  Points should be from 0-1. Works wtih lists of points. Overdraw makes the line go slightly over the border triangle."
                         # a tea-line only needs one of x, y and z
                         # x can either be a numeric, or a list of numerics,
                         # numeric values in the range 0-1 are tied to their variable (x, y or z)
                         # there is also a coding system for the x variable to contain x, y and z
                         ## range [10,11] to indicate x
                         ## range [20,21] to indicate y
                         ## range [30,31] to indicate z
                         # draws nothing if draw=FALSE
                         # returns xy coordinates to make a line-segment
                         if (length(x)>1 && max(x)>1){
                           # handles a single list using just x
                           return(sapply(x, MARGIN = 1, FUN = tea_lines, draw= draw,overDraw = overDraw,col=col,...))
                         }
                         # handles 2-3 lists using x, y or z
                         if (2<=(sum(length (x), length(y), length(z)))){
                           if (!is.null(x)) if (max(x)<=1) x<- x+10
                           if (!is.null(y)) if (max(y)<=1) y<- y+20  # odd structure to avoid warnings
                           if (!is.null(z)) if (max(z)<=1) z<- z+30
                           return(sapply(c(x,y,z), FUN = tea_lines, draw= draw,overDraw = overDraw,col=col,...))
                         }
                         # should only be one x XOR one y XOR one z
                         # feels easier to just feed everything through x
                         if (!is.null(y)){
                           if (y <=1) y<- y+20
                           x <-y
                         }
                         if (!is.null(z)){
                           if (z <=1) z<- z+30
                           x<-z
                         }
                         # no lists should get beyond this point
                         if (x<=1){
                           tPoints<-c(x,1-x,0 ,
                                      x,0  ,1-x)
                         }
                         else if (x>=10 && x<=11){
                           x2<- x-10
                           tPoints<-c(x2,1-x2,0 ,
                                      x2,0   ,1-x2)
                         }
                         else if (x>=20 && x<=21){
                           x2 <- x -20
                           tPoints<-c(1-x2,x2,0 ,
                                      0   ,x2,1-x2)
                         }
                         else if (x>=30 && x<=31){
                           x2 <- x - 30
                           tPoints<-c(1-x2,0   ,x2,
                                      0   ,1-x2,x2)
                         }
                         xy <- tea_segments(tPoints, draw = FALSE,...)

                         if(any(overDraw)){
                           if (x<=1 | (x>=10 && x<=11)){
                             vec <-(myIJK[2,] - myIJK[3,])/40
                             if (overDraw[1]){
                               xy[1:2] <- xy[1:2] + vec
                             }
                             if (overDraw[2]){
                               xy[3:4] <- xy[3:4] - vec
                             }
                           }
                           else if (x >=20 && x<=21){
                             vec <-(myIJK[3,] - myIJK[1,])/40
                             if (overDraw[1]){
                               xy[3:4] <- xy[3:4] + vec
                             }
                             if (overDraw[2]){
                               xy[1:2] <- xy[1:2] - vec
                             }
                           }
                           else if (x >=30 && x<=31){
                             vec <-(myIJK[1,] - myIJK[2,])/40
                             if (overDraw[1]){
                               xy[1:2] <- xy[1:2] + vec
                             }
                             if (overDraw[2]){
                               xy[3:4] <- xy[3:4] - vec
                             }
                           }
                         }
                         if (draw){
                           segments(x0 = xy[1], y0 = xy[2], x1 = xy[3], y1 = xy[4],col=col,...)
                         }
                         return(xy)
                       }, # close tea_lines
                       tea_triangle = function(r=1,xy = myCenter,jCorner = FALSE,raw = FALSE,draw = TRUE, recursing = FALSE,
                                               col = "gradient",alpha = 0.3, ...){
                         "Creates triangle similar to the triangle around the t-plot. r is its size relative to the border-triangle, and xy is its center.  If jCorner=TRUE, instead draws from the j-vector's corner."
                         # in an equilateral triangle, |r| is the length of the sides
                         # negative r flips the triangle vertically and horizontally
                         # jCorner is provided in case you want to draw the triangle from the tip of vector-j
                         # the function still works from the center, so the first step is to move back to the center
                         if(jCorner){
                           center <- center - myIJK[2,]
                         }
                         if (recursing){
                           xy <- r[2:3]
                           r <- r[1]
                         }
                         if (length(r)>1){
                           return(apply(cbind(r,xy),MARGIN=1,FUN=tea_triangle, raw = raw, draw = draw, recursing = TRUE,...))
                         }
                         ijk <- get_ijk(raw = raw)
                         ijk <- ijk * r + matrix(c(xy,xy,xy,xy), nrow = 4, byrow=TRUE)
                         if (draw){
                           if (!is.na(col) && col == "gradient"){
                             col<- tea_xy_gradient_colorize(xy, alpha = alpha)
                           }
                           p_polygon(x = ijk[1:3,1],y=ijk[1:3,2],col = col,...)
                         }
                         return(ijk)
                       }, # close tea_triangle
                       tea_gradient_colorize = function(x=myCoord, alpha = .30){
                         "Returns a list of colors for t-points based on their locations within t-plot"
                         if (is.matrix(x) || is.data.frame(x)){
                           return(apply(x,MARGIN=1,FUN=tea_gradient_colorize,alpha=alpha))
                         }
                         x<-p_normalize_ternary_pt(x)
                         return(rgb(x[1],x[2],x[3], alpha = alpha, max = 1))
                       }, # close tea_graident_colorize
                       tea_xy_gradient_colorize = function(x, alpha = .30, rotate_hue = 0, raw=FALSE){
                         "Returns a list of colors for xy-points based on their locations within t-plot"
                         if (is.matrix(x) || is.data.frame(x)){
                           return(apply(x,MARGIN=1,FUN=tea_xy_gradient_colorize,alpha=alpha))
                         }
                         if (!raw){
                           x<-x - myCenter
                         }
                         x <- (180 / pi * atan2(x[2],x[1]) - 90 + rotate_hue) %% 360
                         return(hcl(x,c=100, l=60, alpha = alpha))
                       }, # close tea_xy_gradient_colorize
                       tea_contrast_colorize = function(x=myCoord, alpha = 0.3){
                         "Returns a list of colors based on sections of provided t-points."
                         if (is.matrix(x) || is.data.frame(x)){
                           return(apply(x, MARGIN=1, FUN = tea_contrast_colorize, alpha = alpha))
                         }
                         else if ((length(x)%%3)!=0){
                           stop("ERROR: tea_contrast_colorize provided wtih incomplete ternary data")
                         }
                         else if (length(x)>3){  ## raw ternary points
                           return(matrix_apply(x,nrow=3, fun=tea_contrast_colorize, alpha = alpha))
                         }
                         x<-p_redistribute_negatives(x)
                         section <- as.character(get_leaf(x)%% 10) # modding by 10 removes complexity, giving just the leaf
                         c<-1
                         key<-c(    rgb(0,0,0,alpha = alpha, max = 1),
                                    rgb(c,0,0,alpha = alpha, max = 1),
                                    rgb(0,c,0,alpha = alpha, max = 1),
                                    rgb(0,0,c,alpha = alpha, max = 1),
                                    rgb(c,c,0,alpha = alpha, max = 1),
                                    rgb(0,c,c,alpha = alpha, max = 1),
                                    rgb(c,0,c,alpha = alpha, max = 1))
                         key <- p_section_key(key)
                         return(key[section])
                       }, # tea_contrast_colorize
                       tea_gradient_background = function(rows = 30, border = NA, alpha=0.1,...){
                         "Creates a gradient of colors in the background of a t-plot. The more rows (of triangles), the finer the gradient."
                         width <- myIJK[3,1] - myIJK[2,1]
                         widthInc <- width / rows
                         triHeight <- 1 /rows

                         heightIncThird <- myIJK[1,2]/rows / 2

                         start <- myIJK[2,] + c(widthInc/2,heightIncThird) + myCenter
                         upDown <- 1;
                         for (c in 1:rows){
                           for (w in 1:(rows - c + 1))
                           {
                             coord <- start + c((w-1) * widthInc,0)
                             tea_triangle(r = triHeight, xy=coord, border = border, alpha=alpha,...)
                           }
                           start <- start + c(widthInc/2,heightIncThird)    # one third to get from upright to inverted center
                           w <- 1
                           while (w <=(rows - c))  ## on last round of parent loop, while loop shouldn't execute
                           {
                             coord <- start + c((w-1) * widthInc,0)
                             tea_triangle(r = 0-triHeight, xy=coord, border = border, alpha=alpha,...)
                             w <- w +1
                           }
                           start <- start + c(0,heightIncThird*2)  # one third to get to the edge, another to the new center
                         }
                         return(0)
                       }, # close tea_gradient_background
                       tea_ticks = function(axes=4, col.axis = "grey"){
                         "Draws tick-lines for all three axes and labels them."
                         if (is.logical(axes)) axes<-4
                         if(length(axes)==1){
                           axes <- 1:(axes-1)/axes
                         }
                         while(max(axes)>1){
                           axes<-axes/10
                         }
                         tea_lines(x=axes,y=axes,z=axes,draw=TRUE,col=col.axis, overDraw = c(TRUE,FALSE))
                         axesLabs <- as.character(round(axes * 100,2))
                         axesLabs <-sapply(axesLabs, paste0, "%")

                         text_line<-function(start, end, axes, axesLabs){
                           upright<-(myIJK[1,2]>0)
                           if(start[2] == end[2])  # bottom
                           {
                             if (upright) pos<-1
                             else         pos<-3
                           }
                           # testing min of start and end versus both j and k, takes care of both normal and flipped varients
                           else if (min(start[1],end[1])>min(myIJK[2:3,1])){
                             pos <- 4  # right
                           }
                           else{
                             pos <- 2  # left
                           }
                           vec <- start - end
                           here <- matrix(rep(start, length(axes)), ncol=2, byrow=TRUE) - matrix(rep(vec, length(axes)), ncol=2, byrow=TRUE) * matrix(rep(axes, each = 2), ncol=2, byrow=TRUE)
                           text(here, labels = axesLabs, pos = pos)
                         }
                         text_line(start = myIJK[2,], end = myIJK[1,], axes = axes, axesLabs = axesLabs)
                         text_line(start = myIJK[3,], end = myIJK[2,], axes = axes, axesLabs = axesLabs)
                         text_line(start = myIJK[4,], end = myIJK[3,], axes = axes, axesLabs = axesLabs)
                       }, # close tea_ticks
                       tea_label_axes = function(axis.labels = "default"){
                         "Labels the axes of a t-plot"
                         if (all(axis.labels=="default")){
                           axis.labels <- colnames(myCoord)
                         }
                         upright<- myStretch[1]>0
                         # pos: 1=below, 2=left, 3=above,4=right
                         if(upright) {jPos <- 1; kPos <- 1} else {jPos <- 3; kPos <- 3}
                         iPos <- 2
                         ijk <- myIJK[1:3,] + matrix(rep(myCenter,3), ncol=2, byrow=TRUE)
                         text(ijk, labels = axis.labels, pos = c(iPos,jPos,kPos))
                       }, # close tea_label_axes
                       tea_frame_plot = function(...){
                         "Draws the border triangle around a t-plot.  May also provide a color to shade the background"
                         bg<-extractFrom...("bg",defaultTo=NA,...)
                         frame.plot<-extractFrom...("frame.plot",defaultTo=TRUE,...)
                         if (frame.plot) border<-"black"
                         else border<-NA
                         if(!is.na(bg) && bg=="gradient"){
                           tea_triangle(r=1,col=NA,border=border,...)
                           tea_gradient_background()
                         }
                         else{
                           tea_triangle(r=1,col=bg,border=border,...)
                         }
                       }, #  close tea_frame_plot
                       tea_plot = function(data = myCoord, dataLabels=NA,  main = "", newPlot = TRUE, bullseye = FALSE,
                                           ticks = NA, axis.labels = "default", col.axis="grey",...){
                         "Creates a t-plot.  If newPlot=FALSE, will instead draw atop current plot."
                         if(newPlot){
                           p_plot(myIJK* 1.15, frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "", type = 'n',...)
                           mtext(text = main, cex = 1.5, side=3)
                         }
                         else{
                           if(myIJK[1,2]>0){
                              here<-myCenter + myIJK[1,]
                              text(here[1],here[2],labels=main,cex=1.5,pos=3)
                           }
                           else{
                             here <- myCenter + c(0,myIJK[2,2])
                             text(here[1],here[2],labels=main,cex=1.5,pos=3)
                           }
                         }
                         tea_frame_plot(...)
                         # uses integer or list of integers for ticks
                         if(!all(is.na(ticks)) && !(is.logical(ticks) && ticks==FALSE)){
                           tea_ticks(ticks, col.axis = col.axis)
                         }
                         # axes?
                         if (!any(is.na(axis.labels)) &&(is.character(axis.labels) || all(is.logical(axis.labels),axis.labels==TRUE))){
                           if (all(axis.labels=="default") || all(is.logical(axis.labels))){
                             tea_label_axes(colnames(data))  # just in case data was inputted rather than using myCoord
                           }
                           else{
                             tea_label_axes(axis.labels)
                           }
                         }
                         if(bullseye){
                           tea_segments(c(1,1,0,1,1,1),col=col.axis)
                           tea_segments(c(0,1,1,1,1,1),col=col.axis)
                           tea_segments(c(1,0,1,1,1,1),col=col.axis)
                         }

                         ## define xy ####
                         xy<-get_xy(data)
                         col <- p_color_center(xy,...)
                         if(all(is.na(dataLabels))){
                           p_points(xy, col = col,...)
                         }
                         else{
                           if (is.logical(dataLabels) && all(dataLabels==TRUE)){
                             dataLabels<-rownames(myCoord)
                           }
                           p_text(xy,labels=dataLabels,col=col,...)
                         }
                       },
                       ## Private Fuctions ####
                       p_plot     = function(...,alpha){
                         "private function, removes variables after ... before calling its namesake"
                         plot(...)},
                       p_segments = function(...,alpha){
                         "private function, removes variables after ... before calling its namesake"
                         segments(...)},
                       p_points   = function(...,alpha){
                         "private function, removes variables after ... before calling its namesake"
                         points(...)},
                       p_polygon  = function(...,alpha){
                         "private function, removes variables after ... before calling its namesake"
                         polygon(...)},
                       p_text     = function(...,alpha){
                         "private function, removes variables after ... before calling its namesake"
                         text(...)},
                       # END WRAPPER PRIVATE FUNCTIONS
                       p_color_center = function(xy,...){
                         "Private method. Works with tea_plot to extract col variable from ..., and handle gradient and contrast colors."
                         col <-extractFrom...("col",defaultTo="gradient",...)
                         alpha <- extractFrom...("alpha",defaultTo=0.3,...)
                         if (length(col)==1){
                           if (!is.na(col) && substr(col,1,8) == 'gradient'){
                             if(nchar(col)>8){
                               rotate_hue <-as.integer(substr(col,9,11))
                             }
                             else rotate_hue<-0
                             col<- tea_xy_gradient_colorize(xy, alpha=alpha, rotate_hue=rotate_hue)
                           }
                           else if (!is.na(col) && col == 'contrast'){
                             col<- tea_contrast_colorize(alpha=alpha)
                           }
                         }
                         return(col)
                       }, # close p_color_center
                       p_center_frame = function(xFrame, yFrame, inplace = FALSE){
                         "Private method. Returns center of xFrame and yFrame.  If inplace=TRUE, sets myCenter<<-center"
                         center <- c(0,0)
                         if (!any(is.na(xFrame) | length(xFrame)<2)) {
                           center[1] <- mean(xFrame)
                         }
                         if (!any(is.na(yFrame) | length(yFrame)<2)){
                           x <-(yFrame[1]+yFrame[1]+yFrame[2])/3  # mean function was returning integers on test...divide more reliable
                           center[2] <- min(yFrame) + (x - min(yFrame)) * sqrt(3)/2
                         } # for frame 0-1, triangle will be sqrt(3)/2 tall. Center of triangle determined by weighted average
                         if (inplace) {
                           myCenter <<- center
                         }
                         return(center)
                       }, # close p_center_frame
                       p_size_frame = function(xFrame = NA, yFrame = NA, inplace = FALSE){
                         "Private method. Returns the difference between xy-values of xFrame and yFrame.  If inplace=TRUE, sets myStretch<<-newStretch"
                         newStretch <-c(NA,NA)
                         if (any(is.na(xFrame)))     { newStretch[1] <- 1 }
                         else if (length(xFrame)>1)  { newStretch[1] <- xFrame[2] - xFrame[1] }
                         else                        { newStretch[1] <- xFrame }
                         if (any(is.na(yFrame)))     { newStretch[2] <- 1 }
                         else if (length(yFrame)>1)  { newStretch[2] <- yFrame[2] - yFrame[1] }
                         else                        { newStretch[2] <- yFrame }

                         if (inplace) { myStretch <<- newStretch}
                         return(newStretch)
                       }, # close p_size_frame
                       p_redistribute_negatives = function(x=myCoord){
                         "Private method. Given ternary-coordinates, redistributes negative values for each axis to the other two axes as half their value."
                         if (is.matrix(x) || is.data.frame(x)){
                           return(apply(x,MARGIN = 1, FUN = p_redistribute_negatives))
                         }
                         these <- x<0
                         add   <- c(0,0,0)
                         if (these[1]){
                           add[2] <- add[2] - x[1]/2   # x is -1, so this works like addition
                           add[3] <- add[3] - x[1]/2
                         }
                         if (these[2]){
                           add[1] <- add[1] - x[2]/2
                           add[3] <- add[3] - x[2]/2
                         }
                         if (these[3]){
                           add[1] <- add[1] - x[3]/2
                           add[2] <- add[2] - x[3]/2
                         }
                         x[these]<- 0
                         x <- x + add
                         return(x)
                       }, # close p_redistribute_negatives
                       p_normalize_ternary_pt = function(x = myCoord){
                         "Private method. For each ternary point, calls p_redistribute negatives, than divides all by their max value."
                         if (is.matrix(x) || is.data.frame(x)){
                           return(apply(x,MARGIN = 1, FUN = p_normalize_ternary_pt))
                         }
                         if (length(x)>3){
                           if (length(x) %% 3 ==0)
                           {
                             x<-matrix(x, nrow=3, byrow = TRUE)
                             return(apply(x, MARGIN = 1, FUN = p_normalize_ternary_pt))
                           }
                           else{
                             stop("ERROR: p_normalize_ternary_pt requires x to be a list of length 3 or multiple of 3")
                           }
                         }
                         x <- p_redistribute_negatives(x)
                         x <- x/sum(x)
                         return(x)
                       }, # close p_normalize_ternary_pt
                       p_section_key = function(x = NULL){
                         "Private method.  Given a list, returns a list with names set to section names"
                         if(any(is.null(x))){
                           x  <- (c("0" = "Origin",
                                    "1" = "i-leaf",
                                    "2" = "j-leaf",
                                    "3" = "k-leaf",
                                    "4" = "i-j-border",
                                    "5" = "i-k-border",
                                    "6" = "j-k-border",
                                    "11"= "inner-i-leaf",
                                    "12"= "inner-j-leaf",
                                    "13"= "inner-k-leaf",
                                    "14"= "inner-i-j-border",
                                    "15"= "inner-i-k-border",
                                    "16"= "inner-j-k-border",
                                    "21"= "inner-border-i-leaf",
                                    "22"= "inner-border-j-leaf",
                                    "23"= "inner-border-k-leaf",
                                    "31"= "outer-triangle-i-leaf",
                                    "32"= "outer-triangle-j-leaf",
                                    "33"= "outer-triangle-k-leaf",
                                    "41"= "outer-border-i-leaf",
                                    "42"= "outer-border-j-leaf",
                                    "43"= "outer-border-k-leaf") )
                         }
                         # create a new key
                         else{
                           labs <- as.character (c(0:6,11:16,21:23,31:33,41:43))
                           names(x) <- labs[1:length(x)]
                         }
                         return(x)
                       } # close p_section_key
                     )) # close methods and fields



## TEASET UTILITY FUNCTIONS ####
#' matrix_apply
#' @description applies the matrix function with nrow to a list before using apply with a function
#' @param x a 1d list of data
#' @param nrow how many rows to put the list into
#' @param fun the function to apply() to the matrix after creating it
#' @param ... Any additional arguments for the function provided
#' @return a list from using the apply function on the new matrix
#' @note length of new list must be evenly divisible by nrow.  This function was built to simplify some code in TeaSet
#' @examples test_data <- 1:16
#'          matrix_apply(test_data, nrow = 4, fun = sum)
#' @export

# for some reason. packages give me errors with <- for functions...
matrix_apply <- function(x, nrow, fun,...){
  x<-matrix(x, nrow = nrow, byrow=TRUE)
  if(!((length(x) %% nrow) == 0)){
    stop("ERROR: in matrix_apply, length of list x not evenly divisible by nrow")
  }
  return(apply(x, MARGIN = 1, FUN = fun,...))
}
#' extractFrom...
#' @description Given a variable name and ..., attempts to extract the variable from ...
#' @param varNames a string with the name of a variable, or a list of such strings
#' @param defaultTo if variable not in ..., return this.  Not designed to work with lists.
#' @param ... Parameters from which to extract variables
#' @return the variable's value, or a list of the variables' values.  NULL if cannot be found
#' @examples \dontrun{extractFrom("col",...)}
#'
#' @export
extractFrom...<-function(varNames,defaultTo=NULL,...){
  if(length(varNames)>1){
    return(sapply(varNames,FUN=extractFrom...,defaultTo,...))
  }
  x<-list(...)[varNames][[1]]
  if(is.null(x)) return(defaultTo)
  return(x)
}
#' normalize
#' @description normalizes or standardizes the columns of a dataframe using one of three methods: max, min-max, or zscore
#' @param data a 2d array, matrix, or data.frame of numeric data.  Defaults to min-max
#' @param type a string indicating how to standardize: max, min-max, or z-score
#' @return normalized data
#' @examples data(mtcars)
#'           normalize(mtcars, type = "z-score")
#' @export
normalize<-function(data, type = "min-max"){
  if (type %in% c("max","Max","MAX","simple")){
    fun <- function(x){
      m <- max(abs(x))
      return(x/m)
    }
  }
  else if(type %in% c("z-score","Z-Score","zscore","zScore","ZScore","z score","Z Score")){
    fun <- function(x){
      ave <- mean(x)
      sDev <- sd(x)
      return((x-ave)/sDev)
    }
  }
  else if (type %in% c("min-max","min max","minmax","Min-Max","Min Max","MinMax","MIN-MAX","MIN MAX","MINMAX","minMax")){
    fun <- function(x){
      MAX<- max(x)
      MIN<- min(x)
      return((x-MIN)/(MAX-MIN))
    }
  }
  else{
    message("WARNING: should input \"z-score\", \"max\", or \"min-max\" or similar inputs, defaulting to \"min-max\". ")
    fun <-function(x){
      MAX<- max(x)
      MIN<- min(x)
      return((x-MIN)/(MAX-MIN))
    }
  }
  return(apply(data,MARGIN=2,FUN=fun))
}
## TEASET ACCESS FUNCTIONS ####
#' @title brew_tea
#' @description returns a TeaSet object
#' @param data a 3-column dataframe, defaults to NA
#' @param ... Additional parameters that work with TeaSet$new()
#' @return a TeaSet object
#' @note an easier way to create a TeaSet object
#' @example \dontrun{brew_tea(matrix(c(rnorm(300)),ncol=3))}
#' @export
brew_tea <- function(data=NA,...){
  teaSet<-TeaSet$new(data,...)
  return(teaSet)
}
#' tea_plot
#' @description a plot function specifically for TeaSets.  Calls teaSet$tea_plot(parameters)
#' @param teaSet a teaSet object. No default.
#' @param dataLabels labels for t-points.  If included, calls text instead of points.
#' @param main the title of the plot
#' @param newPlot whether to start a new plot, or draw atop an old plot. Defaults to TRUE
#' @param bullseye whether to draw a Y-shape separating i-leaf, j-leaf and k-leaf
#' @param ticks Whether to draw tick-lines.  Can either specify the lines (list of numerics between 0 and 1), the number of lines (integer), or whether to draw for quarters (TRUE)
#' @param axis.labels a 3-length list with labels for three axes.  If left as axis.labels="default", will use colnames(myCoord)
#' @param col.axis the color for axis lines and the bullseye.  Defaults to grey
#' @param ... normal graphical parameters not already listed: col, pch, cex, etc.  Can also set col="contrast", "gradient", to rotate gradient, add degrees such as "gradient90"
#' @return NONE
#' @note same functionality as teaSet$tea_plot()
#' @example \dontrun{testData<-testData <- as.data.frame(matrix(rnorm(300, sd=10),ncol=3))
#'                   teaSet<-brew_tea(testData)
#'                   plot(teaSet, main = "My Ternary Plot",col = "gradient180)}
#' @export
tea_plot<-function(teaSet,dataLabels = NA, main ="",newPlot = TRUE, bullseye = FALSE,
            ticks = NA, axis.labels="default",col.axis="grey",...){
  teaSet$tea_plot(dataLabels=dataLabels, main=main,newPlot=newPlot,bullseye=bullseye,
                  ticks=ticks,axis.labels=axis.labels,col.axis=col.axis,...)
}

#' @title ternary_to_xy
#' @description converts ternary-coordinates to xy-coordinates
#' @param data data consists of a matrix or data.frame of ternary points, or a list with the coordinates for one ternary point
#' @param ... Additional parameters that work with TeaSet$new()
#' @return a matrix of xy-coordinates
#' @note does not handle strings, NAs, or NULLs in data
#' @example \dontrun{ternary_to_xy(matrix(c(rnorm(300)),ncol=3))}
#' @export
ternary_to_xy <- function(data,...)
{
  teaSet<-TeaSet$new(data,...)
  return(teaSet$get_xy())
}

#' @title quick_triangle
#' @description creates an equilateral triangle with all sides equal to r centered on a given point
#' @param r the side length
#' @param center the center of the triangle, as an xy-coordinate
#' @param draw a boolean indicating whether or not to actually draw the triangle
#' @param ... Additional parameters that work with polygon
#' @return the coordinates of the triangle's vertices, with the first coordinate relisted at the end
#' @note not as flexible as tea_triangle() in the TeaSet class
#' @example \dontrun{quick_triangle()}
#' @export
quick_triangle <- function(r=1,center = c(0,0),draw = TRUE, ...){
  if(length(r)>1){
    return(sapply(r,quick_triangle,center=center,draw=draw,...))
  }
  ijk     <- matrix(data = c( 0  , sqrt(3)/3,
                              -0.5,-sqrt(3)/6,
                              0.5,-sqrt(3)/6), ncol = 2, nrow = 3, byrow = TRUE) * r + matrix(rep(center, 3), ncol=2, byrow=TRUE)
  if(draw)  polygon(ijk,...)
  return(ijk)
}
#' @title quick_tea_plot
#' @description creates a tea_plot given ternary data
#' @param data the data with which to create the ternary plot
#' @param dataLabels include if you want to plot text instead of symbols.  Defaults to NA
#' @param main the title for your plot
#' @param col takes normal colors plus "contrast", "gradient", "gradient90", "gradient180" and "gradient270"
#' @param cex as in plot, points and text functions, increases size of these.  Defaults to 1.5 here.
#' @param pch as in plot and points functions.  Determines symbol of points plotted.  Defautls to 17 (triangle)
#' @param ... any additional variables that may work wtih TeaSet's tea_plot() method or the points() function
#' @return NONE
#' @note essentially accesses tea_plot without having to create a TeaSet object oneself
#' @example \dontrun{quick_tea_plot(matrix(c(rnorm(300)),ncol=3))}
#' @export
quick_tea_plot <- function(data,dataLabels=NA,main = "",col = "gradient",pch=17,cex=1.5,...){
  teaSet<-brew_tea(data)
  teaSet$tea_plot(axis.labels = colnames(data),ticks=TRUE,main=main,pch=pch,cex=cex,...)
}
