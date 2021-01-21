#'
#'
#' This function plots the signal, time-frequency power and power spectrum of Morlet wavelet transformed information from CIFTI files.
#' @import ggthemes
#' @import gridExtra
#'
#' @param obj morlet transformed image data object
#' @param interpolate interpolation
#' @return return(list(a,fig1,fig2,fig3))
#' @export

drawplot<-function(obj,interpolate=TRUE){
  nt = length(obj$y)
  y.c = scale(obj$y, center=TRUE, scale=FALSE)
  mat = as.matrix(Mod(obj$Power))
  colnames(mat)<-(obj$period)
  rownames(mat)<-obj$x
  perid.x = obj$period
  perid.x0=1/((1:(nt/2-1))/nt)
  y.c.fft = Mod(fft(y.c)[2:(nt/2)])^2/(2*pi*nt)
  lspec.fit = lspec(period = y.c.fft)
  pred = dlspec( 1/perid.x *2*pi,lspec.fit)
  lspec.est = pred$d + pred$m

  mat.long = melt(mat) %>%
    rename(Period=Var2, Time=Var1)
  fig1<-ggplot(data.frame(x=obj$x, y=y.c),aes(x=x,y=y))+
    geom_line() + theme_bw() + xlab('Time')
  fig3<-ggplot(rbind(data.frame(x=perid.x0, y=y.c.fft, g='periodogram'),
                     data.frame(x=perid.x, y=lspec.est,g='lspec')),
               aes(x=x,y=y,colour=g)) + geom_line() +
    theme_bw() + xlab('') + coord_flip() + scale_x_continuous(trans = 'log2',limits=range(obj$period)) +
    theme(legend.position='bottom',
          legend.box.margin=margin(5,5,5,5))+ ggtitle('Power spectrum') + ylab('Power')

  fig2<-ggplot(mat.long, aes(x=Time, y=Period, fill=value)) +
    geom_raster(interpolate=interpolate) +
    scale_y_continuous(trans = 'log2') +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    theme_bw()+ ggtitle('Time-frequecy Power (Mod)') +
    theme(legend.position='bottom')

  lay <- rbind(c(1,1,1,NA),
               c(2,2,2,3),
               c(2,2,2,3),
               c(2,2,2,3))
  a = grid.arrange(fig1,fig2,fig3,
                   layout_matrix = lay)
  return(list(a,fig1,fig2,fig3))
}
