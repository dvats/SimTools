plot.boxx <- function(x,dimn, CIs, mean.color, quan.color, mn, quans,...,range, width, varwidth, notch, outline, names, plot, border,
                      col, ann, horizontal, add) 
{
  boxplot(x,...,range=range,width=width,varwidth=varwidth,notch=notch,outline=outline,names=names,plot=plot,border=border,
          col=col,ann=ann,horizontal=horizontal,add=add)
  for(i in 1:dimn) {
    quansi <- quans[, i]
    qcil = CIs$lower.ci.mat[, i]
    qciu = CIs$upper.ci.mat[, i]
    for(j in 1:length(quansi))
    {
      if(horizontal==TRUE) {
        polygon(c(qcil[j],qcil[j],qciu[j],qciu[j]), c(i-(0.2*min(dimn,2)),i+(0.2*min(dimn,2)),i+(0.2*min(dimn,2)), i-(0.2*min(dimn,2))),  col = quan.color, border = FALSE)
      }else {
        polygon(c(i-(0.2*min(dimn,2)),i+(0.2*min(dimn,2)),i+(0.2*min(dimn,2)), i-(0.2*min(dimn,2))), c(qcil[j],qcil[j],qciu[j],qciu[j]), col = quan.color, border = FALSE)
      } 
    }
    for(j in 1:length(quansi))
    {
      if(horizontal==TRUE) {
        segments(quansi[j],i-(0.2*min(dimn,2)), quansi[j], i+(0.2*min(dimn,2)))
      }else {
        segments(i-(0.2*min(dimn,2)), quansi[j], i+(0.2*min(dimn,2)), quansi[j])
      }
      
    }
  }
}

boxplot.Siid <- function(x, Q = c(0.25, 0.50,0.75), thresh = 0.001, alpha = 0.05,mean.col = 'plum4',
                      quan.col = 'lightsteelblue3',...,range = 1.5, width = NULL, varwidth = FALSE,
                      notch = FALSE, outline = TRUE, names, plot = TRUE,
                      border = par("fg"), col = 'white',
                      ann = !add, horizontal = FALSE, add = FALSE, opaq = .6)
{
  foo3 <- error.est(x,Q, alpha,thresh = thresh,iid = TRUE,...)
  plot.boxx(x, dimn = length(x[1,]), CIs = foo3, mean.color = mean.col, quan.color = adjustcolor(quan.col, alpha.f = opaq), mn = foo3$mean.est,quans = foo3$xi.q
            ,...,range=range,width=width,varwidth=varwidth,notch=notch,outline=outline,names=names,plot=plot,border=border,
            col=col,ann=ann,horizontal=horizontal,add=add)
}
