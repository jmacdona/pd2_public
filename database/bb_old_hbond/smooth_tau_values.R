tab = read.table("tau.dat")
smooth_tab = tab

x=as.numeric(tab$V4)
smoothingSpline = smooth.spline(x, tab$V5, spar=0.4)
plot(x,tab$V5)
#lines(smoothingSpline)
smooth_tab$V5 = smoothingSpline$y
lines(x,smooth_tab$V5)


smoothingSpline = smooth.spline(x, tab$V6, spar=0.4)
plot(x,tab$V6)
#lines(smoothingSpline)
smooth_tab$V6 = smoothingSpline$y
lines(x,smooth_tab$V6)


smoothingSpline = smooth.spline(x, tab$V7, spar=0.4)
plot(x,tab$V7)
#lines(smoothingSpline)
smooth_tab$V7 = smoothingSpline$y
lines(x,smooth_tab$V7)




write.table(smooth_tab, "smoothed_tau.dat", row.names = F, col.names = F, sep="\t", quote=F)

