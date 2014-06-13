##
## test the narrow-band sliding spectrum code
##
## simple FM example
##

## window size
n = 10

## number of samples of data to generate
ns = 1000

## main frequency in spectrum
f = 5

## low frequency to watch
flow = 2 / 1000

## high frequency to watch
fhigh = 8 / 1000

##
## modulation frequency
f2 = 0.03

## modulation amplitude
amp = 0.1

## i and q components
t = (0:(ns - 1)) * (2 * pi / n)

## modulated frequency
fm = cumsum(f * (1 + amp * sin(f2 * t))) / (1:length(t))

## i and q components
xi = cos(fm * t)
xq = sin(fm * t)

## send to testSlidingSpectrum program

tf = tempfile()

p = pipe(sprintf("./testNarrowBandSlidingSpectrum 1 %d %f %f  > %s", n, flow, fhigh, tf), "w")

cat (c(rbind(xi, xq)), file=p)

close(p)

res = read.csv(tf, as.is=TRUE, header=FALSE)

## number of overlapped windows
pwr = matrix(res[,3], nrow = ns - n, byrow=TRUE)

## plot
library(rgl)
material3d(col="black")
persp3d(pwr, ylab="Frequency", xlab="Time", zlab="Power", col="lightblue")


