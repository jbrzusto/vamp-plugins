##
## test the sliding spectrum code
##
## simple FM example
##
## window size 10, padding x 3, overlap 5
##

## window size
n = 10

## number of windows of data to generate
nw = 100

## padding factor (1 means no padding)
pf = 10

## overlap (must be in the range 0... window size - 1)
overlap = 5

## window type (0 = Hamming; 1 = rectangular)
wintype = 0

## main frequency in spectrum
f = 3

## modulation frequency
f2 = 0.05

## modulation amplitude
amp = 0.3

## i and q components
t = (0:(n * nw - 1)) * (2 * pi / n)

## modulated frequency
fm = cumsum(f * (1 + amp * sin(f2 * t))) / (1:length(t))

## i and q components
xi = cos(fm * t)
xq = sin(fm * t)

## send to testSlidingSpectrum program

tf = tempfile()

p = pipe(sprintf("./testSlidingSpectrum %d %d %d %d > %s", n, n * (pf - 1), overlap, wintype, tf), "w")

cat (c(rbind(xi, xq)), file=p)

close(p)

res = read.csv(tf, as.is=TRUE, header=FALSE)

## number of overlapped windows
now = (nw * n - n) / (n - overlap)
pwr = matrix(res[,3], n * pf, now)

## plot
library(rgl)
material3d(col="black")
persp3d(pwr, xlab="Frequency", ylab="Time", zlab="Power", col="lightblue")


