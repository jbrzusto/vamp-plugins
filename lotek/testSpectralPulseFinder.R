##
## test the spectral pulse finder
##
## simple noise with superimposed pulses
##
## window size 10, padding x 3, overlap 5
##

## pulse width (samples)
plen = 100

## window size
n = 10

## number of windows of data to generate
nw = 1000

## padding factor (1 means no padding)
pf = 1

## overlap (must be in the range 0... window size - 1)
overlap = 5

## window type (0 = Hamming; 1 = rectangular)
wintype = 0

## noise amplitude

noiseamp = 0.01

## signal amplitude

sigamp = 0.5

bkgd = noiseamp * (rnorm(nw * n) + 1i * rnorm(nw * n))

## 5 random pulses of different offset frequency

for (i in 1:5) {
    p = sigamp * exp(2*pi*1i*(0:(plen-1))/plen*(2.567 * i))
    rep =  300 + 1000 * (i-1) + seq(along=p)
    bkgd [rep] = p + bkgd[rep]
}

## send to testSlidingSpectrum program

tf = tempfile()

p = pipe(sprintf("./testSpectralPulseFinder %d %d %d %d %d %d %f %f > %s", plen, n, n * (pf - 1), overlap, 1, n * pf, 10, 8, tf), "w")

cat (c(rbind(Re(bkgd),Im(bkgd))), file=p)

close(p)

res = read.csv(tf, as.is=TRUE, header=FALSE)

names(res) = c("count", "bin", "sig", "noise", "Z")

with(res, plot(count, sig, type="l"))

