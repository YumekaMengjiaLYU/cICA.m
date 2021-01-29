# cICA.m
R Package for colored independent component analysis that utilizes Morlet wavelet convolution developed by [Prof. Seonjoo Lee](https://sites.google.com/site/seonjool/).

### Install
with `devtools`:

```S
devtools:install_github('cICA.m')
```

### Use
There are only two functions in this package.

Call `cICA_morlet` to perform colorICA Morlet algorithm on CIFTI files. 

Call `drawplot` to plot the signal, time-frequency power and power spectrum of Morlet wavelet transformed information from CIFTI files.
