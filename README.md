# icewebPy
Python re-implementation of AVO IceWeb application including spectrogram browser

The first 6 examples below focus on the different options for applying the colormap.

## Minimal example 1 - Unscaled, amplitude units

Assumes you have only a Stream object (here called stZ), and that it is in units of m/s.
All we have to do is create an instance of an icewebSpectrogram object, and then call the plot method. Each spectrogram is individually scaled to make best use of the colormap.

<pre>
spobj = IceWeb.icewebSpectrogram(stream=stZ)

sgramfile = 'myspecgram_unscaled.png'
spobj.plot(outfile=sgramfile)
</pre>

## Minimal example 2 - Best overall scale, amplitude units

As in 1, but we want to choose the best overall spectral amplitude scale so all spectrogram plots are scaled (colored) the same:

<pre>
sgramfile = 'myspecgram_scaled.png'
spobj.plot(outfile=sgramfile, equal_scale=True)
</pre>

## Minimal example 3 - Fixed overall scale, amplitude units

As in 2, but we want to provide our own scale (colormap) limits. This is the default for AVO North spectrograms:

<pre>
sgramfile = 'myspecgram_fixed.png'
spobj.plot(outfile=sgramfile, clim=[1e-8, 1e-5])
</pre>

Note that the scale here is in units of m/s/Hz. 


## Minimal example 4 -  Unscaled, decibel (dB) units
<pre>
sgramfile = 'myspecgram_unscaled.png'
spobj.plot(outfile=sgramfile, dbscale=True)
</pre>

![2005-05-01-1344-36S MVO___025_sgram](https://user-images.githubusercontent.com/233816/122980196-60daa980-d366-11eb-9349-b574a09701d5.png)


## Minimal example 5 -  Best overall scale, decibel (dB) units
<pre>
sgramfile = 'myspecgram_unscaled.png'
spobj.plot(outfile=sgramfile, equal_scale=True, dbscale=True)
</pre>

![2005-05-01-1344-36S MVO___025_sgram_scaled](https://user-images.githubusercontent.com/233816/122980222-67692100-d366-11eb-8802-a13cde04744e.png)


## Minimal example 6 -  Fixed overall scale, decibel (dB) units
<pre>
sgramfile = 'myspecgram_unscaled.png'
spobj.plot(outfile=sgramfile, clim=[1e-8, 1e-5], dbscale=True)
</pre>

![2005-05-01-1344-36S MVO___025_sgram_fixed](https://user-images.githubusercontent.com/233816/122980247-6e902f00-d366-11eb-85ef-53ade9fbca4e.png)


## Changing the colormap
To change the colormap, pass the optional cmap name:value pair to the plot method.
spobj.plot(..., cmap=pqlx )
The default colormap is pqlx. Other options are viridis_white_r, obspy_divergent, obspy_sequential

