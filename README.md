# icewebPy
Python re-implementation of AVO IceWeb application including spectrogram browser. Some minimal examples follow.

## Creating an icewebSpectrogram object

Assuming you have this repository at /path/to/repo and an ObsPy Stream object (st), and that it is in units of m/s.

<pre>
import sys
sys.append('/path/to/repo')
import IceWeb

spobj = IceWeb.icewebSpectrogram(stream=st)
</pre>

## Plotting Options:

We plot an icewebSpectrogram object by calling the plot method. Here is the function prototype:

<pre>
    def plot(self, outfile=None, secsPerFFT=None, fmin=0.5, fmax=20.0, log=False, cmap=pqlx, clim=None, \
                      equal_scale=False, title=None, add_colorbar=True, precompute=False, dbscale=True)
</pre>




### 1 - Unscaled, amplitude units

All we have to do is create an instance of an icewebSpectrogram object, and then call the plot method. Each spectrogram is individually scaled to make best use of the colormap.

<pre>
  sgramfile = 'myspecgram_unscaled.png'
  spobj.plot(outfile=sgramfile)
</pre>

![2005-05-01-1344-36S MVO___025_sgram](https://user-images.githubusercontent.com/233816/122985285-fc224d80-d36b-11eb-8c07-36480b9f4234.png)

### 2 - Best overall scale, amplitude units

As in 1, but we want to choose the best overall spectral amplitude scale so all spectrogram plots are scaled (colored) the same:

<pre>
  sgramfile = 'myspecgram_scaled.png'
  spobj.plot(outfile=sgramfile, equal_scale=True)
</pre>

![2005-05-01-1344-36S MVO___025_sgram_scaled](https://user-images.githubusercontent.com/233816/122985327-0a706980-d36c-11eb-8ef8-7e7dbc8eb96a.png)


### 3 - Fixed overall scale, amplitude units

As in 2, but we want to provide our own scale (colormap) limits. This is the default for AVO North spectrograms:

<pre>
sgramfile = 'myspecgram_fixed.png'
spobj.plot(outfile=sgramfile, clim=[1e-8, 1e-5])
</pre>

Note that the scale here is in units of m/s/Hz. 

![2005-05-01-1344-36S MVO___025_sgram_fixed](https://user-images.githubusercontent.com/233816/122985344-0fcdb400-d36c-11eb-9a6e-2d57047f38d4.png)



### 4 -  Unscaled, decibel (dB) units 
<pre>
  sgramfile = 'myspecgram_unscaled.png'
  spobj.plot(outfile=sgramfile, dbscale=True, title='Unscaled')
</pre>

![2005-05-01-1344-36S MVO___025_sgram](https://user-images.githubusercontent.com/233816/122984730-55d64800-d36b-11eb-92a9-40b6d3f7c6c0.png)



### 5 -  Best overall scale, decibel (dB) units
<pre>
  sgramfile = 'myspecgram_unscaled.png'
  spobj.plot(outfile=sgramfile, equal_scale=True, dbscale=True, title='Scaled')
</pre>

![2005-05-01-1344-36S MVO___025_sgram_scaled](https://user-images.githubusercontent.com/233816/122984751-5bcc2900-d36b-11eb-8374-cb812576db0c.png)



### 6 -  Fixed overall scale, decibel (dB) units
<pre>
  sgramfile = 'myspecgram_unscaled.png'
  spobj.plot(outfile=sgramfile, clim=[1e-8, 1e-5], dbscale=True, title='Fixed')
</pre>

![2005-05-01-1344-36S MVO___025_sgram_fixed](https://user-images.githubusercontent.com/233816/122984777-62f33700-d36b-11eb-8f62-10e759cbd092.png)



## Changing the colormap
To change the colormap, pass the optional cmap name:value pair to the plot method.

The default colormap is pqlx. Other options are viridis_white_r, obspy_divergent, obspy_sequential

<pre>
  from obspy.imaging.cm obspy_sequential
  spobj.plot(..., cmap=obspy_sequential )
</pre>

![2005-05-01-1344-36S MVO___025_sgram_fixed](https://user-images.githubusercontent.com/233816/122985952-c5990280-d36c-11eb-8a71-e316f65e5672.png)


## Making the y-scale logarithmic
<pre>
  spobj.plot(..., log=True )
</pre>

![2005-05-01-1344-36S MVO___025_sgram_fixed](https://user-images.githubusercontent.com/233816/122986284-3d672d00-d36d-11eb-9e49-6aad2b0ccf8e.png)


## Changing the frequency limits on the y-axis
<pre>
  spobj.plot(..., fmin = 0.0, fmax = 10.0 )
</pre>

## Adding a title
<pre>
  spobj.plot(..., title = 'Montserrat 2005-05-31 13:46:36' )
</pre>

## Removing colorbars
<pre>
  spobj.plot(..., add_colorbar = False)
</pre>

## Pre-computing the spectral data

If plotting the same seismic data as spectrograms in different ways, it can be useful to pre-compute the spectrograms:

<pre>
  spobj = IceWeb.icewebSpectrogram(stream=st)
  spobj = spobj.precompute()
  
  # plot same data in 4 different ways
  spobj.plot(dbscale=False)
  spobj.plot(dbscale=False, equal_scale=True)
  spobj.plot(dbscale=True)
  spobj.plot(dbscale=True, equal_scale=True)
</pre>


Glenn Thompson, 2021/06/22
