from pylab import *
import dnest4
import corner
from matplotlib import rcParams

import corner

rcParams["font.size"] = 16
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"



import os

import pyfits



# Piecewise linear stretch
def stretch(x):
	y = x.copy()
	y = (y - y.min())/(y.max() - y.min())
	y[y > 0.1] = 0.1 + 0.05*(y[y > 0.1] - 0.1)
	return y

saveFrames = False # For making movies
if saveFrames:
	os.system('rm Frames/*.png')



hdulist = pyfits.open('Data/fitim.fits')
data = hdulist[0].data # assuming the first extension is a table
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
imshow(stretch(data), cmap='copper', interpolation='nearest')
filename = 'test.png'
plt.savefig(filename)
filename = 'test.pdf'
plt.savefig(filename)
filename = 'test.eps'
plt.savefig(filename)

hdulist2 = pyfits.open('Data/sigma.fits')
sig = hdulist2[0].data # assuming the first extension is a table

# Data::get_instance().load("Data/test_metadata.txt", "Data/psfim.fits", "Data/sigma.fits");
#data = loadtxt('Data/test_image.txt')
#sig = loadtxt('Data/test_sigma.txt')


#hdulist = pyfits.open('Data/fitim.fits')
#prihdr = hdulist[0].header
#print(prihdr)
# scidata = hdulist['SCI'].data

#data = hdulist[1].data

#hdulist2 = pyfits.open('Data/sigma.fits')

#sig = hdulist2[1].data
data_shape = data.shape
xsize =  data_shape[0]
ysize =  data_shape[1]

posterior_sample = atleast_2d(dnest4.my_loadtxt('posterior_sample.txt', single_precision=True))

#print( type(posterior_sample[0,:]))
#print( len(posterior_sample[0,:]))
#print( type(posterior_sample.shape[0,:]))


ion()
hold(False)


print ("Posterior shape: ", posterior_sample.shape)
print ("Lenth: ", posterior_sample.shape[1])
print ("Starting Position: ", posterior_sample.shape[1] - xsize*ysize)
starting = posterior_sample.shape[1] - xsize*ysize


for i in range(0, posterior_sample.shape[0]):


	subplot(1, 3, 1)
	imshow(data, cmap='copper', interpolation='nearest')
	title('Original Data')
	gca().set_xticks([-0.5, 49.5, 99.5])
	gca().set_yticks([-0.5, 49.5, 99.5])
	gca().set_xticklabels(['0','0.5', '1'])
	gca().set_yticklabels(['1','0.5', '0'])

	img = posterior_sample[i, starting:starting+xsize*ysize].reshape((xsize, ysize))            # Need to edit this 
	subplot(1, 3, 2)
	imshow(img, cmap='copper', interpolation='nearest')
	title('Model {i}'.format(i=i))
	gca().set_xticks([-0.5, 49.5, 99.5])
	gca().set_yticks([-0.5, 49.5, 99.5])
	gca().set_xticklabels(['0','0.5', '1'])
	gca().set_yticklabels(['1','0.5', '0'])


	subplot(1, 3, 3)
	sigma = sqrt(sig**2 + posterior_sample[i,-1]**2)
	imshow(-(img - data)/sigma, cmap='coolwarm', interpolation='nearest')
	title('Standardised Residuals')
	gca().set_xticks([-0.5, 49.5, 99.5])
	gca().set_yticks([-0.5, 49.5, 99.5])
	gca().set_xticklabels(['0','0.5', '1'])
	gca().set_yticklabels(['1','0.5', '0'])


	draw()

	if saveFrames:
		savefig('Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
		print('Frames/' + '%0.4d'%(i+1) + '.png')

ioff()
show()

print('This is Davids add part to process the data')
# obtain subset without images     [0xcen, 1ycen, 2mag, 3re, 4nser, 5axrat, 6ang, 7box, mag8, rout9,a10, b11, axrat12, ang13, box14] 

#index = [2,3,4,5,6]
#index = [2,3,4,5,6,8,9,10,11,12,13]
subset = posterior_sample[:, 0:starting]
 #labels=["$m$", "$b$", "$\ln\,f$"],
  #                    truths=[m_true, b_true, np.log(f_true)])

# make triangle of disc
#fig = corner.corner(subset,labels=["$mag$", "$re$", "$n$","$q$","$\theta$"],truths=[18.07, 3.85, 0.16, 0.26, 14.9])          ,truths=[18.07,None,None, 0.26,14.9]
#fig = corner.corner(subset,labels=["$mag$", "$re$", "$n$","$q$","$\theta$"])
#truths  = [m,   re ,n, e , theta, f 
#truths1 = [18.07, 3.85, 0.16, 0.26, 14.9, 0.43 ]
#truths2 = [17.07, 2.38, 2.24, 0.65, -45.5, 0.57, 17.15] 



#subset = subset[:, index]
print('Shape:', subset.shape)

#disc
print('Make temporary cornerplot for the disc:')
index = [2,3,4,5,6]
labelz = ["$mag$", "$re$", "$n$","$q$","$\theta$"]
fig = corner.corner(subset[:,index],labels=labelz,quantiles=None,plot_contours=False)
fig.savefig("disc_triangle.png")
close(fig)


print('Make temporary cornerplot for the bar')
index = [8,9,10,11,12,13]
labelz = ["$mag$", "$Rout$", "$a$","$b$","q","$\theta$"]
fig = corner.corner(subset[:,index],labels=labelz,quantiles=None,plot_contours=False)
fig.savefig("bar_triangle.png")






