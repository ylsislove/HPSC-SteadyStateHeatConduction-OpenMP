from matplotlib import pyplot as plt
from matplotlib import cm  # colormaps
import numpy as np


try:
    fname = 'solution.txt'
    x, y, u = np.loadtxt(fname, unpack=True)
    # unpack=True, If True, the returned array is transposed, arrays are returned for each field.
except:
    err_msg = "Could not load data from file %s." % fname \
              + " Did you forget to run the program?"
    raise Exception(err_msg)

# Solution is plotted on n by n grid so length of each vector should be n**2
# Determine n:

n = int(np.sqrt(len(x)))
assert n * n == len(x), "Expected len(x) to be a perfect square, len(x) = %s" % len(x)

X = x.reshape(n, n)
Y = y.reshape(n, n)
U = u.reshape(n, n)

# Pseudocolor plot

plt.figure(1)
plt.clf()  # clear figure
plt.axis('scaled')  # so x- and y-axis scaled the same (square)
plt.pcolor(X, Y, U, cmap=cm.jet)  # pseudo-color plot using colormap "jet"
# https://matplotlib.org/gallery/color/colormap_reference.html
plt.xlim((0., 1.))  # x range
plt.ylim((0., 1.))  # y range
plt.clim(0., 1.)  # colors range from u=0 to u=1
plt.colorbar(ticks=np.linspace(0., 1., 11))  # add a color bar to show temperature scale, 11 samples array(0.1 intervel)
# plt.colorbar(ticks=np.arange(0., 1.1, 0.1))
plt.title('Temperature')

plt.savefig('pcolor.png')
print('Saved pseudocolor plot as pcolor.png')

# Contour plot

plt.figure(2)
plt.clf()

# should be commented
# plt.axis('scaled')

# contour line levels:
clines = np.linspace(0., 1., 26)

# do contour plot:
C = plt.contour(X, Y, U, clines, colors='k')
#     b : blue. g : green. r : red. c : cyan.  m : magenta. y : yellow. k : black. w : white.


# add labels on every other line:
# clines[0] & clines[-1] is removed from the QuadContourSet object C in the new version Matplotlib
# Hence we replace clines[1::2] with C.levels[1::2]
plt.clabel(C, C.levels[1::2], inline=1, fontsize=10)
# levels - A list of level values, that should be labeled.
# a = [1,2,3,4,5,6,7,8,9] a[1::2]=[2, 4, 6, 8]
# inline - If True the underlying contour is removed where the label is placed. Default is True.
# fontsize - Size in points

plt.title('Contours of temperature')

plt.savefig('contour.png')
print('Saved contour plot as contour.png')

plt.close()