from __future__ import division
from imports import *
import nhqm.bases.many_body_spill as mb
from scipy.sparse import lil_matrix
import matplotlib.cm as cm
import time
t=time.time()
Hdict = mb.get_hamilton_dict(15, 4)
print len(Hdict)
print time.time()-t

#ord=len(Hdict)

fig = plt.figure(2)
ax = fig.add_subplot(1,1,1)
#H=lil_matrix((ord,ord))
for (i,j), lista in Hdict.iteritems():
    w=len(lista)
    #print w
    if w==10:
        c='k'
    else:
        c=cm.gray(w/4.,1)
    ax.plot(i,j,'s',markersize=3,color=c)
#plt.pcolor(H.todense())
#plt.gca().invert_yaxis()
ax.invert_yaxis()
ax.xaxis.set_ticks_position('top')
plt.show()
