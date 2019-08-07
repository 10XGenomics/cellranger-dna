import numpy as np

##
## Support functions for clustering
##

def diff( v, window=5 ):
    assert len(v) > 1
    x = np.zeros( len(v) )
    for i in xrange(len(v)):
        low = i-window
        high = i+window
        if low < 0 or high > len(v):
            continue
        x[i] = ( sum( v[i:high] ) - sum( v[low:i] ) )
    return x


def coarse_grain( X, window ):
    C = np.cumsum(X, axis=1)
    Z = np.zeros(X.shape)
    window = 5
    assert window % 2 == 1
    n = X.shape[-1]
    for i in xrange(n):
        imin = max(0, i - window/2)
        imax = min(i + window/2, n-1 )
        Z[:, i] = C[:, imax] - C[:, imin]
    return Z

def sum_window( X, window ):
    assert len(X.shape) == 2
    Y = []
    for i in xrange(X.shape[0]):
        y = []
        for j in xrange( 0, X.shape[1], window ):
            y.append( X[i][j:j+window].sum() )
        Y.append(y)
    Y = np.array(Y)
    return Y

