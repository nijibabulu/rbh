#! /usr/bin/env python

import abc
import math
import itertools
import operator
import matrix

def median(data):
    data = sorted(data)
    n = len(data)
    mid = n//2
    if n % 2 == 1:
        return data[mid]
    else:
        return float(data[mid-1]+data[mid])/2

def median_abs_dev(data,center):
    return median([abs(p)-center for p in data])

def weighted_median_abs_dev(data,weight):
    sw = [k for k,w in sorted(zip(weight,data), key=operator.itemgetter(1))]
    sdata = list(sorted(data))
    z = float(sum(sw))
    s = 0.
    for n,p in enumerate(sw):
        s += p/z
        if s >= 0.5:
            break
    print n
    if p + sw[n]/z > 0.5:
        return sdata[n]/0.6745
    else:
        return (sdata[n] + sdata[n+1])/(2*0.6745)

def psi_bisquare(data,c=4.685):
    return [(1 - min(1,p/c)**2)**2 for p in data]

class LeastSquaresRegression(object):
    def __init__(self,x,y,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)
        self.x = x
        self.y = y
        self.Xm = matrix.Matrix([[1.,xi] for xi in x])
        self.ym = matrix.Matrix(y)
        self._compute_params()
        self.alpha,self.beta = self._params
        self.resid =  [yi-(xi*self.beta+self.alpha) 
                       for xi,yi in zip(self.x,self.y)]

    def _compute_params(self):
        X,y = self.Xm,self.ym
        params_m = X.T().dot(X).inverse().dot(X.T()).dot(y)
        self._params = list(itertools.chain.from_iterable(params_m.m))
      
class WeightedLeastSquaresRegression(LeastSquaresRegression):
    def __init__(self,x,y,w):
        super(WeightedLeastSquaresRegression,self).__init__(x,y,wm=matrix.Matrix.diag(w))

    def _compute_params(self):
        X,y,w = self.Xm,self.ym,self.wm
        params_m = X.T().dot(w).dot(X).inverse().dot(X.T()).dot(w).dot(y)
        self._params = list(itertools.chain.from_iterable(params_m.m))

class ObjectiveFunction(object):
    __metaclass__ = abc.ABCMeta
    @abc.abstractmethod
    def __init__(self,**kwargs): pass
    @abc.abstractmethod
    def psi(self,residuals): pass

class TukeyBiweight(ObjectiveFunction):
    def __init__(self,c=4.685):
        if c < 1.548:
            raise ValueError, 'TukeyBiweight requires a tuning const >= 1.548'
        self.c = c 
    def psi(self,residuals,scale):
        return [((1-min(1,abs(ri/scale/self.c))**2)**2) for ri in residuals]

class RobustLinearRegression(object):
    def __init__(self,x,y,objective=TukeyBiweight(),max_iterations=100,
                 convergence_condition=1e-4):
        self.x = x
        self.y = y
        self.obj = objective
        self.max_iterations = max_iterations
        self.converged = False
        self.wt = None
        self.tests = []
        self.eps = convergence_condition
        self._compute()

    def _compute(self):
        cur = WeightedLeastSquaresRegression(self.x,self.y,[1]*len(self.x))
        for i in range(self.max_iterations):
            scale = median_abs_dev(cur.resid,0)/0.6745
            old_resid = cur.resid
            w = self.obj.psi(cur.resid,scale)
            self.wt = w
            cur = WeightedLeastSquaresRegression(self.x,self.y,self.wt)
            #print old_resid
            #print w
            #print cur.resid
            #print sum([(ro-rn)**2 for ro,rn in zip(old_resid,cur.resid)])
            convi = math.sqrt(
                sum([(ro-rn)**2 for ro,rn in zip(old_resid,cur.resid)]
                       )/max(1e-20,sum(ro**2 for ro in old_resid))
            )
            print convi
            if convi < self.eps:
                self.converged = True
                break
        self.alpha = cur.alpha
        self.beta = cur.beta
        self.resid = cur.resid



if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        raise SystemExit, 'please give me a table'
    xvals = []
    yvals = []
    sds = []
    for line in open(sys.argv[1]):
        xval,yval,sd = line.split()
        try:
            xvals.append(float(xval))
            yvals.append(float(yval))
            sds.append(float(sd)**-2)
        except:
            continue
    print xvals
    print yvals 
    print sds
    ols = LeastSquaresRegression(xvals,yvals)#,sds)
    print ols.alpha,ols.beta,ols.resid
    wls = WeightedLeastSquaresRegression(xvals,yvals,sds)
    print wls.alpha,wls.beta,wls.resid
    rls = RobustLinearRegression(xvals,yvals)
    print rls.alpha,rls.beta
    #print compute_ols_regression(xvals,yvals)
    #print compute_wls_regression(xvals,yvals,sds)
    #print median([1,5,2,3,4])
    #print median([1,5,2,3,4,6])
    #print median_abs_dev([-1,5,-2,-3,4,6],0)
    #print weighted_median_abs_dev([-1,5,-2,-3,4,6,20],[1,2,3,1,2,2,4])
