
def gradient(x, y, a):
    import numpy as np
    grad = np.zeros(4)
    t0 = -x/a[2+0]; t1 = -x/a[2+1]
    e0 = np.exp(-x/a[2+0])
    e1 = np.exp(-x/a[2+1])
    dev = y-a[0]*e0-a[1]*e1
    grad[0] = sum(dev*e0)
    grad[1] = sum(dev*e1)
    grad[2+0] = a[2+0]*sum(dev*x*e0)/a[2+0]**2
    grad[2+1] = a[2+1]*sum(dev*x*e1)/a[2+1]**2
    grad *= 2
    return grad
def secondDerivativeMatrix(x, y, a):
    import numpy as np
    matrix = np.ones((4,4))
    e0 = np.exp(-x/a[2+0])
    e1 = np.exp(-x/a[2+1])
    dev = a[0]*e0+a[1]*e1 - y
    matrix[0,0] = sum(e0**2)
    matrix[0,1] = sum(e0*e1)
    matrix[0,2] = sum(x*e0*(a[0]*e0+dev))/a[2+0]**2
    matrix[0,3] = a[1]*sum(x*e0*e1)/a[2+1]**2
    matrix[1,1] = sum(e1**2)
    matrix[1,2] = a[0]*sum(x*e0*e1)/a[2+0]**2
    matrix[1,3] = sum(x*e1*(a[1]*e1+dev))/a[2+1]**2
    matrix[2,2] = a[0]*sum(x*e0*(x*e0*a[0]/a[2+0]+dev*(x/a[2+0]-2)))/a[2+0]**3
    matrix[2,3] = a[0]*a[1]*sum(x**2*e0*e1)/a[2+0]**2/a[2+1]**2
    matrix[3,3] = a[1]*sum(x*e1*(x*e1*a[1]/a[2+1]+dev*(x/a[2+1]-2)))/a[2+1]**3
    matrix[1,0] = matrix[0,1]
    matrix[2,0] = matrix[0,2]
    matrix[2,1] = matrix[1,2]
    matrix[3,0] = matrix[0,3]
    matrix[3,1] = matrix[1,3]
    matrix[3,2] = matrix[2,3]
    matrix *= 2
    return matrix
