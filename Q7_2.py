import utilities as li
import math as m

#function dy/dx=z
def func(x,y,z):
    return(z)
def function(x,y,z):
    g=-z+1-x
    return(g)
#function analytic
def fun(x):
    j=-1+m.exp(-x)-(x**2)/2+2*x
    return(j)
f=func
g=function
analy=fun
h=[-0.02,0.02,0.01,0.1]
li.RK4(0,2,1,h,f,g,analy,10)
