import utilities as li
import math as m

def function(y,x):
    f=(y*m.log(y))/x
    return(y)
def func1(x):
    s=m.exp(x/2)
    return(s)
fun=func1
func=function
h=[0.01,0.02,0.1,0.5]
li.Explicit_euler(2,2.71828,h,func,fun,10)

#sub question b
def function2(y,x):
    s=6-((2*y)/x)
    return (s)

def func2(x):
    f=2*x-45/x**2
    return(f)
fun2=function2
func2=func2

li.Explicit_euler(3,1,h,fun2,func2,10)


#output: graph plot attached
