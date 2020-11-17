# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 11:35:51 2018

@author: yadav
"""

#code for computation of parameters c1,c0 of Joukowsky transform for CW paper self consistently. See notebook for details
import math
import cmath
import numpy
import contour_integral
contr = numpy.loadtxt("input_output_files//contour//nu_contour_c1=1_c0=0.5_20000points_edit.txt",float)
t = 1    # V(x)=(x**2/2)
c1 = 5
c0 = 3    
jacobian = numpy.empty([2,2],complex)
ff = numpy.empty(2,complex)
c = numpy.empty(2,complex)
c[0] = c1
c[1] = c0
while True:
    def f1(z):    # for J(1,1), See notebook
            return (1.0/(2*(math.pi)*1j*t))*(z)
            
    jacobian[0,0] = contour_integral.contr_intgl(contr,f1) + (1.0/(c[0]**2))
    print ('J(1,1) is', jacobian[0,0])
       
    def f2(z):    # for J(1,2), See notebook
            return (1.0/(2.0*(math.pi)*1j*t))

    jacobian[0,1] = contour_integral.contr_intgl(contr,f2)
    print ('J(1,2) is', jacobian[0,1])            
            
    def f3(z):    # for J(2,1), See notebook
            return (1.0/(2.0*(math.pi)*1j*t))*(z/(z-1.0/2))
            
    jacobian[1,0] = contour_integral.contr_intgl(contr,f3)
    print ('J(2,1) is', jacobian[1,0])        
            
    def f4(z):    # for J(2,2), See notebook
            return (1.0/(2.0*(math.pi)*1j*t))*(1.0/(z-1.0/2))
            
    jacobian[1,1] = contour_integral.contr_intgl(contr,f4)
    print ('J(2,2) is', jacobian[1,1])

    def f5(z):    # for ff1, See notebook
            return (1.0/(2.0*(math.pi)*1j*t))*(c[0]*z + c[1] - cmath.log((z-1.0/2)/(z+1.0/2)))

    ff[0] = contour_integral.contr_intgl(contr,f5) - (1.0/c[0])  

    def f6(z):    # for ff2, See notebook
            return (1.0/(2*(math.pi)*1j*t))*((c[0]*z + c[1] - cmath.log((z-1.0/2)/(z+1.0/2)))/(z-1.0/2))

    ff[1] = contour_integral.contr_intgl(contr,f6) - 1.0

    delta_X = numpy.linalg.solve(jacobian,ff)            
    
    c_next = c - delta_X    
         
    error = c_next - c
    print('error is', error)
    
    c = c_next
    
#    if(numpy.amax(error)<=1e-2):
    if(abs(error[0])<=(1e-4) and abs(error[1])<=1e-10 ):
        break            