# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 23:42:12 2019

@author: yadav
"""
"""
Created on Thu Nov 01 16:35:31 2018

@author: yadav
"""

# code to generate complex mapping J(s) for soft edge

#https://stackoverflow.com/questions/16669428/process-very-large-20gb-text-file-line-by-line
#modified from https://stackoverflow.com/questions/30216573/reading-specific-columns-from-a-text-file-in-python
#https://docs.python.org/3.3/tutorial/inputoutput.html
# code to generate complex mapping J(s)

#https://stackoverflow.com/questions/16669428/process-very-large-20gb-text-file-line-by-line
#modified from https://stackoverflow.com/questions/30216573/reading-specific-columns-from-a-text-file-in-python
#https://docs.python.org/3.3/tutorial/inputoutput.html
import math
import cmath
import numpy

rho = 2.0
#delta = 0.0001
t=1.0
#t = 1.0+delta    #V(x)=x^2/(2t)

c1 = t    # for V(x)=x^2/(2t). See CW pg. no. 1071
c1_short = c1

c0 = t/2.0    # for V(x)=x^2/(2t). See CW pg. no. 1071
c0_short = c0

#c1 = 0.82454510531806569
#c1_short = 0.824

#c0 = 0.27671552249123316
#c0_short = 0.276

gamma = 0.4
iteration = 1.0

n = 100000    # number of points in each quadrant

#c1 = 50.0
#c0 = 25.0

f_in=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu1_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","r")
f_out=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/mapping/mapping_output_nu1_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")

#f_in=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu1_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","r")
#f_out=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/mapping/mapping_output_nu1_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")

lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))
    x_nu1= c1*z+c0-cmath.log((z-1.0/2)/(z+1.0/2))
    #f_out.write(str(x_nu1.real)+'\n')
    f_out.write(str(x_nu1.real)+" "+str(x_nu1.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file

f_in=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu2_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","r")
f_out=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/mapping/mapping_output_nu2_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")

#f_in=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu2_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","r")
#f_out=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/mapping/mapping_output_nu2_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")

lines=f_in.readlines()
#i=0
for data in lines:
    x= data.split(' ')[0]  
    y= data.split(' ')[1]  
    #z(i)=x(i)+y(i)j    gives syntax error
    z=complex(float(x),float(y))
    x_nu2= c1*z+c0-cmath.log((z-1.0/2)/(z+1.0/2))
    #f_out.write(str(x_nu2.real)+'\n')
    f_out.write(str(x_nu2.real)+" "+str(x_nu2.imag)+'\n')
    #i+=1
f_in.close()    # () at the end is necessary to close the file
f_out.close()    # () at the end is necessary to close the file