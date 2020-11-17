# -*- coding: utf-8 -*-
"""
Created on Sat Jun 01 14:14:07 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 03 23:22:24 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 03 12:27:35 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 07 15:35:06 2019

@author: yadav
"""

# code to calculate f2(y) self-consistently using Newton-Rapson method (see notebook or section 6.3.7 in computational physics book by mark newman)
# this code takes value of gamma and theta as user input from command line. It is to be used for HYDRA.
import math
import cmath
import numpy
import contour_integral
import integration_trapezoidal_nonuniformspacing
import derivative_central_difference
import timeit

start = timeit.default_timer()

gamma = float(raw_input("gamma is: "))    # asks for user input from command line. to be used in terminal    
rho = float(raw_input("rho is: "))    # V(x)=rho*x    # parameters from solved C.R. problem    # asks for user input from command line. to be used in terminal
#delta = float(raw_input("delta is: "))    # V(x)=rho*x    # parameters from solved C.R. problem    # asks for user input from command line. to be used in terminal
#alpha = float(raw_input("alpha is: "))    # V(x)=rho*x+alpha*x**(1.0/alpha), alpha>1    # asks for user input from command line. to be used in terminal 
c1 = float(raw_input("c1 is: "))    # from self-consistent calculation using jouwkowsky_parameters_c_selfconsistent_calculation_CR_problem.py for alpha-log(x) potential.
c1_short = float(raw_input("c1_short is: "))
c0 = float(raw_input("c0 is: "))    # from self-consistent calculation using jouwkowsky_parameters_c_selfconsistent_calculation_CR_problem.py for alpha-log(x) potential.
c0_short = float(raw_input("c0_short is: "))

#rho_delta = rho*(1.0+delta)

iteration = float(raw_input("iteration is: "))

#gamma = 0.9
#theta = 1.1    # parameters from solved C.R. problem 
#rho = 2.0    # parameters from solved C.R. problem, V(x)=rho*x
#c = 0.73346142451457241    # from self-consistent calculation using jouwkowsky_parameters_c_selfconsistent_calculation_CR_problem.py for alpha-log(x) potential.
epsi = 1e-4

data1 = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu1_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
data2 = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu2_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
data3 = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/mapping/mapping_output_nu1_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
data4 = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/mapping/mapping_output_nu2_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
contr = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu1_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
f_x_input = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/density/f_y_calculation_newton_rapson_method_2_gamma="+str(gamma)+"_9000points_correction4_iter"+str(iteration-1)+".txt",float)
f_out=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/density/f_y_calculation_newton_rapson_method_2_gamma="+str(gamma)+"_9000points_correction4_iter"+str(iteration)+".txt","w")

#data1 = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu1_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
#data2 = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu2_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
#data3 = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/mapping/mapping_output_nu1_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
#data4 = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/mapping/mapping_output_nu2_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
#contr = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu1_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
#f_x_input = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/density/f_y_delta_calculation_newton_rapson_method_2_gamma="+str(gamma)+"_9000points_correction4_iter"+str(iteration-1)+".txt",float)
#f_out=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/density/f_y_delta_calculation_newton_rapson_method_2_gamma="+str(gamma)+"_9000points_correction4_iter"+str(iteration)+".txt","w")


x = data4[:, 0]    # Array Slicing: Accessing array's first column(i.e. index-zero column), https://jakevdp.github.io/PythonDataScienceHandbook/02.02-the-basics-of-numpy-arrays.html
x1 = data1[:, 0]    # for nu1    # Note- to convert integer type array(x) to float type array(xx) Use: xx = x.astype(float)
y1 = data1[:, 1]    # for nu1
r1 = numpy.sqrt(x1**2+y1**2)
x2 = data2[:, 0]    # for nu2    # Note- to convert integer type array(x) to float type array(xx) Use: xx = x.astype(float)
y2 = data2[:, 1]    # for nu2
r2 = numpy.sqrt(x2**2+y2**2)

hbar_x = x2+y2*1j    # h(x) in formula in notebook
h_x = numpy.conjugate(hbar_x)    # hbar(x) in formula in notebook 
#h_y = x1*(1.0+epsi/r1)+y1*(1.0+epsi/r1)*1j    # Adding small epsilon radially, h(y) in formula in notebook
h_x_plus = (h_x.real)*(1.0+epsi/r2)+(h_x.imag)*(1.0+epsi/r2)*1j    # Adding small epsilon radially, h(y) in formula in notebook
hbar_x_plus = (x2)*(1.0+epsi/r2)+(y2)*(1.0+epsi/r2)*1j    # hbar(y) in formula

deriv_hbar_real = derivative_central_difference.derivative_cent_diff_(x,x2)
deriv_hbar_imag = derivative_central_difference.derivative_cent_diff_(x,y2)
deriv_hbar = deriv_hbar_real + deriv_hbar_imag*1j

deriv_V = (2.0/rho)*x    # V(x)=(x^2)/rho
#deriv_V = (2.0/rho_delta)*x    # V(x)=(x^2)/rho_delta

jacobian = numpy.empty([len(x),len(x)],float)
F = numpy.empty(len(x),float)

#f_x = ((deriv_V*x)/gamma)-((1.0-gamma)/gamma)     # Step(1): Initialization
f_x_2 = f_x_input[:, 1]    # f(x) from previous iteration

# Step(2): Calculation of Jacobian
for i in range(len(x)):    # function len() on array gives no. of rows of array
    for k in range(len(x)):
        if(k==(len(x)-1)):
            funct2 = ((1.0-gamma)/(2*(numpy.pi)*gamma))*((((1.0/(h_x_plus[i]-hbar_x[k]))+(1.0/(hbar_x_plus[i]-hbar_x[k])))*deriv_hbar[k]).imag)*(x[k]-x[k-1])
        else:
            funct2 = ((1.0-gamma)/(2*(numpy.pi)*gamma))*((((1.0/(h_x_plus[i]-hbar_x[k]))+(1.0/(hbar_x_plus[i]-hbar_x[k])))*deriv_hbar[k]).imag)*(x[k+1]-x[k])    
        if(i==k):
            jacobian[i,k] = 1 + funct2
        else:
            jacobian[i,k] = funct2
#        print i,k    

#while True:            
    # Step(3): F(X) calculation
for l in range(len(x)):    # function len() on array gives no. of rows of array
    first_term = f_x_2[l]    # first term in formula for F(X) in notebook
    second_term = -((deriv_V[l])/gamma)    # second term in formula for F(X) in notebook
#    second_term = -((f_x[l])/gamma)+((1.0-gamma)/gamma)    # second term in formula for F(X) in notebook
#    funct1 = ((1.0-gamma)/(2*(numpy.pi)*gamma))*f_x*((((1.0/(h_x_plus[l]-hbar_x))+(1.0/(hbar_x_plus[l]-hbar_x)))*deriv_hbar).imag)
    funct1 = ((1.0-gamma)/(2*(numpy.pi)*gamma))*f_x_2*((((1.0/(h_x_plus[l]-hbar_x))+(1.0/(hbar_x_plus[l]-hbar_x)))*deriv_hbar).imag)
    third_term = integration_trapezoidal_nonuniformspacing.int_trap_non_uniform(x,funct1)    # third term in formula for F(X) in notebook
    F[l] = first_term + second_term + third_term
#        print l

    # Step(4): Solving for delta_X in eqn.(6.108) in computational physics book by mark newman 
delta_X = numpy.linalg.solve(jacobian,F)   

    # Step(5): Next estimate
f_x_next = f_x_2 - delta_X
max_error = numpy.amax(abs(f_x_next - f_x_2))
print('maximum error is', max_error)
f_x_2 = f_x_next
#    if(max_error<=1e-10):
#        break

# Step(5): writing output to a file
for m in range(len(x)):    # function len() on array gives no. of rows of array
    f_out.write(str(x[m])+" "+str(f_x_2[m])+'\n')

f_out.close()    # () at the end is necessary to close the file                

stop = timeit.default_timer()
print("run time is: ", stop - start)            