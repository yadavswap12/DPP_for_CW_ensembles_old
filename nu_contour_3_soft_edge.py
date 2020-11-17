# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 15:30:18 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:36:22 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:14:32 2018

@author: yadav
"""
"""
Created on Fri Nov 16 09:50:09 2018

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=c1s+c0-log((s-1/2)/(s+1/2))
#V(x)=x^2/2 ; so c1=1,c0=1/2
#Due to given from of locus equation we first find max and min values of y(here for x=0)
#then we vary y such that we get counterclockwise contour
#stepsize for each quarter of plot can be chosen
import math
import cmath
import numpy
import matplotlib    #pylab is submodule in matplotlib

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

#r=0.6532711880366251
#c1 = 50.0    # to check hard edge behaviour as c1->inf so let c1=50,c0=25, so V(x)=x^2/100
#c0 = 25.0   
#n = 10000    # number of points in each quadrant

f_out=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")
f_out1=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu1_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")
f_out2=file("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu2_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")

#f_out=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")
#f_out1=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu1_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")
#f_out2=file("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu2_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt","w")

#y1 = numpy.linspace(0,0.7,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  
#y1 = numpy.linspace(0,0.062,n)    #linspace creates array of x values, y_max=0.062 for c1=50,c0=25 
#x = numpy.empty(len(y1),float)
#y = numpy.empty(len(y1),float)

y1 = numpy.linspace(20.0/n,20.0,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  

for i in range(len(y1)):    # function len() on array gives no. of rows of array
    yy = y1[i]
    if ((1.0/4)+(yy/math.tan(c1*yy))-yy**2)<0.0:
        yy1 = y1[i-1]
        yy2 = y1[i] 
#            r = yy+ 2.0/n
#            yy_avg = (y1[i]+y1[i-1])/2.0
#            if ((1.0/4)+(yy_avg/math.tan(c1*yy_avg))-yy_avg**2)<=0.0:
        while True:
            yy3 = (yy1+yy2)/2.0
            if ((1.0/4)+(yy3/math.tan(c1*yy3))-yy3**2)<0.0:
                yy2 = yy3
                    
            if ((1.0/4)+(yy3/math.tan(c1*yy3))-yy3**2)>=0.0:
                yy1 = yy3
                
            r=yy2
            error = yy1-yy2
            if abs(error)<=(1e-15):
                break
        break


for i in range(1,18000,2):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    y = r*math.sin(math.radians(i/100.0))
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out1.write(str(x)+" "+str(y)+'\n')
    
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i>9000:
        x = -cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative    
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out1.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
#y2 = numpy.linspace(0.7,0,n)
 
#y2 = numpy.linspace(0.062,0,n)     
for i in range(17999,0,-2):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    y = -r*math.sin(math.radians(i/100.0))
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out2.write(str(x)+" "+str(y)+'\n')
    
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i>9000:
        x = -cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative    
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out2.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
#y3 = numpy.linspace(0,-0.7,n)
#y3 = numpy.linspace(0,-0.062,n)      
#for i in range(len(y3)):    
#    y = y3[i]
#    x = (-1)*cmath.sqrt((1.0/4)+(y3[i]/math.tan(c1*y3[i]))-y3[i]**2).real
#    f_out.write(str(x)+" "+str(y)+'\n')
#    f_out2.write(str(x)+" "+str(y)+'\n')
##f_out.close()    # () at the end is necessary to close the file
##f_out2.close()    # () at the end is necessary to close the file
#y4 = numpy.linspace(-0.7,0,n)
##y4 = numpy.linspace(-0.062,0,n)      
#for i in range(len(y4)):    
#    y = y4[i]
#    x = cmath.sqrt((1.0/4)+(y4[i]/math.tan(c1*y4[i]))-y4[i]**2).real
#    f_out.write(str(x)+" "+str(y)+'\n')
#    f_out2.write(str(x)+" "+str(y)+'\n')
f_out.close()    # () at the end is necessary to close the file
f_out1.close()    # () at the end is necessary to close the file
f_out2.close()    # () at the end is necessary to close the file
#print len(y1)+len(y2)+len(y3)+len(y4)

data = numpy.loadtxt("input_output_files/V_eff/quad_pot/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)
#data = numpy.loadtxt("input_output_files/V_eff_delta/quad_pot/rho_"+str(rho)+"_delta/gamma_"+str(gamma)+"/contour/nu_contour_delta_c1="+str(c1_short)+"_c0="+str(c0_short)+"_9000points_iter"+str(iteration)+".txt",float)

x = data[:,0]
y = data[:,1]
matplotlib.pylab.plot(x,y)
matplotlib.pylab.ylim(-1.25,1.25)    
