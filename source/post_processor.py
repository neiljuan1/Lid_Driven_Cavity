import os
import sys
import numpy as np
from scipy import interpolate

def contour_plot(x,y,nx,ny,u,v,pres):    
    file = os.path.join(sys.path[0],"solution.dat")
    with open(file,'w') as f:
        f.write('TITLE = \"VELOCITY & PRESSURE PROFILE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"P\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        for i in range(0,nx):
            for j in range(0,ny):
                f.write('{} {} {} {} {} \n'.format(x[i],y[j],u[nx-1-j,i],v[nx-1-j,i],pres[nx-1-j,i])) #start from bottom row!
                #f.write('{} {} {} {} {} \n'.format(x[j],y[i],phi[nx-1-j,ny-1-i],u[nx-1-j,ny-1-i],v[nx-1-j,ny-1-i]))

def X_X(x,y,nx,ny,phi,u,v): 
    file = os.path.join(sys.path[0],'X-X_distance.dat')
    with open(file,'w+') as f:
        f.write('TITLE = \"TEMPERATURE PROFILE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"T\" \"U\" \"V\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        for i in range(0,nx):
            f.write('{} {} \n'.format(np.sqrt(x[i]**2+y[i]**2),np.diag(np.flip(np.flip(phi,1)))[i]))  

def horizontal(x,y,nx,ny,u,v,pres): 
    file = os.path.join(sys.path[0],'X-X_distance.dat')
    with open(file,'w+') as f:
        f.write('TITLE = \"HORIZONTAL LINE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"P\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        midline = int(nx/2)
        for j in range(0,ny):
            f.write('{} {} {} {} \n'.format(x[j],u[midline,j],v[midline,j],pres[midline,j]))                 

def interp_cp(x,y,xf,yf,nxf,nyf,phi,u,v):
    f = interpolate.interp2d(x,y,phi,kind = 'cubic')
    g = interpolate.interp2d(x,y,u,kind = 'cubic')
    h = interpolate.interp2d(x,y,v,kind = 'cubic')
    phi_o = f(yf,xf)
    u_o = g(yf,xf)
    v_o = h(yf,xf)
    file = os.path.join(sys.path[0],"sol_int.dat")
    with open(file,'w') as f:
        f.write('TITLE = \"TEMPERATURE PROFILE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"T\" \"U\" \"V\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nxf,nyf))
        for i in range(0,nxf):
            for j in range(0,nyf):
                f.write('{} {} {} {} {}\n'.format(xf[i],yf[j],phi_o[nxf-1-i,nyf-1-j],u_o[nxf-1-i,nyf-1-j],v_o[nxf-1-i,nyf-1-j]))


