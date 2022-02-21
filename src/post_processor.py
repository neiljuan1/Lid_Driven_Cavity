import os
import sys
import numpy as np
from scipy import interpolate

def contour_plot(x,y,nx,ny,u,v,pres):    
    data_dir = os.path.abspath(__file__ + "/../../data")
    file = os.path.join(data_dir,"solution.dat")
    with open(file,'w') as f:
        f.write('TITLE = \"VELOCITY & PRESSURE PROFILE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"P\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        for i in range(0,nx):
            for j in range(0,ny):
                f.write('{} {} {} {} {} \n'.format(x[i],y[j],u[nx-1-j,i],v[nx-1-j,i],pres[nx-1-j,i])) #start from bottom row!
                #f.write('{} {} {} {} {} \n'.format(x[j],y[i],phi[nx-1-j,ny-1-i],u[nx-1-j,ny-1-i],v[nx-1-j,ny-1-i]))

def X_X(x,y,nx,ny,phi,u,v): 
    data_dir = os.path.abspath(__file__ + "/../../data")
    file = os.path.join(data_dir,'X-X_distance.dat')
    with open(file,'w+') as f:
        f.write('TITLE = \"TEMPERATURE PROFILE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"T\" \"U\" \"V\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        for i in range(0,nx):
            f.write('{} {} \n'.format(np.sqrt(x[i]**2+y[i]**2),np.diag(np.flip(np.flip(phi,1)))[i]))  

def horizontal(x,y,nx,ny,u,v,pres): 
    data_dir = os.path.abspath(__file__ + "/../../data")
    file = os.path.join(data_dir,'X-X_distance.dat')
    with open(file,'w+') as f:
        f.write('TITLE = \"HORIZONTAL LINE\"\n')
        #f.write('VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"P\" \n')
        f.write('VARIABLES = \"X\" \"U\" \"V\" \"P\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        midline = int(nx/2)
        for j in range(0,ny):
            f.write('{} {} {} {} \n'.format(x[j],u[midline,j],v[midline,j],pres[midline,j]))              

def vertical(x,y,nx,ny,u,v,pres): 
    data_dir = os.path.abspath(__file__ + "/../../data")
    file = os.path.join(data_dir,'Y-Y_distance.dat')
    with open(file,'w+') as f:
        f.write('TITLE = \"VERTICAL LINE\"\n')
        #f.write('VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"P\" \n')
        f.write('VARIABLES = \"Y\" \"U\" \"V\" \"P\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nx,ny))
        midline = int(nx/2)
        for j in range(0,ny):
            f.write('{} {} {} {} \n'.format(y[j],u[j,midline],v[j,midline],pres[j,midline]))     

def interp_cp(x,y,xf,yf,nxf,nyf,phi,u,v):
    data_dir = os.path.abspath(__file__ + "/../../data")
    f = interpolate.interp2d(x,y,phi,kind = 'cubic')
    g = interpolate.interp2d(x,y,u,kind = 'cubic')
    h = interpolate.interp2d(x,y,v,kind = 'cubic')
    phi_o = f(yf,xf)
    u_o = g(yf,xf)
    v_o = h(yf,xf)
    file = os.path.join(data_dir,"sol_int.dat")
    with open(file,'w') as f:
        f.write('TITLE = \"TEMPERATURE PROFILE\"\n')
        f.write('VARIABLES = \"X\" \"Y\" \"T\" \"U\" \"V\" \n')
        f.write('ZONE, I = {}, J = {}\n'.format(nxf,nyf))
        for i in range(0,nxf):
            for j in range(0,nyf):
                f.write('{} {} {} {} {}\n'.format(xf[i],yf[j],phi_o[nxf-1-i,nyf-1-j],u_o[nxf-1-i,nyf-1-j],v_o[nxf-1-i,nyf-1-j]))

def write_to_file():
    data_dir = os.path.abspath(__file__ + "/../../data")
    file = os.path.join(data_dir, "text.dat")
    with open(file, "w") as f:
        f.write("test")

write_to_file()
print(os.path.dirname(__file__))
print(os.path.abspath(__file__))