# Normalization: q0 = qp, m0 = mp, L0 = 1m, B0 = 1 T, t0 = /omega_c^-1 = m0/(q0*B0), v0 = L0/t0, E0 = m0*v0/(q0 * t0);

# About the gui
# The first entry is for choosing coordinates which can be cartesian or cylindrical now
# The second one is for getting a fields expression ,which should be a python command that give fields a value
# fields[0:1, :, 0:2]:the first position parameter: 0->E, 1->B;The second one is particles num;The last one is coordinates (r,z,theta) in cylindrical and (x,y,z) in cartesian
# It can be like this:'fields[0,:,1]=1;fields[1,:,2]=1' or it can be a function of position(x) like:'fields[0,:,1]=1;fields[1,:,2]=1/x[:,0]' in which x is particles' position the first parameter is particles num; the second one is coordinates (r,z,theta) in cylindrical and (x,y,z) in cartesian

from sgptc import *
from numpy import *
from mayavi import mlab
import tkinter as tk
# import _thread

x = array([[100., 0., 0.]])

def plotm(func = None, step = 10000, color = (1,0,0)):
    global testp
    global ptrack
    ptrack = []
    if func == None:
        func = testp.leap_frog

    # draw
    for i in range(step):
        ptrack.append(testp.x.copy().flatten())
        func()

    point = array(ptrack)
    if testp.cd.coordname == 'cylindrical':
        pos = zeros(point.shape)
        pos[:, 0] = point[:, 0]*cos(point[:,2])
        pos[:, 1] = point[:, 0]*sin(point[:,2])
        pos[:, 2] = point[:, 1]
    else:
        pos = point
    interv = int(step / 10000)
    figure = mlab.gcf()
    ptd = mlab.points3d(pos[0::interv, 0], pos[0::interv, 1], pos[0::interv, 2], mode = 'point', color = color)
    ptda = ptd.glyph.glyph_source.glyph_source.output.points.to_array()
    def picker_callback(picker):
        if picker.actor in ptd.actor.actors:
            point_id = int(picker.point_id / ptda.shape[0])
            if point_id != -1:
                print(pos[:, 0][point_id], pos[:, 1][point_id], pos[:, 2][point_id])
        else:
            print(picker.pick_position)
    picker = figure.on_mouse_pick(picker_callback)
    picker.tolerance = 0.001
    mlab.xlabel('x')
    mlab.outline()


    
# Curl B drift; E =0, Bphi =1,
def cb():
    global testp
    x = array([[1., 0., 0.]])
    testp = particles(dt = 1e-2, x = x, fiefun = (lambda x: array([rollaxis(array([zeros(x.shape[0]), zeros(x.shape[0]), zeros(x.shape[0])]), 0, 2), rollaxis(array([zeros(x.shape[0]), zeros(x.shape[0]), ones(x.shape[0])]), 0, 2)])), v = array([[0.1, 0.0, 0.01]]), m = 1, q = 1, coordname = t1.get())
    plotm(func = testp.boris, color = (0,1,0))
    mlab.show()



# E cross B drift; E in z(y) B in phi(z)
def eb():
    global testp
    x = array([[1., 0., 0.]])
    testp = particles(dt = 1e-2, x = x, fiefun = (lambda x: array([rollaxis(array([zeros(x.shape[0]), ones(x.shape[0]), zeros(x.shape[0])]), 0, 2), rollaxis(array([zeros(x.shape[0]), zeros(x.shape[0]), ones(x.shape[0])]), 0, 2)])), v = array([[0., 0.0, 0.]]), m = 1, q = 1, coordname = t1.get())
    plotm(func = testp.boris, color = (0,1,0))

    mlab.show()


#Grad B corss B drift; E=0, Bphi = R
def gb():
    global testp
    x = array([[1., 0., 0.]])
    testp = particles(dt = 1e-2, x = x, fiefun = (lambda x: array([rollaxis(array([zeros(x.shape[0]), zeros(x.shape[0]), zeros(x.shape[0])]), 0, 2), rollaxis(array([zeros(x.shape[0]), zeros(x.shape[0]), x[:, 0]]), 0, 2)])), v = array([[0.1, 0.01, 0.0]]), m = 1, q = 1, coordname = t1.get())
    plotm(func = testp.boris, color = (0,1,0))

    mlab.show()
# testp = particles(dt = 1e-11, fiefun = (lambda x: array([rollaxis(array([zeros(x.shape[0]),  ones(x.shape[0]), zeros(x.shape[0])]), 0, 2), rollaxis(array([ones(x.shape[0]), zeros(x.shape[0]), zeros(x.shape[0])]), 0, 2)])), v = array([[0., 0., 0.0]]))


def give_field_func():
    global testp
    x = array([[1, 0., 0.]])
    fsp = list(x.shape)
    fsp.insert(0, 2)
    fields = zeros(fsp)
    def field_func(x):
        exec(code)
        return fields
    try:
        exec(t2.get())
    except:
        print('I cannot read your function')
        exit()
    else:
        code = t2.get()
        testp = particles(dt = 1e-2, x = x, fiefun = field_func, v = array([[0.1, 0.0, 0.01]]), m = 1, q = 1, coordname = t1.get())
        plotm(func = testp.rk4, color = (0,1,0), step = 10000)
        mlab.show()
        return
    
root = tk.Tk()
root.title('sgptc')
b1 = tk.Button(root, activebackground = 'grey', activeforeground = 'green', text = 'E cross B drift', width = 20, command = eb )
b2 = tk.Button(root, activebackground = 'grey', activeforeground = 'green', text = 'Curl B drift', width = 20, command = cb )
b3 = tk.Button(root, activebackground = 'grey', activeforeground = 'green', text = 'Grad B cross B drift', width = 20, command = gb )
b4 = tk.Button(root, activebackground = 'grey', activeforeground = 'green', text = 'load field_func and run', width = 20, command = give_field_func)


# L1 = tk.Label(root, text= 'coordinates')
# L1.pack(side = tk.LEFT)
t1 = tk.StringVar()
t1.set('cylindrical')
e1 = tk.Entry(root, textvariable = t1, justify = 'center', bd=5)
e1.pack()


t2 = tk.StringVar()
t2.set('fields[0,:,0]=1;fields[1,:,2]=1')
e2 = tk.Entry(root, textvariable = t2, justify = 'center', bd = 5)
e2.pack()

b1.pack()
b2.pack()
b3.pack()
b4.pack()
root.mainloop()


# #animation
# @mlab.animate(delay = 10)
# def anim():
#     ani = mlab.points3d(x[:, 0] * cos(x[:, 2]), x[:, 0] * sin(x[:, 2]), x[:, 1], mode = 'sphere', line_width = 11.0)
#     mlab.axes()
#     for i in range(step):
#         ptrack.append(testp.x.copy().flatten())
#         testp.rk4()
#         point = testp.x.copy().flatten()
#         ani.mlab_source.reset(x = point[:, 0] * cos(point[:, 2]), y = point[:, 0] * sin(point[:, 2]), z = point[:, 1])
#         # mlab.points3d(point[0] * cos(point[2]), point[0] * sin(point[2]), point[1], mode = 'sphere', line_width = 11.0)
#         print(i)
#         yield
# anim()
# point = array(ptrack)
# pos = zeros(point.shape)
# pos[:, 0] = point[:, 0]*cos(point[:,2])
# pos[:, 1] = point[:, 0]*sin(point[:,2])
# pos[:, 2] = point[:, 1]
# interv = int(step / 10000)
# mlab.points3d(pos[0::interv, 0], pos[0::interv, 1], pos[0::interv, 2], mode = 'point')
# mlab.show()






# # mesh the psi surface
# dataList=[]
# with open('./psi', 'r') as infile:
#     for line in infile:
#         dataList.append(line)

# zmap = map((lambda x: float(x)), dataList[-1].split())
# rmap = map((lambda x: float(x)), dataList[-2].split())
# r=[]; z=[]
# for i, j in zip(rmap,zmap):
#     r.append(i)
#     z.append(j)

# theta = linspace(0, 2*pi, 100)
# r = array(r)
# z = array(z)
# z = z.reshape(z.size, 1)
# r = r.reshape(r.size, 1)
# theta = theta.reshape(1, theta.size)
# x = matmul(r, cos(theta))
# y = matmul(r, sin(theta))
# z = matmul(z, ones((1, theta.size)))
# mlab.mesh(x, y, z, opacity = 0.1)
# mlab.show()

