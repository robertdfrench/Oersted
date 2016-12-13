import os
import numpy as np
import matplotlib.pyplot as plt

imgtypes = ['.pdf','.png']

def plot_sketch(fname):
    print('Sketch: ' +fname)
    data = np.genfromtxt(fname, dtype=float, delimiter=',')
    plt.plot(data[:,0], data[:,1])
    plt.grid(True)
    plt.axes().set_aspect('equal','datalim')
    for ext in imgtypes:
        plt.savefig(fname+ext)
    plt.clf()

def plot_mesh(fname):
    print('Mesh: '+fname)
    data = np.genfromtxt(fname, dtype=float, delimiter=',')
    x = data[:,[0,1,2,0]].transpose()
    y = data[:,[3,4,5,3]].transpose()
    plt.plot(x,y,'b')
    plt.grid(True)
    plt.axes().set_aspect('equal','datalim')
    for ext in imgtypes:
        plt.savefig(fname+ext)
    plt.clf()

root = os.getcwd()
root += '/build/test/output/'

for path, _, files in os.walk(root):
    for name in files:
        _, fext = os.path.splitext(name)
        fpath = os.path.join(path, name)
        if fext == '.oesk':
            plot_sketch(fpath)
        elif fext == '.oeme':
            plot_mesh(fpath)