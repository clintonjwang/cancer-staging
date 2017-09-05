# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 13:58:12 2016

@author: jyoung
"""

import collections
import matplotlib.pyplot as plt
import numpy as np
import SimpleITK as sitk


class Seeds(object):
    """For a stack of images, capture and save coordinates from mouse clicks.
    Either an arrow or a circle can be displayed upon each click."""

    def __init__(self, ax, X, shape):
        self.ax = ax
        self.X = X
        self.shape = shape
        self.slices, numrows, numcols = X.shape
        self.idx = self.slices//2
        self.coords = list()
        self.slice2shapes = collections.defaultdict(list)

        self.im = ax.imshow(self.X[self.idx,:,:], cmap=plt.cm.Greys_r)
        self.update()

    def onscroll(self, event):
        try:
            for rec in self.slice2shapes[self.idx]:
                rec.remove()
        except:
            pass
        if event.button == 'up':
            self.idx = np.clip(self.idx+1, 0, self.slices-1)
        else:
            self.idx = np.clip(self.idx-1, 0, self.slices-1)
        if self.slice2shapes[self.idx] != []:
            for rec in self.slice2shapes[self.idx]:
                plt.gca().add_artist(rec)
        self.update()

    def onclick(self, event):
        ix, iy = int(round(event.xdata)), int(round(event.ydata))
        self.coords.append((ix, iy, int(self.idx)))
        if type(self.shape) is int:
            self.slice2shapes[self.idx].append(plt.Circle((ix, iy), self.shape, 
                color='r'))
        else:
            self.slice2shapes[self.idx].append(plt.annotate('', xy=(ix,iy-4), 
                xytext=(ix,iy+10), arrowprops=dict(arrowstyle='->', color='r')))
        self.im.axes.figure.canvas.draw()

    def update(self):
        self.im.set_data(self.X[self.idx,:,:])
        self.ax.set_title('slice %s' %self.idx)
        self.im.axes.figure.canvas.draw()


def get_seeds(imgstack, shape='arrow'):
    """Show stack of images with mouse scrolling. If 'shape' is an integer, a
    circle with a radius of that size will be drawn at each click. Otherwise,
    an arrow will be displayed by default."""
    X = sitk.GetArrayFromImage(imgstack)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    seedstore = Seeds(ax, X, shape)
    fig.canvas.mpl_connect('scroll_event', seedstore.onscroll)
    fig.canvas.mpl_connect('button_press_event', seedstore.onclick)
    plt.show()

    return seedstore.coords
