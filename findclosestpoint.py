# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:33:18 2021
#script from https://codereview.stackexchange.com/questions/28207/finding-the-closest-point-to-a-list-of-points
used to find closest point location within a list of points
@author: Tonian
"""
import numpy as np
def closestpoint(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

#To make this script callable as a function add the line below and do a "practice run" code below   
x = [1,2,3,4,5,6,8]
y = [1,2,3,4,5,6,8]
a = np.stack((x,y))
a = a.T
pt = [3.12, 3.12]
if __name__ == "__main__":
    w = closestpoint(pt,a)
    print(w)
    