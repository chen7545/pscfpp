#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 04:23:14 2023

@author: kexinchen
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


with open('out/c') as f:
    content = f.readlines()

content = [x.strip() for x in content]

#concentration matrix 
c = []
for line in content:
    mat =[]
    s = line.split(' ')
    #Remove iteration number
    for j in range(len(s)):
        if(isfloat(s[j])):
            mat.append(float(s[j]))
    c.append(mat)

#Remove nx and nm value
c = c[2:]

#Cound the summation of concentration at each grid differ from 1


def column(matrix, i):
    return [row[i] for row in matrix]

index = column(c, 0);
A = column(c, 1)
B = column(c, 2)

plt.plot(index,A)
plt.plot(index,B)
        