# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd

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
    s = s[1:]
    for j in range(len(s)):
        if(isfloat(s[j])):
            mat.append(float(s[j]))
    c.append(mat)

#Remove nx and nm value
c = c[2:]

#Cound the summation of concentration at each grid differ from 1
diff = 0
for i in range(len(c)):
    sum_c = 0
    for j in range(len(c[i])):
        #sum the the concentation at grid i
        sum_c += c[i][j]
    diff += abs(sum_c -1)

print("The sum of inaccuracy of concentration on each grid is ")
print(diff)
        



