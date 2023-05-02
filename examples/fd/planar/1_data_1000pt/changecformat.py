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
data = c[2:]

output_file = "output.txt"

# Write the 2D numpy array to the text file with a space before each line
with open(output_file, "w") as f:
    for row in data:
        formatted_row = " ".join(map(lambda x: f"{x:.11e}", row))
        f.write(f" {formatted_row}\n")

        



