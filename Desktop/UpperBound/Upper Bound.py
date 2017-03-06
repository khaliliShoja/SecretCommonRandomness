# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 16:38:02 2017

@author: mohammad
"""
'''
import test_1
import test_2
a=test_1.jam(10,10)
print(a)
a=test_2.jam(8,8)
print(a)
'''

import gamma_xyz
import gamma_xz
import gamma_yz
import gamma_z
l_1=gamma_xyz.lyapunov_list
print(l_1)
l_2=gamma_xz.lyapunov_list
print(l_2)
l_3=gamma_yz.lyapunov_list
print(l_3)
l_5=gamma_z.lyapunov_list
print(l_5)
l_4=[]
for i, j in enumerate(l_1):
    l_4.append(l_1[i]-l_2[i]-l_3[i]+l_5[i])

print(l_4)

