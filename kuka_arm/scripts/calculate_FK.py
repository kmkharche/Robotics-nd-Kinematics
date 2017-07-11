#!/usr/bin/env python

import numpy as np
from sympy import symbols,pi,cos,sin,simplify,pprint
from sympy.matrices import Matrix

# Joint variables

q1,q2,q3,q4,q5,q6,q7=symbols('q1:8')
a0,a1,a2,a3,a4,a5,a6=symbols('a0:7')
alpha0,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6=symbols('alpha0:7')
d1,d2,d3,d4,d5,d6,d7=symbols('d1:8')

# DH parameter values
	
dh = {alpha0 : 0,     a0 : 0,                    d1 : 0.75,
      alpha1 : -pi/2, a1 : 0.35,   q2 : q2-pi/2, d2 : 0,
      alpha2 : 0,     a2 : 1.25,                 d3 : 0,
      alpha3 : -pi/2, a3 : -0.054,               d4 : 1.50,
      alpha4 : pi/2,  a4 : 0,                    d5 : 0,
      alpha5 : -pi/2, a5 : 0,                    d6 : 0,
      alpha6 : 0,     a6 : 0,      q7 : 0,       d7 : 0.303}


# Transform Matrices

T01=Matrix([[ cos(q1),             -sin(q1),            0,            a0],
            [ sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0)*d1],
            [ sin(q1)*sin(alpha0), cos(q1)*sin(alpha0), cos(alpha0),  cos(alpha0)*d1],
            [ 0,                   0,                   0,            1]])

T01=T01.subs(dh)

print ("T01 =")
pprint (T01)
print('\n')
   
T12=Matrix([[ cos(q2),             -sin(q2),            0,            a1],
            [ sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1)*d2],
            [ sin(q2)*sin(alpha1), cos(q2)*sin(alpha1), cos(alpha1),  cos(alpha1)*d2],
            [ 0,                   0,                   0,            1]]) 

T12=T12.subs(dh)  

print ("T12 =")
pprint (T12)
print('\n')

T23=Matrix([[ cos(q3),             -sin(q3),            0,            a2],
            [ sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2)*d3],
            [ sin(q3)*sin(alpha2), cos(q3)*sin(alpha2), cos(alpha2),  cos(alpha2)*d3],
            [ 0,                   0,                   0,            1]]) 

T23=T23.subs(dh)

print ("T23 =")
pprint (T23)
print('\n')
  
T34=Matrix([[ cos(q4),             -sin(q4),            0,            a3],
            [ sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3)*d4],
            [ sin(q4)*sin(alpha3), cos(q4)*sin(alpha3), cos(alpha3),  cos(alpha3)*d4],
            [ 0,                   0,                   0,            1]])  

T34=T34.subs(dh) 

print ("T34 =")
pprint (T34)
print('\n')

T45=Matrix([[cos(q5),             -sin(q5),            0,            a4],
            [sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4)*d5],
            [sin(q5)*sin(alpha4), cos(q5)*sin(alpha4), cos(alpha4),  cos(alpha4)*d5],
            [0,                   0,                   0,            1]]) 

T45=T45.subs(dh)

print ("T45 =")
pprint (T45)
print('\n')
  
T56=Matrix([[ cos(q6),             -sin(q6),            0,            a5],
            [ sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5)*d6],
            [ sin(q6)*sin(alpha5), cos(q6)*sin(alpha5), cos(alpha5),  cos(alpha5)*d6],
            [ 0,                   0,                   0,            1]])
 
T56=T56.subs(dh)

print ("T56 =")
pprint (T56)
print('\n')
  
T67=Matrix([[ cos(q7),             -sin(q7),            0,            a6],
            [ sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
            [ sin(q7)*sin(alpha6), cos(q7)*sin(alpha6), cos(alpha6),  cos(alpha6)*d7],
            [ 0,                   0,                   0,            1]]) 

T67=T67.subs(dh)

print ("T67 =")
pprint (T67)
print('\n')

# Transforms with respect to link 0

print ("T01 =")
pprint (T01)
print('\n')
   
T02=simplify ( T01 * T12)
print ("T02 =")
pprint (T02)
print('\n')

T03=simplify ( T02 * T23)
print ("T03 =")
pprint (T03)
print('\n')

T04=simplify ( T03 * T34)
print ("T04 =")
pprint (T04)
print('\n')

T05=simplify ( T04 * T45)
print ("T05 =")
pprint (T05)
print('\n')

T06=simplify ( T05 * T56)
print ("T06 =")
pprint (T06)
print('\n')

T07=simplify ( T06 * T67)
print ("T07 =")
pprint (T07)
print('\n')

# Correction from DH to URDF

#	Body fixed rotation about z-axis by 180 degrees

Rz = Matrix([[ cos(pi), -sin(pi), 0, 0 ],
             [ sin(pi), cos(pi),  0, 0 ],
             [ 0,       0,        1, 0 ],
             [ 0,       0,        0, 1 ]])

#	Body fixed rotation about y-axis by -90 degrees

Ry = Matrix([[ cos(-pi/2),  0, sin(-pi/2), 0 ],
             [ 0,           1, 0,          0 ],
             [ -sin(-pi/2), 0, cos(-pi/2), 0 ],
             [ 0,           0, 0,          1]])

Rcorrect=simplify ( Rz * Ry )

T0G=simplify(T07*Rcorrect)

print ("T0G =")
pprint (T0G)
print ('\n')

# Direct homogenous transform if we know orientaion and position

r,p,y,px,py,pz=symbols("r p y px py pz")

Rzy=Matrix([[cos(y),-sin(y),0,0],
           [sin(y),cos(y),0,0,],
           [0,0,1,0],
           [0,0,0,1]])

Ryp=Matrix([[cos(p),0,sin(p),0],
            [0,1,0,0],
            [-sin(p),0,cos(p),0],
            [0,0,0,1]])

Rxr=Matrix([[1,0,0,0],
            [0,cos(r),-sin(r),0],
            [0,sin(r),cos(r),0],
            [0,0,0,1]])

T=Matrix([[1,0,0,px],
          [0,1,0,py],
          [0,0,1,pz],
          [0,0,0,1]])

H=T*Rzy*Ryp*Rxr

print ("Base to gripper =")
pprint (H)
print('\n')

# Matrix between links 3 and 6

T36=simplify(T34*T45*T56)

print ("T36 =")
pprint (T36)
print ("\n")

# Final Transform

#T_total = simplify ( T07 * Rcorrect )

#print T_total
#print (T03.evalf(subs={q1:1.156578117790, q2:0.272831686871, q3:0.615138780493}))
#print (T06.evalf(subs={q1:1.156578117790, q2:0.272831686871, q3:0.615138780493,q4:0,q5:0,q6:0}))

#T36=simplify(T34*T45*T56)
#print(T36)
