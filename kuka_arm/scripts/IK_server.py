#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []

        for x in xrange(0, len(req.poses)):

            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Define DH param symbols
            q1,q2,q3,q4,q5,q6,q7=symbols('q1:8')
    
            
            # Joint angle symbols
            a0,a1,a2,a3,a4,a5,a6=symbols('a0:7')
            alpha0,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6=symbols('alpha0:7')
            d1,d2,d3,d4,d5,d6,d7=symbols('d1:8')
    
      
            # Modified DH params
            dh = {alpha0 : 0,     a0 : 0,                    d1 : 0.75,
                  alpha1 : -pi/2, a1 : 0.35,   q2 : q2-pi/2, d2 : 0,
                  alpha2 : 0,     a2 : 1.25,                 d3 : 0,
                  alpha3 : -pi/2, a3 : -0.054,               d4 : 1.50,
                  alpha4 : pi/2,  a4 : 0,                    d5 : 0,
                  alpha5 : -pi/2, a5 : 0,                    d6 : 0,
                  alpha6 : 0,     a6 : 0,      q7 : 0,       d7 : 0.303}
    
            
            # Define Modified DH Rotation matrix
            
            R01=Matrix([[ cos(q1),             -sin(q1),            0,          ],
                        [ sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0)],
                        [ sin(q1)*sin(alpha0), cos(q1)*sin(alpha0), cos(alpha0)]])
    
            R01=R01.subs(dh) 
    
            R12=Matrix([[ cos(q2),             -sin(q2),            0,          ],
                        [ sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1)],
                        [ sin(q2)*sin(alpha1), cos(q2)*sin(alpha1), cos(alpha1)]])
    
            R12=R12.subs(dh)
    
            R23=Matrix([[ cos(q3),             -sin(q3),            0,          ],
                        [ sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2)],
                        [ sin(q3)*sin(alpha2), cos(q3)*sin(alpha2), cos(alpha2)]])
    
            R23=R23.subs(dh)
            R03=simplify(R01*R12*R23)

            #Correction from DH to URDF
    
            #       Body fixed rotation about z-axis by 180 degrees
            Rz = np.matrix([[ cos(pi), -sin(pi), 0, 0 ],
                            [ sin(pi), cos(pi),  0, 0 ],
                            [ 0,       0,        1, 0 ],
                            [ 0,       0,        0, 1 ]])
    
            #       Body fixed rotation about y-axis by -90 degrees
            Ry = np.matrix([[ cos(-pi/2),  0, sin(-pi/2), 0 ],
                            [ 0,           1, 0,          0 ],
                            [ -sin(-pi/2), 0, cos(-pi/2), 0 ],
                            [ 0,           0, 0,          1]])
    
    
            Rcorrect=np.dot( Rz, Ry )


            # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            # Calculate joint angles using Geometric IK method

            # Calculate rotation matrix from base to gripper link based on end effector orientation values

            Rrpy=np.dot(tf.transformations.quaternion_matrix([req.poses[x].orientation.x, req.poses[x].orientation.y,req.poses[x].orientation.z, req.poses[x].orientation.w]),Rcorrect.getT())
	    Rrpy=np.delete(Rrpy,(3),axis=0)
            Rrpy=np.delete(Rrpy,(3),axis=1)

            # Calulate position of wrist center by moving along end effector z XIS

	    xc=float(px-(dh[d6]+dh[d7])*Rrpy[0,2])
	    yc=float(py-(dh[d6]+dh[d7])*Rrpy[1,2])
            zc=float(pz-(dh[d6]+dh[d7])*Rrpy[2,2])

            # Variables needed for calculations (described in write-up)
	    r=np.sqrt(xc**2+yc**2)
	    l=np.sqrt((dh[a3])**2+(dh[d4])**2)
            n=np.sqrt((zc-dh[d1])**2+(r-dh[a1])**2)

	    cosA=((dh[a2])**2+n**2-l**2)/(2*dh[a2]*n)
            cosB=(l**2+(dh[a2])**2-n**2)/(2*l*dh[a2])

            # Calculate theta1, theta2 and theta3

	    theta1=np.arctan2(yc,xc)
	    theta2=np.pi/2-np.arctan2(np.sqrt(1-cosA**2),cosA)-np.arctan2(zc-dh[d1],r-dh[a1])
            theta3=np.arctan2(dh[d4],-dh[a3])-np.arctan2(np.sqrt(1-cosB**2),cosB)

            # Calculate R03 and R36

            R03=R03.evalf(subs={q1:theta1,q2:theta2,q3:theta3})
            R36=(R03.T)*Matrix(Rrpy)

            #Calculate theta4, theta5 and theta6 using R36

            cosq5=float(R36[1,2])

            theta5=np.arccos(cosq5)

            if (theta5 == 0):
		theta4 = 0
		theta6 = np.arctan2(-float(R36[2,0]),float(R36[0,0]))
            elif (theta5 == np.pi):
		theta4 = 0
		theta6 =-np.arctan2(float(R36[2,0]),-float(R36[0,0]))
            else:
		theta4 = np.arctan2(float(R36[2,2]),-float(R36[0,2]))
		theta6 = np.arctan2(-float(R36[1,1]),float(R36[1,0]))

#            print (theta1,theta2,theta3,theta4,theta5,theta6)

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        return CalculateIKResponse(joint_trajectory_list)

def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
