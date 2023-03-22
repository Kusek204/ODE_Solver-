#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:25:03 2023

@author: riebel

"""
"""
This program is a first order differential equations solver. The user need to enter the initial x and y values, dy/dx,
step size, and the final y values to initialize the problem. Then the program provide Runge Kutta 2, Runge Kutta 4, and
Adams Bashforth Multon to approximate the solution. If the user has the exact solution, the program will print a table
that demonstrates global error; else it will return the approximate answer for the user. 
"""
import sympy as sp
import numpy as np
import math
#%%
def Initialization(): #user input to initialize the problem
    f = input("enter your dy/dx: ")
    x0 = float(input("enter your x0: "))
    y0 = float(input("enter your y0: "))
    yf = float(input("enter your yf: "))
    h = float(input("enter your step size: "))
    n = int((yf-y0)/h) #calculating how many steps are needed
    check = input("Is there an exact solution? (enter Yes ot No)")
    if check == "Yes":
        check = True
    elif check == "No":
        check = False
    return f,x0,y0,yf,h,n,check
#%%
def F(x,y,f): #calculate the approximate y value 
    return eval(f)
#%%
def Ans(x,ans): # calculate the true answer by using the solution user provided
    return eval(ans)
#%%
def true_ans(x0,n,h,ans): #calculate and store the exact solution of each step
    trueans = np.zeros(n+1,dtype=float)
    for i in range(0,len(trueans)):
        trueans[i] = Ans(x0,ans)
        x0+=h
    return trueans
#%%
def RK2(x0,y0,h,n,f): #calculate x and y by using the Runge Kutta 2 method 
    x = np.zeros(n+1,dtype=float)
    x[0] = x0
    y = np.zeros(n+1,dtype=float)
    y[0] = y0
    for i in range(1,len(y)):
        k1 = F(x[i-1],y[i-1],f)
        k2 = F(x[i-1]+h,y[i-1]+h*k1,f)
        y[i] = y[i-1] + h/2 *(k1+k2)
        x[i] = x[i-1] + h
    return x,y
#%%
def RK4(x0,y0,h,n,f): #calculate x and y by using the Runge Kutta 4 method
    x = np.zeros(n+1,dtype=float)
    x[0] = x0
    y = np.zeros(n+1,dtype=float)
    y[0] = y0
    for i in range(1,len(y)):
        k1 = F(x[i-1],y[i-1],f)
        k2 = F(x[i-1]+1/2*h,y[i-1]+1/2*h*k1,f)
        k3 = F(x[i-1]+1/2*h,y[i-1]+1/2*h*k2,f)
        k4 = F(x[i-1]+h,y[i-1]+h*k3,f)
        y[i] = y[i-1]+h*(k1+ 2*k2+ 2*k3+k4)/6
        x[i] = x[i-1]+h
    return x,y
#%%
def ABM(x,y,h,f): #calculate x and y by using the ABM method
    y_ = np.zeros(len(y),dtype=float) #this is y'
    for i in range(len(y)):
        y_[i] = F(x[i],y[i],f)
    ys = y[-1] + h/24 * (55*y_[-1]-59*y_[-2]+37*y_[-3]-9*y[-4]) #this is y*
    y_n = F(x[-1]+h,ys,f)
    y_.tolist().append(y_n)
    y_ = np.asarray(y_)
    fa = y[-1]+h/24 * (9*y[-1]+19*y[-2]-5*y[-3]+y[-4])
    y.tolist().append(fa)
    y = np.asarray(y)
    return y
#%%
def Error(y,trueans,n): #calculate and store the absolute and percent error 
    abs_err = np.zeros(n+1,dtype=float)
    percent_err = np.zeros(n+1,dtype=float)
    for i in range(0,len(y)):
        abs_err[i] = trueans[i]-y[i]
        percent_err[i] = abs_err[i]/trueans[i]*100
    return abs_err, percent_err
#%%
def Table(x,y,abs_err,percent_err,trueans): #print a table with x values, approximate y values, actual y values, absolute error, and percent error 
    from prettytable import PrettyTable
    t = PrettyTable()
    t.field_names = ["Xn","Yn","Actual Value","Abs. Error","% Rel. Error"]
    for i in range(0,len(x)):
        t.add_row([x[i],y[i],trueans[i],abs_err[i],percent_err[i]])
    print(t)
#%%
def Answer(y): #print the approximate answer 
    print("Your answer is approximately to be: %.5s" % y[-1])
#%%
def Execute():
    f,x0,y0,yf,h,n,check = Initialization()
    meth = input("What method do you wnat to use (choose one from RK2, RK4, or ABM)\n:")
    if meth == "RK2":
        x,y = RK2(x0,y0,h,n,f) #execute Runge Kutta 2 method
    elif meth == "RK4":
        x,y = RK4(x0,y0,h,n,f) #execute Runge Kutta 4 method
    elif meth == "ABM":
        x,y = RK4(x0,y0,h,n,f) #execute Runge Kutta 4 in order to do ABM
        y = ABM(x,y,h,f)
    if check == True: #if there are actual solutions
        ans = input("Enter the equation here (use '**' for exponential and math.e for e)\n:")
        trueans = true_ans(x0,n,h,ans) 
        abs_err, percent_err = Error(y,trueans,n)
        Table(x,y,abs_err,percent_err,trueans)
    elif check == False: #if there are no actual solutions
        Answer(y)