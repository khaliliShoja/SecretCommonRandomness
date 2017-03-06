# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 00:26:57 2017

@author: mohammad
"""

# -*- coding: utf-8 -*-
"""


@author: mohammad
"""
import numpy as np
from numpy import linalg as LA
P=[]
a_t=0
b_t=0
for i in range(0,19):
    a_t=a_t+.05
    b_t=b_t+.05
    u=[[1-a_t, a_t],[b_t,1-b_t]]
    P.append(u)
  
lyapunov_list=[]
for i in P:
    
        
    P_transmission=i
    N=200
    delta=(2*N)
    
    M=np.zeros(shape=(2*(N),2*(N)))
    #Defining Transmission Probability
    alpha=P_transmission[0][1]
    beta=P_transmission[1][0]
    
    #Defining Emission Probability for x and y
    P_Emission_x=[[.6, .4], [.1, .9]]
    P_Emission_y=[[.75, .25], [.35, .65]]
    P_Emission_z=[[.7, .3], [.4, .6]]  
    #print(P_Emission_y)
    # Defining Q_0 for sequence XY
    def Q_z_0(w,x):
        a=0
        if(x==1):
            a=1
        elif(x==0):
            a=0
        else:
        
            z=(x/(1-x))*(alpha * w + (1-beta)* (1-w))/((1-alpha)*w + beta *(1-w))
        
            a_1=((P_Emission_z[0][0])/(P_Emission_z[1][0]))
            if(a_1<= z):
                a=a+P_Emission_z[0][0]
        
            a_3=((P_Emission_z[0][1])/(P_Emission_z[1][1]))
            if(a_3<= z):
                a=a+ P_Emission_z[0][1]
    
    
        
        return a    
    
    # Defining Q_1 for sequence XY
    def Q_z_1(w,x):
        a=0
        if(x==1):
            a=1
        elif(x==0):
            a=0
        else:
            z=(x/(1-x))*(alpha * w + (1-beta)* (1-w))/((1-alpha)*w + beta *(1-w))
      
            
            a_1=((P_Emission_z[0][0])/(P_Emission_z[1][0]))
            if(a_1<= z):
                a=a+P_Emission_z[1][0]
        
    
            a_3=((P_Emission_z[0][1])/(P_Emission_z[1][1]))
            if(a_3<= z):
                a=a+ P_Emission_z[1][1]
    
    
    
        return a    
    
    # Defining Derivatives
    def der_Q_0(w,x):
        #print(delta)
        #der=(Q_xy_0(w,x+delta)-Q_xy_0(w,x-delta))/(2*delta)
        der=(Q_z_0(w,x+(2/delta))-Q_z_0(w,x))/(2/delta)     
        return der
        
    def der_Q_1(w,x):
        
        #der=(Q_xy_1(w,x+delta)-Q_xy_1(w,x-delta))/(2*delta)
        der=(Q_z_1(w,x+(2/delta))-Q_z_1(w,x))/(2/delta)
        return der    
    
    # Defining matrix M
    for x_i in range(0,N):
        for w_i in range (0,N):
            if(w_i==0):
                M[x_i][w_i]=(1-alpha)*(1/N)*der_Q_0((w_i/N),(x_i/N))
                #print((1-alpha)*(1/N)*der_Q_0((w_i/N),(x_i/N)))
                #print("ali", x_i/N)
            else:
                M[x_i][w_i]=(1-alpha)*(1/(2*N))*(der_Q_0((w_i/N),(x_i/N))+der_Q_0(((w_i+1)/N),(x_i/N)))
                #print((1-alpha)*(1/(2*N))*(der_Q_0((w_i/N),(x_i/N))+der_Q_0(((w_i+1)/N),(x_i/N))))
        for w_i in range (0,N):
            if(w_i==0):
                M[x_i][w_i+(N)]=(beta)*(1/N)*der_Q_0((w_i/N),(x_i/N))
                #print((beta)*(1/N)*der_Q_0((w_i/N),(x_i/N)))
            else:
                M[x_i][w_i+(N)]=(beta)*(1/(2*N))*(der_Q_0((w_i/N),(x_i/N))+der_Q_0(((w_i+1)/N),(x_i/N)))
                #print((beta)*(1/(2*N))*(der_Q_0((w_i/N),(x_i/N))+der_Q_0(((w_i+1)/N),(x_i/N))))
        
    
    for x_i in range(0,N):
        for w_i in range (0,N):
            if(w_i==0):
                M[x_i+(N)][w_i]=(alpha)*(1/N)*der_Q_1((w_i/N),(x_i/N))
                #print((alpha)*(1/N)*der_Q_1((w_i/N),(x_i/N)))
            else:
                M[x_i+(N)][w_i]=(alpha)*(1/(2*N))*(der_Q_1((w_i/N),(x_i/N))+der_Q_1(((w_i+1)/N),(x_i/N)))
                #print((alpha)*(1/(2*N))*(der_Q_1((w_i/N),(x_i/N))+der_Q_1(((w_i+1)/N),(x_i/N))))
        for w_i in range (0,N):
            if(w_i==0):
                M[x_i+(N)][w_i+(N)]=(1-beta)*(1/N)*der_Q_1((w_i/N),(x_i/N))
                #print((1-beta)*(1/N)*der_Q_1((w_i/N),(x_i/N)))
            else:
                M[x_i+(N)][w_i+(N)]=(1-beta)*(1/(2*N))*(der_Q_1((w_i/N),(x_i/N))+der_Q_1(((w_i+1)/N),(x_i/N)))
                #print((1-beta)*(1/(2*N))*(der_Q_1((w_i/N),(x_i/N))+der_Q_1(((w_i+1)/N),(x_i/N))))
           
    # End of Defining Matrix M
                
                
    #print("########################")            
    #print(Q_xy_1(1/N, 1))            
    #print(M)            
    #print("salam")
    #val, vector= LA.eig(M)
    #print(val)
    # To find the corresponding eigen vector wih one
    def Calculator_Eigen_1(M,N):
        val, vector= LA.eig(M)
        index_get="error"
        for index, value in enumerate(val):
            
            if(abs(value.real-1)<=.0001 and value.imag==0):
                index_get=index
        vector_get=np.zeros(shape=(2*N,1))
        for i in range(0,2*N):
            #print("PPPPPPPPPPPPP",P_transmission)
            #print(index_get)
            vector_get[i]=vector[i][index_get].real
        #print("aliiiiiiii", index_get)
        return vector_get 
    
    vector_get=Calculator_Eigen_1(M,N)
    #print(vector_get) 
    #val, vector= LA.eig(M)
    
    ##print("eigen value 1 Start")
    ##print(vector_get)
    ##print("eigen value 1 End")
    #vector_get_new=2*vector_get
    #print(vector_get_new)
    
    def Normaliz_Eigen_To_one(vector_get):
        sum_column=0
        for i in vector_get:
            sum_column=sum_column+i
        #print("ssssssssuuuuuuuummmmm", sum_column)
        vector_get_new=(1/sum_column)*vector_get
        return vector_get_new
    
       
      
    
    def m_0(vector_get):
        m_0=np.zeros(shape=(N,1))
        for i in range(0,N):
            m_0[i]=vector_get[i]
        return m_0
    
        
    def m_1(vector_get):
        m_1=np.zeros(shape=(N,1))
        for i in range(0,N):
            m_1[i]=vector_get[i+N]
        return m_1
    
    
    m_0=m_0(vector_get)
    m_1=m_1(vector_get)
    
    m_0=Normaliz_Eigen_To_one(m_0)
    m_1=Normaliz_Eigen_To_one(m_1)
    
    ##print(m_0)
    ##print(m_1)
    
    ''' Test Of normalization
    print(Normaliz_Eigen_To_one(vector_get))
    suuuum=0
    for i in Normaliz_Eigen_To_one(vector_get):
        suuuum=suuuum+i
    print(suuuum) 
    ''' #End Of Testing Normalization
    
       
    #print(val)
    '''
    # For writing the eigen values and eigen vectors in a file
    for i in val:
        g.writelines(str(i))
    for i in vector:
        f.writelines(str(i))
    #print("EiGGGGGGGGGGGGGen VVVVVVvector")
    
        
    #print(vector)    
    f.close()
    g.close()
    # End of writing the eigen values and eigen vectors in a file
    '''
    # defining G_0
    def G_0(w_i):
        S_G_0=0
        for i in range (0,2):
            
                i_1=((1-alpha)*w_i + beta*(1-w_i)) * P_Emission_z[0][i]
                i_2=(alpha * w_i + (1-beta)*(1-w_i)) * P_Emission_z[1][i]
                S_G_0 = S_G_0 + (np.log10(i_1+i_2) * P_Emission_z[0][i])
                
        return S_G_0
        
    # defining G_1
    def G_1(w_i):
        S_G_1=0
        for i in range (0,2):
                i_1=((1-alpha)*w_i + beta*(1-w_i)) * P_Emission_y[0][i]
                i_2=(alpha * w_i + (1-beta)*(1-w_i)) *P_Emission_y[1][i]
                S_G_1 = S_G_1 + (np.log10(i_1+i_2) *P_Emission_y[1][i])
                
        return S_G_1
    
        
    # Calculating Lyapanov
    
    def lyapanov_1(m_0, m_1):
        i_1=0
        for x_i in range(0,N):
            i_1=i_1 + ((((1-alpha) * G_0(x_i/N)) + (alpha * G_1(x_i/N))) * m_0[x_i])
            
        i_2=0
        for x_i in range(0,N):
            i_2=i_2 + ((((beta) * G_0(x_i/N)) + ((1-beta) * G_1(x_i/N))) * m_1[x_i]) 
        
        return(i_1 + i_2)
        
        
        
    lyap=lyapanov_1(m_0, m_1)
    lyapunov_list.append(lyap)
'''    
for i in lyapunov_list:
    print(i)
print(lyapunov_list[0])
'''