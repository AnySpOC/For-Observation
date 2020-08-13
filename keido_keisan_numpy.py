
#TLEから直下点緯度経度への計算のプログラム。関数集
from math import *
import numpy as np

def calcu_semimajoraxis(n):
    GE=2.975537*(10 ** 15)
    abunsuu=GE/float((4*(pi ** 2)*(n ** 2)))
    a=(abunsuu) ** (1.0/3.0)
    #a=42241.09773/(n **(2/3))
    return(a)

def calcu_risinkintenkaku(M,e):
    E=M
    while True :
        E2=E-(E-e*sin(E)-M)/float((1-e*cos(E)))
        if abs(E2-E) < 0.00000000001:
            break
        E=E2
    return(E)

def calcu_SecularPerturbation(J2,ae,a,e,i,n,delutaT):
    #dsomega=-(3/4.0)*J2*((ae/float(a*(1-e*e)))**2)*(1-5*cos(i)*cos(i))*n*2*pi*delutaT
    #dlomega=(-3/2.0)*J2*((ae/float(a*(1-e*e)))**2)*cos(i)*n*2*pi*delutaT

    dsomega=((3*n*J2*ae*ae/float(2*a*a*(1-e*e)*(1-e*e)))*(2-(5/2.0)*sin(i)*sin(i))*delutaT)*pi/180.0
    dlomega=(((-3*J2*ae*ae*n)/float(2*a*a*(1-e*e)*(1-e*e)))*cos(i)*delutaT)*pi/180.0

    #dsomega=((180*0.174*(2-2.5*sin(i)*sin(i)))/(pi*((a/ae)**3.5)))*delutaT*pi/180
    #dlomega=((180*0.174*cos(i))/(pi*((a/ae)**3.5)))*delutaT*pi/180
    
    SecularPerturbation=np.matrix([dsomega,dlomega])
    SecularPerturbation=SecularPerturbation.T
    return(SecularPerturbation)
    
def calcu_kidoumenzahyou(E,a,e):
    U=a*(cos(E)-e)
    V=sqrt(1-(e**2))*(a*sin(E))
    v_kidoumen=np.matrix([U,V,0])
    return(v_kidoumen)

def calcu_sekidoutyokkou(lomega,somega,i,v_kidoumen):
    M1=np.matrix([[cos(somega),-sin(somega),0],[sin(somega),cos(somega),0],[0,0,1]])
    M2=np.matrix([[1,0,0],[0,cos(i),-sin(i)],[0,sin(i),cos(i)]])
    M3=np.matrix([[cos(lomega),-sin(lomega),0],[sin(lomega),cos(lomega),0],[0,0,1]])
    v_kidoumen=v_kidoumen.T
    v_sekidou=M3*M2*M1*v_kidoumen
    return(v_sekidou)

def calcu_Gtisintyokkou(v_sekidou,thG):
    M4=np.matrix([[cos(-thG),-sin(-thG),0],[sin(-thG),cos(-thG),0],[0,0,1]])
    v_Gtisin=M4*v_sekidou
    return(v_Gtisin)

def calcu_long(v_Gtisin):
    hoseiti=np.matrix([0.136,-0.521,-0.681])
    hoseiti=hoseiti.T
    v_Jtisin=v_Gtisin+hoseiti
    return(v_Jtisin)

def calcu_keido(v_Jtisin):
    X=v_Jtisin[0]
    Y=v_Jtisin[1]
    Z=v_Jtisin[2]
    keido=(atan2(Y,X))*180.0/pi
    return(keido)

def calcu_ido(v_Jtisin):
    X=v_Jtisin[0]
    Y=v_Jtisin[1]
    Z=v_Jtisin[2]
    bunbo=sqrt(X*X+Y*Y)
    #bunsu=Z/bunbo
    ido=atan2(Z,bunbo)*180.0/pi
    return(ido)

def calcu_ido_daen(v_Jtisin,a,e):
    X=v_Jtisin[0]
    Y=v_Jtisin[1]
    Z=v_Jtisin[2]
    bunbo=sqrt(X*X+Y*Y)
    fai=0
    i=0
    while(1):
        bunb=sqrt(1+(1-e*e)*tan(fai)*tan(fai))
        fai2=atan2(Z*bunb+a*e*e*tan(fai),sqrt(X*X+Y*Y)*bunb)
        i=i+1
        subfai=abs(tan(fai)-tan(fai2))
        if(subfai<=0.000000001):
            break
        else:
            fai=fai2
    ido=fai2*180.0/pi
    return ido

def calcu_koudo(keido,ido,ae,Z):
    e2=0.006674372230614#ベッセルの楕円体の離心率
    N=ae/sqrt(1-e2*sin(ido)*sin(keido))#東西線曲率半径
    ido=ido*pi/180
    w=Z
    h=sqrt(1+tan(ido)*tan(ido))*((w/tan(ido))-((ae*(1-e2))/sqrt(1+(1-e2)*tan(ido)*tan(ido))))
    return h
