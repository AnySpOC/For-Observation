# coding: shift_jis
#AzElから直下点緯度経度関数

from math import *
from geo import *
import numpy as np

def calcu_ramuda(Az,k_ido):
    bunsi=sin(Az)*sin(Az)-1
    bunbo=sin(Az)*sin(Az)*cos(k_ido)*cos(k_ido)-1
    bunsu=sqrt(bunsi/bunbo)
    ramuda=acos(bunsu)
    return ramuda

def calcu_Sat_distance(ramudad,k_ido,ae,h):
    d=sqrt(h*h+2*ae*(ae+h)*(1-cos(k_ido)*cos(ramudad)))
    return d

def calcu_xyzh(d,El,Az):
    xyzh=VECTOR()
    xyzh[0]=d*cos(El)*cos(Az)
    xyzh[1]=d*cos(El)*sin(Az)
    xyzh[2]=d*sin(El)
    return xyzh

def calcu_kasima_zahyo(ae,k_ido,k_keido,th_k):
    xyzk=VECTOR()#地球半径はaeでなく鹿島の半径を書くべき
    xyzk[0]=ae*cos(k_ido)*cos(th_k)
    xyzk[1]=ae*cos(k_ido)*sin(th_k)
    xyzk[2]=ae*sin(k_ido)
    return xyzk

def calcu_ido(a,b,c,th_g):
    ido=atan2(c,sqrt(a*a+b*b))
    arufa=atan2(b,a)
    ramudaE=arufa-th_g
    idokeido=VECTOR()
    idokeido=(ido,ramudaE)
    return (idokeido)

#ここから方向余弦版で用いる関数
def calcu_tihei_sekidou(k_ido,th_k,Az,El):
    #方向余弦を用いて地平座標から赤道座標へ
    A=Az-pi
    l=cos(El)*cos(A)
    m=-cos(El)*sin(A)
    n=sin(El)
    lmn=np.matrix([l,m,n])
    lmn=lmn.T
    #地平座標から赤道座標へ回転
    M1=np.matrix([[cos(th_k),-sin(th_k),0],[sin(th_k),cos(th_k),0],[0,0,1]])
    M2=np.matrix([[sin(k_ido),0,cos(k_ido)],[0,1,0],[-cos(k_ido),0,sin(k_ido)]])
    M3=M1*M2
    LMN=M3*lmn
    return LMN

def calcu_sekidou(th_g,k_ido,th_k,arufa,siguma):
    #方向余弦を用いて赤道座標から地平座標へ
    L=cos(arufa)*cos(siguma)#cos(siguma)*cos(arufa)(arufa+th_g)
    M=sin(arufa)*cos(siguma)#cos(siguma)*sin(arufa)(arufa+th_g)
    N=sin(siguma)#sin(siguma)
    LMN=np.matrix([L,M,N])
    LMN=LMN.T
    M2=np.matrix([[cos(th_k),sin(th_k),0],[-sin(th_k),cos(th_k),0],[0,0,1]])
    M1=np.matrix([[sin(k_ido),0,-cos(k_ido)],[0,1,0],[cos(k_ido),0,sin(k_ido)]])
    M3=M1*M2
    lmn=M3*LMN
    l=lmn[0]
    m=lmn[1]
    n=lmn[2]
    A=atan2(-m,l)
    h=asin(n)
    if A>=0:
        Az=degrees(A+2*pi)
    else:
        Az=degrees(2*pi-abs(A))
    El=degrees(h)
    AzEl=[Az,El]
    return AzEl
#赤道座標から直下点緯度経度に変換
def calcu_Gtisintyokkou(LMN,th_g):
    M=np.matrix([[cos(th_g),sin(th_g),0],[-sin(th_g),cos(th_g),0],[0,0,1]])
    LMN2=M*LMN
    return(LMN2)

#地心視差の計算
def calcu_TisinSisaAz(R,k_ido,Az):
    deluta=42164.15
    a=(R*cos(k_ido)*sin(abs(pi-Az)))/deluta
    TisinSisa=asin(a)
    return TisinSisa
def calcu_TisinSisaEl(R,El):
    deluta=42164.15
    a=(R*sin(abs(pi/2-El)))/deluta
    TisinSisa=asin(a)
    return TisinSisa
