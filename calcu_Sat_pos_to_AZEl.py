
#G系地心直交座標系からAz,El計算プログラム

from math import *
import numpy as np
lo=140.662899*pi/180#局の経度
la=35.954793*pi/180#局の緯度
h=25#局の高度
a=6377.39715500#赤道半径
e2=0.006674372230614#ベッセルの楕円体の離心率
    
def Cal_ground_station_pos(thG):
    e2=0.006674372230614#ベッセルの楕円体の離心率
    N=a/sqrt(1-e2*sin(la)*sin(la))#東西線曲率半径
    #J系直行座標系の計算
    u=(N+h)*cos(la)*cos(lo)
    v=(N+h)*cos(la)*sin(lo)
    w=(N*(1-e2)+h)*sin(la)
    Juvw=np.matrix([u,v,w])
    Juvw=Juvw.T
    #G系直交座標系に変換
    hoseiti=np.matrix([-0.136,+0.521,+0.681])
    hoseiti=hoseiti.T
    Guvw=Juvw+hoseiti
    #地心赤道直行座標系に変換
    #M1=np.matrix([[cos(thG),-sin(thG),0],[sin(thG),cos(thG),0],[0,0,1]])
    uvw=Guvw#M1*Guvw
    return uvw

def Cal_TiheizahyouHenkan(v_Gtisin,uvw):
    dudvdw=v_Gtisin-uvw#平行移動
    #座標回転
    M1=np.matrix([[sin(la),0,-cos(la)],[0,1,0],[cos(la),0,sin(la)]])
    M2=np.matrix([[cos(lo),sin(lo),0],[-sin(lo),cos(lo),0],[0,0,1]])
    xyz=M1*M2*dudvdw
    return xyz
