print("Helo world");# coding: shift_jis
#TLEから衛星位置and直下点緯度経度の計算プログラム_2017/1/19/14:57(UTC)を仮定(1/1から18.62152778日目)(L20,L21により日時変更可)

from math import *
from keido_keisan_numpy import *
from calcu_Sat_pos_to_AZEl import *

#----------------------------------------定数----------------------------------------------
th0=100.837942 #2017/1/1_平均恒星時_国立天文台のサイトから取得(100.837942度,0.2801004382rev/day)
GE=2.975537*(10 ** 15) #gravity_constant
I=0 #ループ回数初期値
ae=6378.137#赤道半径
J2=0.0010826#地球重力ポテンシャルJ2項

#------------------------------------メインループ-------------------------------------------
#create_env()#グラフィクス準備

#TLE読み込み--------------------------------------------------------------------------------
for line2 in open("TLE201717to18.txt"):
    data = line2.split()
    if data[0]=="1":
        delutaT=(17018.62291667-float(data[3]))#設定時間と元期の差(日)
        thG_deg=th0+0.004178074622295*(18.62291667*24*60*60)#343.443570	平均恒星時計算(度)
        thG=(thG_deg%360)*pi/180.0#恒星時をラジアンに変換
        delutas=delutaT*24*60*60
        Sat_number=data[1]
        ndot=float(line2[35:43])/100000000.0
        if data[33:34]=='-':
            ndot=(-1)*ndot
    elif data[0]=="2":
        n=float(data[7])#平均運動
        n=n+ndot*delutaT#観測時刻での平均運動
        M0=float(data[6])*pi/180.0#平均近点角
        M=M0+((n*2*pi*delutaT)%(2*pi)) #nは一日当たりの回転数なので角度に変換してる.観測時刻でのM
        e=float(data[4])/10000000.0#離心率
        lomega=float(data[3])*pi/180.0#昇交点赤経
        somega=float(data[5])*pi/180.0#近地点引数
        i_do=float(data[2])#軌道傾斜角(deg)
        i=float(data[2])*pi/180.0
#ここから要素計算---------------------------------------------------------------------------       
        a=calcu_semimajoraxis(n)#軌道長半径計算
        SecularPerturbation=calcu_SecularPerturbation(J2,ae,a,e,i,n,delutaT)#近地点引数、昇交点赤経の摂動計算
        somega=somega+SecularPerturbation[0]#近地点引数に摂動加える
        lomega=lomega+SecularPerturbation[1]#昇交点赤経に摂動加える
        E=calcu_risinkintenkaku(M,e)#離心近点角計算
        v_kidoumen=calcu_kidoumenzahyou(E,a,e)#軌道面座標計算
        v_sekidou=calcu_sekidoutyokkou(lomega,somega,i,v_kidoumen)#地心赤道直交座標に変換
        v_Gtisin=calcu_Gtisintyokkou(v_sekidou,thG)#G系地心直行座標系に変換
        v_Jtisin=calcu_long(v_Gtisin)#J系地心直交座標に変換
        X=v_Jtisin[0]
        Y=v_Jtisin[1]
        Z=v_Jtisin[2]
        r=sqrt(X*X+Y*Y+Z*Z)
        keido=calcu_keido(v_Jtisin)
        ido=calcu_ido_daen(v_Jtisin,a,e)
        keido2keta=round(keido,2)
        ido2keta=round(ido,2)
        h=calcu_koudo(keido,ido,ae,Z)

#地平座標AZEl(局から見た座標)に変換------------------------------------------------------------
        uvw=Cal_ground_station_pos(thG)#局の座標計算(鹿島)
        u=uvw[0]
        v=uvw[1]
        w=uvw[2]
        xyz=Cal_TiheizahyouHenkan(v_Jtisin,uvw)#地平座標に変換
        Az=atan2(-xyz[1],xyz[0])
        El=asin(xyz[2]/sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]))
        Az=Az*180/pi
        if Az<=0:
            Az=360-abs(Az)
        El=El*180/pi

#結果表示-------------------------------------------------------------------------------------
        if r>=30000 and abs(i) <=15:#only_geo
            print(name,Sat_number,keido,ido,Az,El,h[0],thG_deg%360)
        
        #sphereI = sphere(pos=(X, Y, Z), radius=100)#グラフィックス表示
        #I+1

    else :
        name=data
print("Finish!")
