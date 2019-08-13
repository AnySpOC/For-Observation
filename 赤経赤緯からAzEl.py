# coding: shift_jis
#赤経、赤緯からAzElへの変換メインプログラム----静止軌道のみ-----

from math import *
from geo import *
from AzEl_Tyokkaten_calcu import *
import numpy as np

#定数
k_ido=radians(35.954793) #鹿島緯度rad
k_keido=radians(140.662899) #鹿島経度rad
th_g0=118.579594#1/19日00:00における平均恒星時国立天文台webより取得
ae=6378.140 #地球半径赤道半径を仮定km
h=35786 #高度静止軌道を仮定km
R=6378.15#地球半径_静止衛星設計入門より

print "観測[時]を入力してください(UTC)"
T_h=input(">>>")
print "観測[分]を入力してください(UTC)"
T_m=input(">>>")
print "観測[秒]を入力してください(UTC)"
T_s=input(">>>")
print "赤経を入力してください(deg)"
arufa=radians(input(">>>"))
print "赤緯を入力してください(deg)"
siguma=radians(input(">>>"))

subT=T_h*60*60+T_m*60+T_s#恒星時の基準時間との時間差s
th_g=th_g0+subT*0.004178074622295#恒星時計算
th_g=radians(th_g)#恒星時
th_k=th_g+k_keido#鹿島恒星時


#ここからmain関数
#---------------------方向余弦計算によるAzEL計算-------------------------
AzEl=calcu_sekidou(th_g,k_ido,th_k,arufa,siguma)
Az=AzEl[0]
El=AzEl[1]
TisinSisaAz=calcu_TisinSisaAz(R,k_ido,Az)
TisinSisaEl=calcu_TisinSisaEl(R,El)
print Az-TisinSisaAz,El-TisinSisaEl
