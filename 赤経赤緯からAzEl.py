# coding: shift_jis
#�Ԍo�A�Ԉ܂���AzEl�ւ̕ϊ����C���v���O����----�Î~�O���̂�-----

from math import *
from geo import *
from AzEl_Tyokkaten_calcu import *
import numpy as np

#�萔
k_ido=radians(35.954793) #�����ܓxrad
k_keido=radians(140.662899) #�����o�xrad
th_g0=118.579594#1/19��00:00�ɂ����镽�ύP���������V����web���擾
ae=6378.140 #�n�����a�ԓ����a������km
h=35786 #���x�Î~�O��������km
R=6378.15#�n�����a_�Î~�q���݌v������

print "�ϑ�[��]����͂��Ă�������(UTC)"
T_h=input(">>>")
print "�ϑ�[��]����͂��Ă�������(UTC)"
T_m=input(">>>")
print "�ϑ�[�b]����͂��Ă�������(UTC)"
T_s=input(">>>")
print "�Ԍo����͂��Ă�������(deg)"
arufa=radians(input(">>>"))
print "�Ԉ܂���͂��Ă�������(deg)"
siguma=radians(input(">>>"))

subT=T_h*60*60+T_m*60+T_s#�P�����̊���ԂƂ̎��ԍ�s
th_g=th_g0+subT*0.004178074622295#�P�����v�Z
th_g=radians(th_g)#�P����
th_k=th_g+k_keido#�����P����


#��������main�֐�
#---------------------�����]���v�Z�ɂ��AzEL�v�Z-------------------------
AzEl=calcu_sekidou(th_g,k_ido,th_k,arufa,siguma)
Az=AzEl[0]
El=AzEl[1]
TisinSisaAz=calcu_TisinSisaAz(R,k_ido,Az)
TisinSisaEl=calcu_TisinSisaEl(R,El)
print Az-TisinSisaAz,El-TisinSisaEl
