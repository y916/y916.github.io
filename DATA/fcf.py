# -*- coding: utf-8 –*-

import os
import shutil
import sys

WP0 = r'H:\ABAQUS'
MXA = ['S35H35H35H35H']
MDS = []
GG = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G']
HG = {'0':0,'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'A':10,
        'B':11,'C':12,'D':13,'E':14,'F':15,'G':16}
for i in GG:
    for j in GG:
        MDS.append('T'+i+j)
for i in GG:
    for j in GG:
        MDS.append('M'+i+j)
Density0 = 1.5233e-09
Material0 = (50500.0, 50500.0, 10100.0, 0.3, 0.32, 0.32, 4000.0, 4800.0, 4800.0)
NUM0 = 6
T0 = 30.6
LEN0 = NUM0*T0
W0 = 3.0
W1 = round(T0*NUM0/16.0/W0)
TFTX = {}
SFTX = {}
FST0 = [0.1,0.3,0.7,1.1,1.2,1.3,1.4,1.5]
STW = 1.0
STH = 0.0
sys.path.append(WP0+r'\DATA')

import calculator

from Builder import Builder
for MDL0 in MXA:

    if MDL0[0] is 'T':
        Builder(MDL0,WP0,NUM0,T0,W0,TFTX,SFTX,FST0,STW,STH).tbuilder()
    elif MDL0[0] is 'S':
        Builder(MDL0,WP0,NUM0,T0,W0,TFTX,SFTX,FST0,STW,STH).sbuilder()

    WP2 = WP0 + '\\'+MDL0[0] + 'model\\'+MDL0

    '''
    import jzt
    ONM = 'JZT'
    WP = WP2+'\\'+ONM
    if not os.path.exists(WP):
        os.makedirs(WP)
    os.chdir(WP)
    shutil.copy(WP2+'\\'+MDL0+'.cae',ONM+'.cae')
    shutil.copy(WP2+'\\'+MDL0+'.jnl',ONM+'.jnl')
    jzt.jzt(ONM,LEN0,Density0,Material0)
    calculator.calculator(ONM,MDS)
    '''

    import jz
    ONM = 'JZ'
    WP = WP2+'\\'+ONM
    if not os.path.exists(WP):
        os.makedirs(WP)
    os.chdir(WP)
    shutil.copy(WP2+'\\'+MDL0+'.cae',ONM+'.cae')
    shutil.copy(WP2+'\\'+MDL0+'.jnl',ONM+'.jnl')
    jz.jz(ONM,LEN0,Density0,Material0)
    calculator.calculator(ONM,MDS)
