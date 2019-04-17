#!/usr/bin/python
# -*- coding: utf-8 –*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import os
import Tkinter, tkFileDialog

for i in MDS:
    j = i[1]; k = i[2]
    MBJD = int(HG[j]*(W1*16.0+1.0)*W1+HG[k]*W1)*2 + 1
    if i[0] is 'M':
        x0 = session.xyDataObjects['A:Magnitude   复数: 大小 PI: MB1-'+str(LXZ)+' N: ' + str(MBJD)]
        session.writeXYReport(fileName=MDL0 + '_' + i + '.rpt', appendMode=OFF, xyData=(x0, ))
        del session.xyDataObjects['A:Magnitude   复数: 大小 PI: MB1-'+str(LXZ)+' N: ' + str(MBJD)]
    else:
        x0 = session.xyDataObjects['A:Magnitude   复数: 大小 PI: MB1-0 N: ' + str(MBJD + 1)]
        session.writeXYReport(fileName=MDL0 + '_' + i + '.rpt', appendMode=OFF, xyData=(x0, ))
        del session.xyDataObjects['A:Magnitude   复数: 大小 PI: MB1-0 N: ' + str(MBJD + 1)]

session.odbs[ONM].close()
print('SUCCESSFUL! The data of '+MDL0+' has been saved in '+WP)
