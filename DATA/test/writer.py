#!/usr/bin/python
# -*- coding: mbcs -*-

from abaqus import *

aac = session.xyDataObjects.items()
for i in range(len(aac)):
    x0 = session.xyDataObjects[aac[i][0]]
    session.writeXYReport(fileName=MDS[i] + '.rpt', appendMode=OFF, xyData=(x0, ))
    del session.xyDataObjects[aac[i][0]]

session.odbs[ONM].close()
