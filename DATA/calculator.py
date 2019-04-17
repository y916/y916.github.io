#!/usr/bin/python
# -*- coding: utf-8 –*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import os

def calculator(ONM,MDS):
    mdb.Job(name=ONM, model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=2, 
        numDomains=2, numGPUs=0)
    mdb.jobs[ONM].submit()
    mdb.jobs[ONM].waitForCompletion()
    mdb.save()

    odb = session.openOdb(name=ONM+'.odb')
    myViewport = session.viewports['Viewport: 1']
    myViewport.setValues(displayedObject=odb)
    myViewport.odbDisplay.display.setValues(plotState=(UNDEFORMED, ))
    session.printToFile(fileName='image', format=PNG,canvasObjects=(myViewport,))

    session.XYDataFromHistory(name='Frequency', odb=odb, 
        outputVariableName='Eigenfrequency: EIGFREQ for Whole Model', steps=('Step-1', ), )
    x0 = session.xyDataObjects['Frequency']
    session.writeXYReport(fileName='Frequency.rpt', xyData=(x0, ))
    del session.xyDataObjects['Frequency']
    session.odbData[ONM+'.odb'].setValues(activeFrames=(('Step-2', ('0:-1', )), ))
    session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('A', 
        NODAL, ((INVARIANT, 'Magnitude'), )), ), numericForm=COMPLEX_MAGNITUDE, nodeSets=MDS)

    aac = session.xyDataObjects.items()
    for i in range(len(MDS)):
        x0 = session.xyDataObjects[aac[i][0]]
        session.writeXYReport(fileName=MDS[i] + '.rpt', appendMode=OFF, xyData=(x0, ))
        del session.xyDataObjects[aac[i][0]]

    session.odbs[ONM+'.odb'].close()
