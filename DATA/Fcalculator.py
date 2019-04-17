#!/usr/bin/python
# -*- coding: utf-8 -*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import os
import Tkinter, tkFileDialog

MDY = 'LR04'
MDZ = 'P0M88'
MDS = ['M88','T88','M68','T68','M86','T86','M48','T48','M84','T84','M28','T28','M82','T82','M80','T80']

Frequency0 = 1
MaxFrequency0 = 5000.0
FrequencyPoint0 = 20
ncs=1
LGA=1
TM = {'6':30.6}
HG = {'0':0,'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'A':10,
        'B':11,'C':12,'D':13,'E':14,'F':15,'G':16}

root = Tkinter.Tk()
root.withdraw()
filepath0 = tkFileDialog.askopenfilename()
filepath = filepath0.encode('gbk')
for i in range(len(filepath)):
    if filepath[-(i+1)] is '/':
        break
WP0=filepath[0:-i-1]
WP=WP0+'\\'+MDY+MDZ
if not os.path.exists(WP):
    os.makedirs(WP)
MDL0 = filepath[-i:-4]
LXZ = (len(MDL0)-5)/3
T0 = TM[MDL0[-3]]
NUM0 = HG[MDL0[-4]]
LEN0 = T0*NUM0
W0 = 3.0
W1 = round(LEN0/16.0/W0)
W = LEN0/16.0/W1
os.chdir(WP)

session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200.0, 
    height=100.0)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)

openMdb(pathName = filepath)
session.viewports['Viewport: 1'].setValues(displayedObject=None)

#1
mdb.models['Model-1'].FrequencyStep(name='Step-1', previous='Initial', maxEigen=MaxFrequency0, numEigen=Frequency0)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
#2
mdb.models['Model-1'].SteadyStateModalStep(name='Step-2', previous='Step-1', 
    frequencyRange=((1.0, MaxFrequency0, FrequencyPoint0, 3.0), ), directDamping=((1, Frequency0, 0.01), ))
TF = open(MDL0+'-MoTai.txt','wt')
TF.write('numEigen = '+str(Frequency0)+'\nmaxEigen = '+str(MaxFrequency0)+'\n')
TF.write('frequencyRange = '+str((1.0, MaxFrequency0, FrequencyPoint0, 3.0))+'\n')
TF.write('directDamping = '+str((1, Frequency0, 0.01))+'\n')

session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-2')
mdb.models['Model-1'].fieldOutputRequests['F-Output-2'].setValues(variables=('S', 'U', 'V', 'A'))

session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)

a = mdb.models['Model-1'].rootAssembly
FM = float(2**(HG[MDZ[1]]))
v1 = a.instances['MB1-0'].vertices
v11 = a.instances['MB1-' + str(LXZ)].vertices
a = mdb.models['Model-1'].rootAssembly
if MDZ[0] is 'P':
    if MDZ[2] is 'M':
        FM = -FM
    mdb.models['Model-1'].ConcentratedForce(name='Load-F', createStepName='Step-2', 
    region=a.sets[MDS[0]], cf2=FM+0j, distributionType=UNIFORM, field='', localCsys=None)
else:
    if MDZ[2] is 'T':
        s1 = a.instances['MB1-0'].faces
    else:
        s1 = a.instances['MB1-' + str(LXZ)].faces
        FM = -FM
    side1Faces1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
    region = a.Surface(side1Faces=side1Faces1, name='Surf-FM')
    mdb.models['Model-1'].SurfaceTraction(name='Load-M', createStepName='Step-2', 
    region=region, magnitude=FM+0j, directionVector=(v1[6], v11[6]), 
    distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL)

YS={}; ys0 = []; h = 0
for i in range(LXZ+1):
    YS['MB1-' + str(i)] = []
for i in range(len(MDY)):
    if i is len(MDY)-1 or (not MDY[i] in ['R','L','N','S','M','T']
                                and MDY[i+1] in ['R','L','N','S','M','T']):
        ys0.append(MDY[h:i+1])
        h = i + 1
for i in ys0:
    if i[-2] in ['R','L','N','S','M','T']:
        for j in i[:-1]:
            YS['MB1-' + i[-1]].append(j)
    else:
        for j in range(HG[i[-2]],HG[i[-1]]+1):
            for k in i[:-2]:
                YS['MB1-' + str(j)].append(k)
for key,value in YS.items():
    if value == []:
        del YS[key]
    else:
        YS[key] = sorted(set(value))
del ys0,h
TF.write('Restriction = '+str(YS)+'\n')
TF.close()

session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
yvk = {'R':4,'L':1,'N':32,'S':16,'M':2,'T':8}
for key,value in YS.items():
    yvx=0
    for yv in value:
        yvx = yvx + yvk[yv]
    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances[key].faces
    faces1 = f1.getSequenceFromMask(mask=('[#' + hex(yvx)[2:] + ' ]', ), )
    region = a.Set(faces=faces1, name='YS-' + key)
    mdb.models['Model-1'].DisplacementBC(name='BC-' + key, createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.Job(name=MDL0, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=ncs, 
    numDomains=ncs, numGPUs=0)

mdb.jobs[MDL0].submit()
mdb.jobs[MDL0].waitForCompletion()
mdb.saveAs(pathName=WP + '\\' + MDL0)
ONM = WP + '\\' + MDL0 + '.odb'
odb = session.openOdb(name=ONM)
myViewport = session.viewports['Viewport: 1']
myViewport.setValues(displayedObject=odb)
myViewport.odbDisplay.display.setValues(plotState=(UNDEFORMED, ))
session.printToFile(fileName=MDL0, format=PNG,canvasObjects=(myViewport,))

session.XYDataFromHistory(name='Frequency', odb=odb, 
    outputVariableName='Eigenfrequency: EIGFREQ for Whole Model', steps=('Step-1', ), )
x0 = session.xyDataObjects['Frequency']
session.writeXYReport(fileName=MDL0 + '_Frequency.rpt', xyData=(x0, ))
del session.xyDataObjects['Frequency']
session.odbData[ONM].setValues(activeFrames=(('Step-2', ('0:-1', )), ))
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('A', 
    NODAL, ((INVARIANT, 'Magnitude'), )), ), numericForm=COMPLEX_MAGNITUDE, nodeSets=MDS)

for i in range(len(WP0)):
    if WP0[-(i+1)] is '/':
        break
WT0=WP0[0:-i]
WT=WT0+'writer.py'
execfile(WT)
