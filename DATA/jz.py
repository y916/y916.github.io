# -*- coding: utf-8 –*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *
from material import *
import os

def jz(ONM,LEN0,Density0,Material0):

    openMdb(pathName=ONM+'.cae')
    mdb.models['Model-1'].materials['FC'].Density(table=((Density0, ), ))
    mdb.models['Model-1'].materials['FC'].Elastic(type=ENGINEERING_CONSTANTS, table=(Material0, ))

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].FrequencyStep(name='Step-1', previous='Initial', maxEigen=5000.0)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    mdb.models['Model-1'].SteadyStateModalStep(name='Step-2', previous='Step-1', 
        frequencyRange=((0.0, 5000.0, 20, 3.0), ), directDamping=((1, 100, 0.005), ))
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-2')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON, adaptiveMeshConstraints=OFF)
    
    a = mdb.models['Model-1'].rootAssembly
    a.ReferencePoint(point=(LEN0/2.0, -1.0, LEN0/2.0))

    a = mdb.models['Model-1'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[326], )
    region1=regionToolset.Region(referencePoints=refPoints1)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['MB1-0'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
    region2=regionToolset.Region(side1Faces=side1Faces1)
    mdb.models['Model-1'].Coupling(name='Constraint-O', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, 
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, interactions=OFF, constraints=OFF, engineeringFeatures=OFF)
    
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-2')
    a = mdb.models['Model-1'].rootAssembly
    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].ConcentratedForce(name='Load-1', createStepName='Step-2', 
        region=region, cf2=0+1j, distributionType=UNIFORM, field='', localCsys=None)

    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
        region=region, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    mdb.save()
