# -*- coding: utf-8 –*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *
from material import *
import os

def jzt(ONM,LEN0,Density0,Material0):

    openMdb(pathName=ONM+'.cae')
    a = mdb.models['Model-1'].rootAssembly
    a.ReferencePoint(point=(LEN0/2.0, -2.0, LEN0/2.0))

    s = mdb.models['Model-1'].ConstrainedSketch(name='JZT', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, -3.0), point2=(LEN0, -1.0))
    p = mdb.models['Model-1'].Part(name='JZT', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['JZT']
    p.BaseSolidExtrude(sketch=s, depth=LEN0)
    s.unsetPrimaryObject()

    mdb.models['Model-1'].materials['FC'].Density(table=((Density0, ), ))
    mdb.models['Model-1'].materials['FC'].Elastic(type=ENGINEERING_CONSTANTS, table=(Material0, ))
    mdb.models['Model-1'].Material(name='Material-2')
    mdb.models['Model-1'].materials['Material-2'].Density(table=((7.8e-09, ), ))
    mdb.models['Model-1'].materials['Material-2'].Elastic(table=((200000.0, 0.3), 
        ))
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section-2', 
        material='Material-2', thickness=None)
    p = mdb.models['Model-1'].parts['JZT']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(cells=cells, name='Set-1')
    p = mdb.models['Model-1'].parts['JZT']
    p.SectionAssignment(region=region, sectionName='Section-2', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        adaptiveMeshConstraints=OFF)
    a = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts['JZT']
    a.Instance(name='JZT-1', part=p, dependent=ON)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].FrequencyStep(name='Step-1', previous='Initial', maxEigen=5000.0)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    mdb.models['Model-1'].SteadyStateModalStep(name='Step-2', previous='Step-1', 
        frequencyRange=((0.0, 5000.0, 20, 3.0), ), directDamping=((1, 100, 0.005), ))
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-2')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON, adaptiveMeshConstraints=OFF)
    
    '''
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=FRICTIONLESS)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=OFF, 
        constraintEnforcementMethod=DEFAULT)
    mdb.models['Model-1'].ContactStd(name='Int-1', createStepName='Initial')
    mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, 'IntProp-1'), ))
    '''
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['JZT-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
    region1=regionToolset.Region(side1Faces=side1Faces1)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['MB1-0'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
    region2=regionToolset.Region(side1Faces=side1Faces1)
    mdb.models['Model-1'].Tie(name='Constraint-tie', master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['JZT-1'].cells
    cells1 = c1.getSequenceFromMask(mask=('[#1 ]', ), )
    region2=a.Set(cells=cells1, name='JZT')
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[326], )
    region1=regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].RigidBody(name='Constraint-JZT', refPointRegion=region1, 
        bodyRegion=region2)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, interactions=OFF, constraints=OFF, 
        engineeringFeatures=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-2')
    a = mdb.models['Model-1'].rootAssembly

    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].ConcentratedForce(name='Load-1', createStepName='Step-2', 
        region=region, cf2=0+1j, distributionType=UNIFORM, field='', localCsys=None)

    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
        region=region, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-2')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
        bcs=OFF, predefinedFields=OFF, connectors=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(meshTechnique=ON)

    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
        engineeringFeatures=OFF, mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    p = mdb.models['Model-1'].parts['JZT']
    p.seedPart(size=LEN0/32, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['JZT']
    p.generateMesh()
    a = mdb.models['Model-1'].rootAssembly
    a1 = mdb.models['Model-1'].rootAssembly
    a1.regenerate()
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)

    mdb.save()