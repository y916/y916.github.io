# -*- coding: utf-8 –*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *
from material import *
import os

MXA = ['T35H35H35H35H66']
MDS = ['M88','T88','M68','T68','M86','T86','M48','T48','M84','T84','M28','T28','M82','T82','M80','T80']
WP0 = 'H'
Density0 = 1.5233e-09
Material0 = (50500.0, 50500.0, 10100.0, 0.3, 0.32, 0.32, 4000.0, 4800.0, 4800.0)

TM = {'6':30.6}
LH = {'0':0.0,'1':1.0}
FTX = {'1':pi/12.0,'2':pi/6.0,'3':pi/4.0,'4':pi/3.0}
HG = {'0':0,'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'A':10,
        'B':11,'C':12,'D':13,'E':14,'F':15,'G':16}
XG = {'0':0.0,'1':0.25,'2':0.5,'3':0.75,'4':1.0,'5':1.25,'6':1.5,'7':1.75,'8':2.0,'9':2.25,'A':2.5,
        'B':2.75,'C':3.0,'D':3.25,'E':3.5,'F':3.75,'G':4.0}
if not os.path.exists(WP0+':'):
    if os.path.exists('G:'):
        WP0 = 'G'
    else:
        WP0 = 'D'

def dzero(num0):
    A = num0[:]
    while A[0] is '0':
        A = A[1:]
        if len(A) == 0:
            A = '0'
            break
    return A

def dmask(mask0):
    MASK0 = mask0[:]
    h = len(MASK0)//8
    MASK1 = ''
    for i in range(h):
        MASK1 = MASK1 + '#'+dzero(MASK0[-8:])+' '
        MASK0 = MASK0[:-8]
    if len(MASK0)>0:
        MASK1 = '['+MASK1 + '#'+MASK0+' ]'
    else:
        MASK1 = '['+MASK1 +']'
    return MASK1

myViewport = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200.0, height=100.0)
myViewport.makeCurrent()
myViewport.maximize()
executeOnCaeStartup()
myViewport.partDisplay.geometryOptions.setValues(referenceRepresentation=ON)

for MDL0 in MXA:
    T0 = TM[MDL0[-1]]
    NUM0 = HG[MDL0[-2]]
    LEN0 = T0*NUM0
    W0 = 3.0
    W1 = round(LEN0/16.0/W0)
    W = LEN0/16.0/W1

    MDL = MDL0[:-2]
    WP = WP0 + ':\\ABAQUS\\' + MDL0
    if not os.path.exists(WP):
        os.makedirs(WP)
    os.chdir(WP)
    TF = open(MDL0+'.txt','wt')
    TF.write(MDL0+'\n\n\n'+'\nWork Path: '+WP+'\n')
    TF.write('Density = '+str(Density0)+'\n'+'Material = '+str(Material0)+'\n')
    TF.write('Cycle = '+str(T0)+'\n'+'Number = '+str(NUM0)+'\n'+'Length = '+str(LEN0)+'\n')
    TF.write('Grid Length = '+str(W)+'\n'+'Number = '+str(W1)+'\n'+'Panel Number = '+str(int(W1*16)**2)+'\n\n')

    XZ = []; XZK = {}
    for i in range(len(MDL[1:])/3):
        XZ.append(MDL[3*i+1:3*i+3])
        XZK[MDL[3*i+1:3*i+3]] = 'T' + MDL[3*i+1:3*i+3]
    MX = sorted(set(XZ))

    Mdb()
    myViewport.setValues(displayedObject=None)

    s = mdb.models['Model-1'].ConstrainedSketch(name='MB1', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0, -1.0), point2=(T0*NUM0, 0.0))
    p = mdb.models['Model-1'].Part(name='MB1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['MB1']
    p.BaseSolidExtrude(sketch=s, depth=LEN0)
    s.unsetPrimaryObject()

    for i in MX:
        ttx = XG[i[1]]/2.0*tan(FTX[i[0]]/2.0)
        tty = XG[i[1]]
        s = mdb.models['Model-1'].ConstrainedSketch(name=XZK[i], sheetSize=200.0)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        s.Line(point1=(0.0, 0.0), point2=((T0-20)/4.0+ttx, 0.0)) 
        s.HorizontalConstraint(entity=g[2], addUndoState=False)
        s.Line(point1=((T0-20)/4.0+ttx, 0.0), point2=((T0+20)/4.0+ttx, 10.0*tan(FTX[i[0]])))
        s.Line(point1=((T0+20)/4.0+ttx, 10.0*tan(FTX[i[0]])), point2=(T0/2.0, 10.0*tan(FTX[i[0]])))
        s.HorizontalConstraint(entity=g[4], addUndoState=False)
        s.Line(point1=(T0/2.0, 10.0*tan(FTX[i[0]])), point2=(T0/2.0, 10.0*tan(FTX[i[0]])+tty))
        s.VerticalConstraint(entity=g[5], addUndoState=False)
        s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
        s.Line(point1=(T0/2.0, 10.0*tan(FTX[i[0]])+tty), point2=((T0+20)/4.0-ttx, 10.0*tan(FTX[i[0]])+tty))
        s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
        s.HorizontalConstraint(entity=g[6], addUndoState=False)
        s.Line(point1=((T0+20)/4.0-ttx, 10.0*tan(FTX[i[0]])+tty), point2=((T0-20)/4.0-ttx, 0.0+tty))
        s.Line(point1=((T0-20)/4.0-ttx, 0.0+tty), point2=(0.0, 0.0+tty))
        s.HorizontalConstraint(entity=g[8], addUndoState=False)
        s.copyMirror(mirrorLine=g[5], objectList=(g[2], g[3], g[4], g[6], g[7], g[8]))
        s.delete(objectList=(g[5], ))
        if NUM0>1:
            s.linearPattern(geomList=(g[2], g[3], g[4], g[6], g[7], g[8], g[9], g[10], 
            g[11], g[12], g[13], g[14]), vertexList=(), number1=NUM0, spacing1=T0, angle1=0.0, 
            number2=1, spacing2=20.0, angle2=90.0)
        s.Line(point1=(0.0, 0.0), point2=(0.0, 0.0+tty))
        s.VerticalConstraint(entity=g[12*NUM0+3], addUndoState=False)
        s.PerpendicularConstraint(entity1=g[2], entity2=g[12*NUM0+3], addUndoState=False)
        s.PerpendicularConstraint(entity1=g[8], entity2=g[12*NUM0+3], addUndoState=False)
        s.Line(point1=(LEN0, 0.0), point2=(T0*NUM0, 0.0+tty))
        s.VerticalConstraint(entity=g[12*NUM0+4], addUndoState=False)
        s.PerpendicularConstraint(entity1=g[12*NUM0+2], entity2=g[12*NUM0+4], addUndoState=False)
        s.PerpendicularConstraint(entity1=g[12*NUM0-3], entity2=g[12*NUM0+4], addUndoState=False)
        p = mdb.models['Model-1'].Part(name=XZK[i], dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts[XZK[i]]
        p.BaseSolidExtrude(sketch=s, depth=LEN0)
        s.unsetPrimaryObject()

    myViewport.partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON)
    myViewport.partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)

    mdb.models['Model-1'].Material(name='fc3')
    mdb.models['Model-1'].materials['fc3'].Density(table=((Density0, ), ))
    mdb.models['Model-1'].materials['fc3'].Elastic(
        type=ENGINEERING_CONSTANTS, table=(Material0, ))

    mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material='fc3', thickness=None)

    p = mdb.models['Model-1'].parts['MB1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(cells=cells, name='Set-1')
    p = mdb.models['Model-1'].parts['MB1']
    p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

    for i in MX:
        p = mdb.models['Model-1'].parts[XZK[i]]
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = p.Set(cells=cells, name='Set-1')
        p = mdb.models['Model-1'].parts[XZK[i]]
        p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

    p = mdb.models['Model-1'].parts['MB1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(cells=cells)
    p = mdb.models['Model-1'].parts['MB1']
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#2 ]', ), )
    normalAxisRegion = p.Surface(side1Faces=side1Faces, name='Surf-1')
    p = mdb.models['Model-1'].parts['MB1']
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#20 ]', ), )
    primaryAxisRegion = p.Set(edges=edges, name='Set-2')
    mdb.models['Model-1'].parts['MB1'].MaterialOrientation(region=region, 
        orientationType=DISCRETE, axis=AXIS_1, normalAxisDefinition=SURFACE, 
        normalAxisRegion=normalAxisRegion, flipNormalDirection=False, 
        normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
        primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_1, 
        flipPrimaryDirection=False, additionalRotationType=ROTATION_NONE, 
        angle=0.0, additionalRotationField='', stackDirection=STACK_3)

    mask0 = '1'+'f'*NUM0
    mask1 = '[#0:'+str(NUM0)+' #10 ]'
    for i in MX:
        p = mdb.models['Model-1'].parts[XZK[i]]
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = regionToolset.Region(cells=cells)
        p = mdb.models['Model-1'].parts[XZK[i]]
        s = p.faces
        side1Faces = s.getSequenceFromMask(mask=(dmask(mask0), ), )
        normalAxisRegion = p.Surface(side1Faces=side1Faces, name='Surf-1')
        p = mdb.models['Model-1'].parts[XZK[i]]
        e = p.edges
        edges = e.getSequenceFromMask(mask=(mask1, ), )
        primaryAxisRegion = p.Set(edges=edges, name='Set-2')
        mdb.models['Model-1'].parts[XZK[i]].MaterialOrientation(region=region, 
        orientationType=DISCRETE, axis=AXIS_1, normalAxisDefinition=SURFACE, 
        normalAxisRegion=normalAxisRegion, flipNormalDirection=False, 
        normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
        primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_1, 
        flipPrimaryDirection=False, additionalRotationType=ROTATION_NONE, 
        angle=0.0, additionalRotationField='', stackDirection=STACK_3)

    h = 0
    DT = [0]
    for i in XZ:
        h = h + 10*tan(FTX[i[0]]) + XG[i[1]]
        DT.append(h)
    TF.write('DT = '+str(DT)+'\n')

    myViewport.assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    for i in range(len(XZ)+1):
        a = mdb.models['Model-1'].rootAssembly
        p = mdb.models['Model-1'].parts['MB1']
        a.Instance(name='MB1-'+str(i), part=p, dependent=ON)
    a.DatumCsysByDefault(CARTESIAN)

    for i in range(len(XZ)):
        a = mdb.models['Model-1'].rootAssembly
        p = mdb.models['Model-1'].parts[XZK[XZ[i]]]
        a.Instance(name=XZK[XZ[i]] + '-' + str(i+1), part=p, dependent=ON)

    for i in range(len(XZ)):
        a = mdb.models['Model-1'].rootAssembly
        a.translate(instanceList=('MB1-' + str(i+1), ), vector=(0.0, DT[i+1]+i+1, 0.0))

    for i in range(len(XZ)):
        a = mdb.models['Model-1'].rootAssembly
        if i > 0:
            a.translate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), vector=(0.0, DT[i]+i, 0.0))
        if MDL[3*i+3] is 'J':
            a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), 
                axisPoint=(LEN0/2.0, (DT[i]+DT[i+1])/2.0+i, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)
        elif LEN0 == T0*NUM0:
            if MDL[3*i+3] is 'I':
                a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), axisPoint=(LEN0/2.0, 0.0, LEN0/2.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
            elif MDL[3*i+3] is 'K':
                a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), axisPoint=(LEN0/2.0, 0.0, LEN0/2.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
                a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), 
                    axisPoint=(LEN0/2.0, (DT[i]+DT[i+1])/2.0+i, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)

    a = mdb.models['Model-1'].rootAssembly
    myViewport.setValues(displayedObject=a)
    myViewport.view.setValues(session.views['Iso'])
    session.printToFile(fileName=MDL0, format=PNG,canvasObjects=(myViewport,))
    myViewport.setValues(displayedObject=None)
    myViewport.assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON)

    mask2 = '4'*(NUM0+1)+'0'*NUM0
    mask3 = '4'*NUM0
    for i in range(len(XZ)):
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances['MB1-' + str(i)].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
        region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-' + str(2*i+1))
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances[XZK[XZ[i]] + '-' + str(i+1)].faces
        if MDL[3*i+3] in ['H','I']:
            side1Faces1 = s1.getSequenceFromMask(mask=(dmask(mask2), ), )
        else:
            side1Faces1 = s1.getSequenceFromMask(mask=(dmask(mask3), ), )
        region2=a.Surface(side1Faces=side1Faces1, name='s_Surf-' + str(2*i+1))
        mdb.models['Model-1'].Tie(name='Constraint-' + str(2*i+1), master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    for i in range(len(XZ)):
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances['MB1-' + str(i+1)].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
        region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-' + str(2*i+2))
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances[XZK[XZ[i]] + '-' + str(i+1)].faces
        if MDL[3*i+3] in ['H','I']:
            side1Faces1 = s1.getSequenceFromMask(mask=(dmask(mask3), ), )
        else:
            side1Faces1 = s1.getSequenceFromMask(mask=(dmask(mask2), ), )
        region2=a.Surface(side1Faces=side1Faces1, name='s_Surf-' + str(2*i+2))
        mdb.models['Model-1'].Tie(name='Constraint-' + str(2*i+2), master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    myViewport.assemblyDisplay.setValues(mesh=ON, 
        interactions=OFF, constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
    myViewport.assemblyDisplay.meshOptions.setValues(
        meshTechnique=ON)

    elemType1 = mesh.ElemType(elemCode=C3D8I, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    #1
    p = mdb.models['Model-1'].parts['MB1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    pickedRegions =(cells, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
    p = mdb.models['Model-1'].parts['MB1']
    c = p.cells
    pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=HEX)
    p = mdb.models['Model-1'].parts['MB1']
    p.seedPart(size=W, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['MB1']
    p.generateMesh()

    for i in MX:
        p = mdb.models['Model-1'].parts[XZK[i]]
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        p = mdb.models['Model-1'].parts[XZK[i]]
        c = p.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        p.setMeshControls(regions=pickedRegions, elemShape=HEX_DOMINATED)
        p = mdb.models['Model-1'].parts[XZK[i]]
        p.seedPart(size=W, deviationFactor=0.1, minSizeFactor=0.1)
        p = mdb.models['Model-1'].parts[XZK[i]]
        p.generateMesh()

    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()

    for JHD in MDS:
        j = JHD[1]; k = JHD[2]
        JDMOD = 10**int((HG[j]*(W1*16.0+1.0)*W1+HG[k]*W1)%16//2)
        if (HG[j]*(W1*16.0+1.0)*W1+HG[k]*W1)%16%2 == 1:
            JDMOD = JDMOD*4
        a = mdb.models['Model-1'].rootAssembly
        if JHD[0] is 'M':
            mask4 = '[#0:' + str(int(HG[j]*(W1*16.0+1.0)*W1+HG[k]*W1)//16) + ' #' + str(JDMOD) + ' ]'
            n1 = a.instances['MB1-' + str(len(XZ))].nodes
        else:
            mask4 = '[#0:' + str(int(HG[j]*(W1*16.0+1.0)*W1+HG[k]*W1)//16) + ' #' + str(2*JDMOD) + ' ]'
            n1 = a.instances['MB1-0'].nodes
        nodes1 = n1.getSequenceFromMask(mask=(mask4, ), )
        a.Set(nodes=nodes1, name=JHD)

    myViewport.assemblyDisplay.setValues(mesh=OFF, adaptiveMeshConstraints=ON)
    myViewport.assemblyDisplay.meshOptions.setValues(meshTechnique=OFF)
    mdb.saveAs(pathName=WP + '\\' + MDL0)
    TF.write('Output point = '+str(MDS)+'\n')
    TF.close()