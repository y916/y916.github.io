# -*- coding: utf-8 –*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import *
from material import *
import os

GG = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G']
HG = {'0':0,'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'A':10,
        'B':11,'C':12,'D':13,'E':14,'F':15,'G':16}
XG = {'0':0.0,'1':0.25,'2':0.5,'3':0.75,'4':1.0,'5':1.25,'6':1.5,'7':1.75,'8':2.0,'9':2.25,'A':2.5,
        'B':2.75,'C':3.0,'D':3.25,'E':3.5,'F':3.75,'G':4.0}
MDS0 = []
for i in GG:
    for j in GG:
        MDS0.append('T'+i+j)
        MDS0.append('M'+i+j)

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

class Builder(object):
    def __init__(self,MDL0,WP0,NUM0,T0,W0,TFTX,SFTX,FST0,STW,STH):
        self.MDL0 = MDL0
        self.WP0 = WP0
        self.NUM0 = NUM0
        self.T0 = T0
        self.W0 = W0
        self.TFTX = TFTX
        self.SFTX = SFTX
        self.FST0 = FST0
        self.STW = STW
        self.STH = STH
        self.LEN0 = NUM0*T0
        self.W1 = round(self.LEN0/16.0/self.W0)
        self.W = self.LEN0/16.0/self.W1
        self.MBH = XG['4']

    def tbuilder(self):

        MDL = self.MDL0[:]
        FTX0 = {'1':pi/12.0,'2':pi/6.0,'3':pi/4.0,'4':pi/3.0}
        FTX = FTX0
        for key, value in self.TFTX.items():
            FTX[key] = value

        myViewport = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200.0, height=100.0)
        myViewport.makeCurrent()
        myViewport.maximize()
        executeOnCaeStartup()
        myViewport.partDisplay.geometryOptions.setValues(referenceRepresentation=ON)

        WP = self.WP0 + '\\Tmodel\\' + MDL
        if not os.path.exists(WP):
            os.makedirs(WP)
        os.chdir(WP)
        TF = open(MDL+'.txt','wt')
        TF.write(MDL+'\n\n\n'+'\nWork Path: '+WP+'\n')
        TF.write('Cycle = '+str(self.T0)+'\n'+'Number = '+str(self.NUM0)+'\n'+'Length = '+str(self.LEN0)+'\n')
        TF.write('Grid Length = '+str(self.W)+'\n'+'Number = '+str(self.W1)+'\n'+'Panel Number = '+str(int(self.W1*16)**2)+'\n\n')

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
        s.rectangle(point1=(0, -self.MBH), point2=(self.T0*self.NUM0, 0.0))
        p = mdb.models['Model-1'].Part(name='MB1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts['MB1']
        p.BaseSolidExtrude(sketch=s, depth=self.LEN0)
        s.unsetPrimaryObject()

        for i in MX:
            ttx = XG[i[1]]/2.0*tan(FTX[i[0]]/2.0)
            tty = XG[i[1]]
            s = mdb.models['Model-1'].ConstrainedSketch(name=XZK[i], sheetSize=200.0)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=STANDALONE)
            s.Line(point1=(0.0, 0.0), point2=((self.T0-20)/4.0+ttx, 0.0)) 
            s.HorizontalConstraint(entity=g[2], addUndoState=False)
            s.Line(point1=((self.T0-20)/4.0+ttx, 0.0), point2=((self.T0+20)/4.0+ttx, 10.0*tan(FTX[i[0]])))
            s.Line(point1=((self.T0+20)/4.0+ttx, 10.0*tan(FTX[i[0]])), point2=(self.T0/2.0, 10.0*tan(FTX[i[0]])))
            s.HorizontalConstraint(entity=g[4], addUndoState=False)
            s.Line(point1=(self.T0/2.0, 10.0*tan(FTX[i[0]])), point2=(self.T0/2.0, 10.0*tan(FTX[i[0]])+tty))
            s.VerticalConstraint(entity=g[5], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
            s.Line(point1=(self.T0/2.0, 10.0*tan(FTX[i[0]])+tty), point2=((self.T0+20)/4.0-ttx, 10.0*tan(FTX[i[0]])+tty))
            s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
            s.HorizontalConstraint(entity=g[6], addUndoState=False)
            s.Line(point1=((self.T0+20)/4.0-ttx, 10.0*tan(FTX[i[0]])+tty), point2=((self.T0-20)/4.0-ttx, 0.0+tty))
            s.Line(point1=((self.T0-20)/4.0-ttx, 0.0+tty), point2=(0.0, 0.0+tty))
            s.HorizontalConstraint(entity=g[8], addUndoState=False)
            s.copyMirror(mirrorLine=g[5], objectList=(g[2], g[3], g[4], g[6], g[7], g[8]))
            s.delete(objectList=(g[5], ))
            if self.NUM0>1:
                s.linearPattern(geomList=(g[2], g[3], g[4], g[6], g[7], g[8], g[9], g[10], 
                g[11], g[12], g[13], g[14]), vertexList=(), number1=self.NUM0, spacing1=self.T0, angle1=0.0, 
                number2=1, spacing2=20.0, angle2=90.0)
            s.Line(point1=(0.0, 0.0), point2=(0.0, 0.0+tty))
            s.VerticalConstraint(entity=g[12*self.NUM0+3], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[2], entity2=g[12*self.NUM0+3], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[8], entity2=g[12*self.NUM0+3], addUndoState=False)
            s.Line(point1=(self.T0*self.NUM0, 0.0), point2=(self.T0*self.NUM0, 0.0+tty))
            s.VerticalConstraint(entity=g[12*self.NUM0+4], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[12*self.NUM0+2], entity2=g[12*self.NUM0+4], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[12*self.NUM0-3], entity2=g[12*self.NUM0+4], addUndoState=False)
            p = mdb.models['Model-1'].Part(name=XZK[i], dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p = mdb.models['Model-1'].parts[XZK[i]]
            p.BaseSolidExtrude(sketch=s, depth=self.LEN0)
            s.unsetPrimaryObject()

        myViewport.partDisplay.setValues(sectionAssignments=ON, 
            engineeringFeatures=ON)
        myViewport.partDisplay.geometryOptions.setValues(
            referenceRepresentation=OFF)

        mdb.models['Model-1'].Material(name='FC')

        mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material='FC', thickness=None)

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

        mask0 = '1'+'f'*self.NUM0
        mask1 = '[#0:'+str(self.NUM0)+' #10 ]'
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
            a.translate(instanceList=('MB1-' + str(i+1), ), vector=(0.0, DT[i+1]+(i+1)*self.MBH, 0.0))

        for i in range(len(XZ)):
            a = mdb.models['Model-1'].rootAssembly
            if i > 0:
                a.translate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), vector=(0.0, DT[i]+i*self.MBH, 0.0))
            if MDL[3*i+3] is 'J':
                a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), 
                    axisPoint=(self.LEN0/2.0, (DT[i]+DT[i+1])/2.0+i*self.MBH, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)
            elif self.LEN0 == self.T0*self.NUM0:
                if MDL[3*i+3] is 'I':
                    a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), axisPoint=(self.LEN0/2.0, 0.0, self.LEN0/2.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
                elif MDL[3*i+3] is 'K':
                    a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), axisPoint=(self.LEN0/2.0, 0.0, self.LEN0/2.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
                    a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), 
                        axisPoint=(self.LEN0/2.0, (DT[i]+DT[i+1])/2.0+i*self.MBH, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)

        a = mdb.models['Model-1'].rootAssembly
        myViewport.setValues(displayedObject=a)
        myViewport.view.setValues(session.views['Iso'])
        session.printToFile(fileName=MDL, format=PNG,canvasObjects=(myViewport,))
        myViewport.setValues(displayedObject=None)
        myViewport.assemblyDisplay.setValues(interactions=ON, 
            constraints=ON, connectors=ON, engineeringFeatures=ON)

        mask2 = '4'*(self.NUM0+1)+'0'*self.NUM0
        mask3 = '4'*self.NUM0
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
        p.seedPart(size=self.W, deviationFactor=0.1, minSizeFactor=0.1)
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
            p.seedPart(size=self.W, deviationFactor=0.1, minSizeFactor=0.1)
            p = mdb.models['Model-1'].parts[XZK[i]]
            p.generateMesh()

        a = mdb.models['Model-1'].rootAssembly
        a.regenerate()

        for JHD in MDS0:
            j = JHD[1]; k = JHD[2]
            JDMOD = 10**int((HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)%16//2)
            if (HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)%16%2 == 1:
                JDMOD = JDMOD*4
            a = mdb.models['Model-1'].rootAssembly
            if JHD[0] is 'M':
                mask4 = '[#0:' + str(int(HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)//16) + ' #' + str(JDMOD) + ' ]'
                n1 = a.instances['MB1-' + str(len(XZ))].nodes
            else:
                mask4 = '[#0:' + str(int(HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)//16) + ' #' + str(2*JDMOD) + ' ]'
                n1 = a.instances['MB1-0'].nodes
            nodes1 = n1.getSequenceFromMask(mask=(mask4, ), )
            a.Set(nodes=nodes1, name=JHD)

        myViewport.assemblyDisplay.setValues(mesh=OFF, adaptiveMeshConstraints=ON)
        myViewport.assemblyDisplay.meshOptions.setValues(meshTechnique=OFF)
        mdb.saveAs(pathName=WP + '\\' + MDL)
        TF.write('Output point = '+str(MDS0)+'\n')
        TF.close()

    def sbuilder(self):

        MDL = self.MDL0[:]
        FTX0 = {'1':0.275,'2':0.592,'3':1.0267,'4':1.778}
        FTX = FTX0
        for key, value in self.SFTX.items():
            FTX[key] = value

        myViewport = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200.0, height=100.0)
        myViewport.makeCurrent()
        myViewport.maximize()
        executeOnCaeStartup()
        myViewport.partDisplay.geometryOptions.setValues(referenceRepresentation=ON)

        WP = self.WP0 + '\\Smodel\\' + MDL
        if not os.path.exists(WP):
            os.makedirs(WP)
        os.chdir(WP)
        TF = open(self.MDL0+'.txt','wt')
        TF.write(self.MDL0+'\n\n\n'+'\nWork Path: '+WP+'\n')
        TF.write('self.STW = '+str(self.STW)+'\n'+'self.STH = '+str(self.STH)+'\n')
        TF.write('Cycle = '+str(self.T0)+'\n'+'Number = '+str(self.NUM0)+'\n'+'Length = '+str(self.LEN0)+'\n')
        TF.write('Grid Length = '+str(self.W)+'\n'+'Number = '+str(self.W1)+'\n'+'Panel Number = '+str(int(self.W1*16)**2)+'\n\n')

        FST = self.FST0[:]
        for i in range(len(self.FST0)):
            FST.append(self.FST0[-i-1])
        for i in range(len(FST)-1):
            FST[i+1] = FST[i+1] + FST[i]
        for i in range(len(FST)):
            FST[i] = FST[i]*self.T0/2.0/FST[-1]
        FST.insert(0,0)

        XZ = []; XZK = {}
        for i in range(len(MDL[1:])/3):
            XZ.append(MDL[3*i+1:3*i+3])
            XZK[MDL[3*i+1:3*i+3]] = 'S' + MDL[3*i+1:3*i+3]
        MX = sorted(set(XZ))

        Mdb()
        myViewport.setValues(displayedObject=None)

        s = mdb.models['Model-1'].ConstrainedSketch(name='MB1', sheetSize=200.0)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(0, -self.MBH), point2=(self.T0*self.NUM0, 0.0))
        p = mdb.models['Model-1'].Part(name='MB1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts['MB1']
        p.BaseSolidExtrude(sketch=s, depth=self.LEN0)
        s.unsetPrimaryObject()

        for i in MX:
            FSA = [(self.STW,FTX[i[0]]*self.T0/2.0*sin(pi*2.0*self.STW/self.T0-pi/2.0)/pi-XG[i[1]]/2.0)]
            FSB = []
            for j in FST:
                ls0 = FTX[i[0]]*cos(pi*j*2/self.T0-pi/2)
                lsa = j+XG[i[1]]*ls0/(ls0**2+1)**0.5/2.0
                if lsa>self.STW:
                    FSA.append((lsa,
        FTX[i[0]]*self.T0/2.0*sin(pi*j*2/self.T0-pi/2.0)/pi-XG[i[1]]/(ls0**2+1)**0.5/2.0))
                lsb = j-XG[i[1]]*ls0/(ls0**2+1)**0.5/2.0
                if lsb<self.T0/2.0-self.STW:
                    FSB.append((lsb,
        FTX[i[0]]*self.T0/2.0*sin(pi*j*2/self.T0-pi/2.0)/pi+XG[i[1]]/(ls0**2+1)**0.5/2.0))
            FSB.append((self.T0/2.0-self.STW,FTX[i[0]]*self.T0/2.0*sin(2*pi*(self.T0/2.0-self.STW)/self.T0-pi/2.0)/pi+XG[i[1]]/2.0))
            TF.write('FSA = '+str(FSA)+'\n'+'FSB = '+str(FSB)+'\n')

            h = FTX[i[0]]*self.T0/2.0/pi+XG[i[1]]/2.0+self.STH
            s = mdb.models['Model-1'].ConstrainedSketch(name=XZK[i], sheetSize=200.0)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=STANDALONE)
            s.Spline(points=FSA)
            s.Spline(points=FSB)
            s.Line(point1=(0, -h), point2=(self.STW, -h))
            s.HorizontalConstraint(entity=g[4], addUndoState=False)
            s.Line(point1=(self.STW, -h), point2=FSA[0])
            s.VerticalConstraint(entity=g[5], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
            s.Line(point1=FSB[-1], point2=(self.T0/2.0-self.STW, h))
            s.VerticalConstraint(entity=g[6], addUndoState=False)
            s.Line(point1=(self.T0/2.0-self.STW, h), point2=(self.T0/2.0, h))
            s.HorizontalConstraint(entity=g[7], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[6], entity2=g[7], addUndoState=False)
            s.Line(point1=(self.T0/2.0, h), point2=(self.T0/2.0, -h))
            s.VerticalConstraint(entity=g[8], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[7], entity2=g[8], addUndoState=False)
            s.copyMirror(mirrorLine=g[8], objectList=(g[2], g[3], g[4], g[5], g[6], g[7]))
            s.delete(objectList=(g[8], ))
            if self.NUM0>1:
                s.linearPattern(geomList=(g[2], g[3], g[4], g[5], g[6], g[7], g[9], g[10], 
        g[11], g[12], g[13], g[14]), vertexList=(), number1=self.NUM0, spacing1=self.T0, angle1=0.0, 
        number2=1, spacing2=20.0, angle2=90.0)
            s.Line(point1=(0, -h), point2=(0, -h+XG[i[1]]+self.STH))
            s.VerticalConstraint(entity=g[12*self.NUM0+3], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[4], entity2=g[12*self.NUM0+3], addUndoState=False)
            s.Line(point1=(self.T0*self.NUM0, -h), point2=(self.T0*self.NUM0, -h+XG[i[1]]+self.STH))
            s.VerticalConstraint(entity=g[12*self.NUM0+4], addUndoState=False)
            s.PerpendicularConstraint(entity1=g[12*self.NUM0-1], entity2=g[12*self.NUM0+4], addUndoState=False)
            p = mdb.models['Model-1'].Part(name=XZK[i], dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p = mdb.models['Model-1'].parts[XZK[i]]
            p.BaseSolidExtrude(sketch=s, depth=self.LEN0)
            s.unsetPrimaryObject()

        myViewport.partDisplay.setValues(sectionAssignments=ON, 
            engineeringFeatures=ON)
        myViewport.partDisplay.geometryOptions.setValues(
            referenceRepresentation=OFF)

        mdb.models['Model-1'].Material(name='FC')

        mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material='FC', thickness=None)

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

        MASK0 = 0x44
        h = 0x44
        for i in range(self.NUM0-1):
            h = h*32
            MASK0 = MASK0 + h
        MASK0 = hex(MASK0)[2:]
        if MASK0[-1] is 'L':
            MASK0 = MASK0[:-1]
        TF.write('Material Mask = '+MASK0+'\n')

        for i in MX:
            p = mdb.models['Model-1'].parts[XZK[i]]
            c = p.cells
            cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
            region = regionToolset.Region(cells=cells)
            p = mdb.models['Model-1'].parts[XZK[i]]
            s = p.faces
            side1Faces = s.getSequenceFromMask(mask=(dmask(MASK0), ), )
            normalAxisRegion = p.Surface(side1Faces=side1Faces, name='Surf-1')
            p = mdb.models['Model-1'].parts[XZK[i]]
            e = p.edges
            edges = e.getSequenceFromMask(mask=('[#0 #10 ]', ), )
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
            h = h + FTX[i[0]]*self.T0/pi + XG[i[1]] + 2.0*self.STH
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
            a.translate(instanceList=('MB1-' + str(i+1), ), vector=(0.0, DT[i+1]+(i+1)*self.MBH, 0.0))

        for i in range(len(XZ)):
            a = mdb.models['Model-1'].rootAssembly
            a.translate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), vector=(0.0, (DT[i]+DT[i+1])/2.0+i*self.MBH, 0.0))
            if MDL[3*i+3] is 'J':
                a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), 
                    axisPoint=(self.LEN0/2.0, (DT[i]+DT[i+1])/2.0+i*self.MBH, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)
            elif self.LEN0 == self.T0*self.NUM0:
                if MDL[3*i+3] is 'I':
                    a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), axisPoint=(self.LEN0/2.0, 0.0, self.LEN0/2.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
                elif MDL[3*i+3] is 'K':
                    a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), axisPoint=(self.LEN0/2.0, 0.0, self.LEN0/2.0), axisDirection=(0.0, 1.0, 0.0), angle=90.0)
                    a.rotate(instanceList=(XZK[XZ[i]] + '-' + str(i+1), ), 
                        axisPoint=(self.LEN0/2.0, (DT[i]+DT[i+1])/2.0+i*self.MBH, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=180.0)

        a = mdb.models['Model-1'].rootAssembly
        myViewport.setValues(displayedObject=a)
        myViewport.view.setValues(session.views['Iso'])
        session.printToFile(fileName=self.MDL0, format=PNG,canvasObjects=(myViewport,))
        myViewport.setValues(displayedObject=None)
        myViewport.assemblyDisplay.setValues(interactions=ON, 
            constraints=ON, connectors=ON, engineeringFeatures=ON)

        MASK01 = 1
        for i in range(self.NUM0-1):
            MASK01 = MASK01*2
            h = len(hex(MASK01))-2
            if hex(MASK01)[-1] is 'L':
                h = h-1
            MASK01 = MASK01+int(hex(MASK01)[2],16)*2*16**h
        MASK01 = hex(MASK01)[2:]
        if MASK01[-1] is 'L':
            MASK01 = MASK01[:-1]
        MASK01 = MASK01 + '0'*self.NUM0+'1'
        TF.write('Constraint Down Mask = '+MASK01+'\n')

        MASK02 = 0x10
        h = 0x10
        for i in range(self.NUM0-1):
            h = h*32
            MASK02 = MASK02 + h
        MASK02 = hex(MASK02)[2:]
        if MASK02[-1] is 'L':
            MASK02 = MASK02[:-1]
        TF.write('Constraint Up Mask = '+MASK02+'\n')

        for i in range(len(XZ)):
            a = mdb.models['Model-1'].rootAssembly
            s1 = a.instances['MB1-' + str(i)].faces
            side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
            region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-' + str(2*i+1))
            a = mdb.models['Model-1'].rootAssembly
            s1 = a.instances[XZK[XZ[i]] + '-' + str(i+1)].faces
            if MDL[3*i+3] in ['H','I']:
                side1Faces1 = s1.getSequenceFromMask(mask=(dmask(MASK01), ), )
            else:
                side1Faces1 = s1.getSequenceFromMask(mask=(dmask(MASK02), ), )
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
                side1Faces1 = s1.getSequenceFromMask(mask=(dmask(MASK02), ), )
            else:
                side1Faces1 = s1.getSequenceFromMask(mask=(dmask(MASK01), ), )
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
        p.seedPart(size=self.W, deviationFactor=0.1, minSizeFactor=0.1)
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
            p.seedPart(size=self.W, deviationFactor=0.1, minSizeFactor=0.1)
            p = mdb.models['Model-1'].parts[XZK[i]]
            p.generateMesh()

        a = mdb.models['Model-1'].rootAssembly
        a.regenerate()

        for JHD in MDS0:
            j = JHD[1]; k = JHD[2]
            JDMOD = 10**int((HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)%16//2)
            if (HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)%16%2 == 1:
                JDMOD = JDMOD*4
            a = mdb.models['Model-1'].rootAssembly
            if JHD[0] is 'M':
                MASK0 = '[#0:' + str(int(HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)//16) + ' #' + str(JDMOD) + ' ]'
                n1 = a.instances['MB1-' + str(len(XZ))].nodes
            else:
                MASK0 = '[#0:' + str(int(HG[j]*(self.W1*16.0+1.0)*self.W1+HG[k]*self.W1)//16) + ' #' + str(2*JDMOD) + ' ]'
                n1 = a.instances['MB1-0'].nodes
            nodes1 = n1.getSequenceFromMask(mask=(MASK0, ), )
            a.Set(nodes=nodes1, name=JHD)
        del JHD,JDMOD,MASK0
        myViewport.assemblyDisplay.setValues(mesh=OFF, adaptiveMeshConstraints=ON)
        myViewport.assemblyDisplay.meshOptions.setValues(meshTechnique=OFF)
        mdb.saveAs(pathName=WP + '\\' + self.MDL0)
        TF.write('Output point = '+str(MDS0)+'\n')
        TF.close()