import os
import numpy as np
import struct
#%%

bool_str = lambda x: '.true.' if x else '.false.'

class EnvironmentIn():
    def __init__(self,nbSourceContainers =1,nbObserverContainers = 1, pressureFolderName='pressure/',SPLFolderName= 'spl/',sigmaFolderName= 'sigma/',debugLevel= 1,ASCIIOutputFlag= False,OASPLdBFlag= False,OASPLdBAFlag= False,spectrumFlag= False,
                 SPLdBFlag= False,SPLdBAFlag= False,pressureGradient1AFlag= False, acousticPressureFlag= True,thicknessNoiseFlag= False,loadingNoiseFlag= True,totalNoiseFlag= False,
                 sigmaFlag= False,loadingNoiseSigmaFlag= True,thicknessNoiseSigmaFlag= True,totalNoiseSigmaFlag= True,normalSigmaFlag= True,machSigmaFlag= True,observerSigmaFlag= True,velocitySigmaFlag= True,
                 accelerationSigmaFlag= True,densitySigmaFlag= True,momentumSigmaFlag= True,pressureSigmaFlag= True,loadingSigmaFlag= True,areaSigmaFlag= True,MdotrSigmaFlag= True,iblankSigmaFlag = True):
        self.nbSourceContainers = nbSourceContainers
        self.nbObserverContainers = nbObserverContainers
        self.pressureFolderName = pressureFolderName
        self.SPLFolderName = SPLFolderName
        self.sigmaFolderName = sigmaFolderName
        self.debugLevel = debugLevel
        self.ASCIIOutputFlag = ASCIIOutputFlag
        self.OASPLdBFlag = OASPLdBFlag
        self.OASPLdBAFlag = OASPLdBAFlag
        self.spectrumFlag = spectrumFlag
        self.SPLdBFlag = SPLdBFlag
        self.SPLdBAFlag = SPLdBAFlag
        self.pressureGradient1AFlag = pressureGradient1AFlag
        self.acousticPressureFlag = acousticPressureFlag
        self.thicknessNoiseFlag = thicknessNoiseFlag
        self.loadingNoiseFlag = loadingNoiseFlag
        self.totalNoiseFlag = totalNoiseFlag
        self.sigmaFlag = sigmaFlag
        self.loadingNoiseSigmaFlag = loadingNoiseSigmaFlag
        self.thicknessNoiseSigmaFlag = thicknessNoiseSigmaFlag
        self.totalNoiseSigmaFlag = totalNoiseSigmaFlag
        self.normalSigmaFlag = normalSigmaFlag
        self.machSigmaFlag = machSigmaFlag
        self.observerSigmaFlag = observerSigmaFlag
        self.velocitySigmaFlag = velocitySigmaFlag
        self.accelerationSigmaFlag = accelerationSigmaFlag
        self.densitySigmaFlag = densitySigmaFlag
        self.momentumSigmaFlag = momentumSigmaFlag
        self.pressureSigmaFlag = pressureSigmaFlag
        self.loadingSigmaFlag = loadingSigmaFlag
        self.areaSigmaFlag = areaSigmaFlag
        self.MdotrSigmaFlag = MdotrSigmaFlag
        self.iblankSigmaFlag = iblankSigmaFlag

class EnvironmentConstants():
    def __init__(self,rho = 1.225,c = 342,gamma = 1.4,mu = 1.51e-5,nu  = 1.51e-5/1.225,P_ref = 2e-5,RelHumidity = 70,gasConstant = 287.04):
        
        self.rho = rho
        self.c = c
        self.gamma = gamma
        self.mu = mu
        self.nu = nu
        self.P_ref = P_ref
        self.RelHumidity = RelHumidity
        self.gasConstant = gasConstant

class ObserverIn():
    def __init__(self,nt,Title = 'observer',tMin = None,tMax = None,nbBase = 0,attachedTo = None,nbBaseObsContFrame = 0,nbBaseLocalFrame = 0,fileName = None,
                highPassFrequency = None,lowPassFrequency= None,xLoc = None,yLoc = None ,zLoc= None,nbx= None,nby= None,nbz= None,xMin= None,xMax= None,yMin= None,yMax= None,zMin= None,zMax= None,
                radius = None, nbTheta = None,nbPsi = None,thetaMin = None,thetaMax = None,psiMin = None,psiMax = None,indexSwap = False):
        
        self.nt = nt
        self.Title = Title
        if tMin is not None:
            self.tMin = tMin
            self.tMax = tMax
        self.nbBase = nbBase
        if attachedTo is not None:
            self.attachedTo = attachedTo
            self.nbBaseObsContFrame = nbBaseObsContFrame
            self.nbBaseLocalFrame = nbBaseLocalFrame
        if fileName is not None:
            self.fileName = fileName
        if highPassFrequency is not None:
            self.highPassFrequency = highPassFrequency
        if lowPassFrequency is not None:
            self.lowPassFrequency = lowPassFrequency
        
        if xLoc is not None:
            self.xLoc = xLoc
            self.yLoc = yLoc
            self.zLoc = zLoc

        if nbx is not None:
            self.nbx = nbx
            self.nby = nby
            self.nbz = nbz
            self.xMin = xMin
            self.xMax = xMax
            self.yMin = yMin
            self.yMax = yMax
            self.zMin = zMin
            self.zMax = zMax

        if radius is not None:
            self.radius = radius
            self.nbTheta = nbTheta
            self.nbPsi = nbPsi
            self.thetaMin = thetaMin
            self.thetaMax = thetaMax
            self.psiMin = psiMin
            self.psiMax = psiMax
            self.indexSwap = indexSwap


class ContainerIn():
    def __init__(self,Title,nbContainer = 0,nbBase = 0,patchGeometryFile = None,patchLoadingFile = None,PeggNoiseFlag = False,BPMNoiseFlag = False,dtau = None,periodicKeyOffset = None):
        
        self.Title = Title
        self.nbContainer = nbContainer
        self.nbBase = nbBase
        if patchGeometryFile is not None:
            self.patchGeometryFile = patchGeometryFile
        if patchLoadingFile is not None:
            self.patchLoadingFile = patchLoadingFile
        if Title == 'rotor':
            self.PeggNoiseFlag = PeggNoiseFlag
            self.BPMNoiseFlag = BPMNoiseFlag
        if dtau is not None:
            self.dtau = dtau
        if periodicKeyOffset is not None:
            self.periodicKeyOffset = periodicKeyOffset

class CB():
    def __init__(self,Title = 'cb',Rotation = False,AxisType='TimeIndependent',TranslationType = 'TimeIndependent',AngleType = 'TimeIndependent',AH = None,VH = None, Y0 = None, 
                 Omega = 0,Psi0 = 0, AxisValue = [0,0,1], TranslationValue =None,AngleValue = 0):
        
        self.Title = Title
        self.Rotation = Rotation

        if Rotation or AngleValue is not None :
            self.AxisType = AxisType
            self.AngleType = AngleType
            self.Omega = Omega
            self.Psi0 = Psi0
            self.AxisValue = AxisValue
            self.AngleValue = AngleValue

        if VH is not None:
            self.TranslationType = TranslationType
            self.AH = AH
            self.VH = VH
            self.Y0 = Y0
            self.TranslationValue = TranslationValue

class caseName():
    def __init__(self,globalFolderName = './',caseNameFile = 'case.nam'):
        self.globalFolderName = globalFolderName+'/'
        self.caseNameFile = caseNameFile

def write_nml(file,nml):
    file.write(f'&{type(nml).__name__}\n')
    for k,v in vars(nml).items():
        if type(v) is bool:
            file.write(f"\t{k}={bool_str(v)}\n")
        elif type(v) is str:
            file.write(f"\t{k}='{v}'\n")
        elif type(v) is list:
            file.write(f"\t{k}={str(v)[1:-1]}\n")
        else:
            file.write(f"\t{k}={v}\n")
    file.write("/\n\n")

def write_nml_file(file_path,nml):
    with open(file_path, 'w') as f:
        for nml_iter in nml:
            write_nml(f,nml_iter)

def write_observer_file(file_path, observer):
    # observer should be an array of coordinates of shape (N_observer x 3)
    # file needs to be in ascii format
    iMax = int(len(observer))
    jMax = 1
    kMax = 1
    with open(file_path,'w') as f_ascii:
        f_ascii.write(f'{iMax} \t {jMax} \t {kMax} \n')
        for i in range(3):
            f_ascii.write(f'{str(observer[:,i])[1:-1]}\n')



def constant_compact_loading_write(file_path, loads,ascii = False):
    
    # loads: array of shape [jMax x 3]

    magic_number = 42                #4-byte signed
    version_number = [1,0]
    comments = "Compact loading file"
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 1                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    data_type = 2                    # data is surface pressure (1), surface loading vector (2), flow parameters (3)
    ref_frame = 3                    # reference frame is a stationary (1), rotating ground-fixed frame (2), patch-fixed frame (3). Note that this has no effect on pressure data, and that “2” and “3” are equivalent for load vectors.
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    dataZones = [1,-1]               # number of zones with data, zone designation (negative to skip thickness calc).

    zoneName = "LiftingLine"
    iMax = 1                         # number of chordwise elements
    jMax = loads.shape[0]                    # number of spanwise elements
    
    with open(file_path,'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        comments_bin = struct.pack(f'<{len(comments)}s', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', 2))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geom_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', data_type))
        f_bin.write(struct.pack('<i', ref_frame))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<i', 0))
        [f_bin.write(struct.pack('<i', zone)) for zone in dataZones]
        f_bin.write(struct.pack('<32s', bytes(zoneName, encoding='ascii')))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))
        for i in range(3):
            f_bin.write(struct.pack(f'<{jMax}f', *loads[:,i]))

    if ascii:
        with open(file_path[:-4]+'.ascii','w') as f_ascii:
            f_ascii.write(f'{magic_number}\n')
            f_ascii.write(f'{version_number[0]}\n')
            f_ascii.write(f'{version_number[1]}\n')
            f_ascii.write(f'{comments}\n')
            f_ascii.write(f'{2}\n')
            f_ascii.write(f'{Nzones}\n')
            f_ascii.write(f'{grid_type}\n')
            f_ascii.write(f'{geom_type}\n')
            f_ascii.write(f'{vector_centering}\n')
            f_ascii.write(f'{data_type}\n')
            f_ascii.write(f'{ref_frame}\n')
            f_ascii.write(f'{precision}\n')
            f_ascii.write(f'{0}\n')
            f_ascii.write(f'{0}\n')
            [f_ascii.write(f'{zone}\n') for zone in dataZones]
            f_ascii.write(f'{zoneName}\n')
            f_ascii.write(f'{iMax}\n')
            f_ascii.write(f'{jMax}\n')
            for i in range(3):
                f_ascii.write(f'{loads[:,i]}\n')

def periodic_compact_loading_write(file_path,keys,period, loads,ascii = False):
    
    # loads: array of shape [N_steps x jMax x 3]

    magic_number = 42                #4-byte signed
    version_number = [1,0]
    comments = "Compact loading file"
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 2                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    data_type = 2                    # data is surface pressure (1), surface loading vector (2), flow parameters (3)
    ref_frame = 3                    # reference frame is a stationary (1), rotating ground-fixed frame (2), patch-fixed frame (3). Note that this has no effect on pressure data, and that “2” and “3” are equivalent for load vectors.
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    dataZones = [1,-1]               # number of zones with data, zone designation (negative to skip thickness calc).

    zoneName = "LiftingLine"
    N_steps = len(keys)            # number of time steps [sec]
    iMax = 1                         # number of chordwise elements
    jMax = loads.shape[1]                    # number of spanwise elements
    
    with open(file_path,'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        comments_bin = struct.pack(f'<{len(comments)}s', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', 2))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geom_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', data_type))
        f_bin.write(struct.pack('<i', ref_frame))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<i', 0))
        [f_bin.write(struct.pack('<i', zone)) for zone in dataZones]
        f_bin.write(struct.pack('<32s', bytes(zoneName, encoding='ascii')))
        f_bin.write(struct.pack('<f', period))
        f_bin.write(struct.pack('<i', N_steps))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))

        for i,k_iter in enumerate(keys):
            f_bin.write(struct.pack('<f', k_iter))
            for ii in range(3):
                f_bin.write(struct.pack(f'<{jMax}f', *loads[i][:,ii]))
    if ascii:
        with open(file_path[:-4]+'.ascii','w') as f_ascii:
            f_ascii.write(f'{magic_number}\n')
            f_ascii.write(f'{version_number[0]}\n')
            f_ascii.write(f'{version_number[1]}\n')
            f_ascii.write(f'{comments}\n')
            f_ascii.write(f'{2}\n')
            f_ascii.write(f'{Nzones}\n')
            f_ascii.write(f'{grid_type}\n')
            f_ascii.write(f'{geom_type}\n')
            f_ascii.write(f'{vector_centering}\n')
            f_ascii.write(f'{data_type}\n')
            f_ascii.write(f'{ref_frame}\n')
            f_ascii.write(f'{precision}\n')
            f_ascii.write(f'{0}\n')
            f_ascii.write(f'{0}\n')
            [f_ascii.write(f'{zone}\n') for zone in dataZones]
            f_ascii.write(f'{zoneName}\n')
            f_ascii.write(f'{period}\n')
            f_ascii.write(f'{N_steps}\n')
            f_ascii.write(f'{iMax}\n')
            f_ascii.write(f'{jMax}\n')
            for i,k_iter in enumerate(keys):
                f_ascii.write(f'{k_iter}\n')
                for ii in range(3):
                    f_ascii.write(f'{loads[i][:,ii]}\n')

def aperiodic_compact_loading_write(file_path,t, loads,ascii = False):
    
    # loads: array of shape [N_steps x jMax x 3]

    magic_number = 42                #4-byte signed
    version_number = [1,0]
    comments = "Compact loading file"
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 3                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    data_type = 2                    # data is surface pressure (1), surface loading vector (2), flow parameters (3)
    ref_frame = 3                    # reference frame is a stationary (1), rotating ground-fixed frame (2), patch-fixed frame (3). Note that this has no effect on pressure data, and that “2” and “3” are equivalent for load vectors.
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    dataZones = [1,-1]               # number of zones with data, zone designation (negative to skip thickness calc).

    zoneName = "Lifting Line"
    N_steps = loads.shape[0]            # number of time steps [sec]
    iMax = 1                         # number of chordwise elements
    jMax = loads.shape[1]                    # number of spanwise elements
    
    with open(file_path,'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        comments_bin = struct.pack(f'<{len(comments)}s', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', 2))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geom_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', data_type))
        f_bin.write(struct.pack('<i', ref_frame))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<i', 0))
        [f_bin.write(struct.pack('<i', zone)) for zone in dataZones]
        f_bin.write(struct.pack('<32s', bytes(zoneName, encoding='ascii')))
        f_bin.write(struct.pack('<i', N_steps))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))

        for i,t_iter in enumerate(t):
            f_bin.write(struct.pack('<f', t_iter))
            for ii in range(3):
                f_bin.write(struct.pack(f'<{jMax}f', *loads[i][:,ii]))
    if ascii:
        with open(file_path[:-4]+'.ascii','w') as f_ascii:
            f_ascii.write(f'{magic_number}\n')
            f_ascii.write(f'{version_number[0]}\n')
            f_ascii.write(f'{version_number[1]}\n')
            f_ascii.write(f'{comments}\n')
            f_ascii.write(f'{2}\n')
            f_ascii.write(f'{Nzones}\n')
            f_ascii.write(f'{grid_type}\n')
            f_ascii.write(f'{geom_type}\n')
            f_ascii.write(f'{vector_centering}\n')
            f_ascii.write(f'{data_type}\n')
            f_ascii.write(f'{ref_frame}\n')
            f_ascii.write(f'{precision}\n')
            f_ascii.write(f'{0}\n')
            f_ascii.write(f'{0}\n')
            [f_ascii.write(f'{zone}\n') for zone in dataZones]
            f_ascii.write(f'{zoneName}\n')
            f_ascii.write(f'{N_steps}\n')
            f_ascii.write(f'{iMax}\n')
            f_ascii.write(f'{jMax}\n')
            for i,t_iter in enumerate(t):
                f_ascii.write(f'{t_iter}\n')
                for ii in range(3):
                    f_ascii.write(f'{loads[i][:,ii]}\n')

def constant_compact_geometry_write(file_path, nodes,norms,ascii = False):

    # nodes: array of lifting line coordinates where loading verctors are perscribed [iMax x jMax x 3]
    # norms: array of surface normal vetors of the lifting line [iMaxx jMax x 3]

    magic_number = 42                #4-byte signed
    version_number = [1,0]
    units = 'Pa'  
    comments = "Compact geometry patch file"
    geometryFile = 1                 # If this is a geometry file = 1, -1 for psuedo geometry
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 1                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    iblank = 0                       # iblank values are included (1) or not included (0)
    zoneName = "Lifting line"
    iMax = nodes.shape[0]             # number of chordwise elements
    jMax = nodes.shape[1]             # number of spanwise elements
    
    with open(file_path,'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        f_bin.write(struct.pack('<32s', bytes(units, encoding='ascii')))
        comments_bin = struct.pack('<' + str(len(comments)) + 's', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', geometryFile))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geom_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', iblank))
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<32s', bytes(zoneName, encoding='ascii')))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))
        for ii in range(3):
            f_bin.write(struct.pack(f'<{int(iMax*jMax)}f', *nodes[:,:,ii].flatten(order = 'F')))
        for ii in range(3):
            f_bin.write(struct.pack(f'<{int(iMax*jMax)}f', *norms[:,:,ii].flatten(order = 'F')))

    if ascii:
        with open(file_path[:-4]+'.ascii','w') as f_ascii:
            f_ascii.write(f'{magic_number}\n')
            f_ascii.write(f'{version_number[0]}\n')
            f_ascii.write(f'{version_number[1]}\n')
            f_ascii.write(f'{comments}\n')
            f_ascii.write(f'{geometryFile}\n')
            f_ascii.write(f'{Nzones}\n')
            f_ascii.write(f'{grid_type}\n')
            f_ascii.write(f'{geom_type}\n')
            f_ascii.write(f'{vector_centering}\n')
            f_ascii.write(f'{precision}\n')
            f_ascii.write(f'{iblank}\n')
            f_ascii.write(f'{0}\n')
            f_ascii.write(f'{zoneName}\n')
            f_ascii.write(f'{iMax}\n')
            f_ascii.write(f'{jMax}\n')
            for ii in range(3):
                f_ascii.write(f'{nodes[:,:,ii].squeeze()}\n')
            for ii in range(3):
                f_ascii.write(f'{norms[:,:,ii].squeeze()}\n')


def aperiodic_compact_geometry_write(file_path,t, nodes,norms,ascii = False):

    # nodes: array of lifting line coordinates where loading verctors are perscribed [N_steps x jMax x 3]
    # norms: array of surface normal vetors of the lifting line [N_steps x jMax x 3]

    magic_number = 42                #4-byte signed
    version_number = [1,0]
    units = 'Pa'  
    comments = "Compact geometry patch file"
    geometryFile = 1                 # If this is a geometry file = 1, -1 for psuedo geometry
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 3                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    iblank = 0                       # iblank values are included (1) or not included (0)
    zoneName = "Lifting line"
    N_steps = nodes.shape[0]          # number of time steps [sec]
    iMax = 1                         # number of chordwise elements
    jMax = nodes.shape[1]             # number of spanwise elements
    
    with open(file_path,'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        f_bin.write(struct.pack('<32s', bytes(units, encoding='ascii')))
        comments_bin = struct.pack('<' + str(len(comments)) + 's', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', geometryFile))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geom_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', iblank))
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<32s', bytes(zoneName, encoding='ascii')))
        f_bin.write(struct.pack('<i', N_steps))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))
        for i,t_iter in enumerate(t):
            f_bin.write(struct.pack('<f', t_iter))
            for ii in range(3):
                f_bin.write(struct.pack(f'<{jMax}f', *nodes[i][:,ii]))
            for ii in range(3):
                f_bin.write(struct.pack(f'<{jMax}f', *norms[i][:,ii]))

    if ascii:
        with open(file_path[:-4]+'.ascii','w') as f_ascii:
            f_ascii.write(f'{magic_number}\n')
            f_ascii.write(f'{version_number[0]}\n')
            f_ascii.write(f'{version_number[1]}\n')
            f_ascii.write(f'{comments}\n')
            f_ascii.write(f'{geometryFile}\n')
            f_ascii.write(f'{Nzones}\n')
            f_ascii.write(f'{grid_type}\n')
            f_ascii.write(f'{geom_type}\n')
            f_ascii.write(f'{vector_centering}\n')
            f_ascii.write(f'{precision}\n')
            f_ascii.write(f'{iblank}\n')
            f_ascii.write(f'{0}\n')
            f_ascii.write(f'{zoneName}\n')
            f_ascii.write(f'{N_steps}\n')
            f_ascii.write(f'{iMax}\n')
            f_ascii.write(f'{jMax}\n')
            for i,t_iter in enumerate(t):
                f_ascii.write(f'{t_iter}\n')
                for ii in range(3):
                    f_ascii.write(f'{nodes[i][:,ii]}\n')
                for ii in range(3):
                    f_ascii.write(f'{norms[i][:,ii]}\n')


