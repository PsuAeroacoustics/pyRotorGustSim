
from wopwop_input_generator import *
from AnalyzeDegenGeom import *
from ProcessGeom import *

#%%
def wopwop_input_configure(geom_params,input_params,res_param,observer_params,acs_params,saved_params):

    # initializes namelist and environmental object
    nml = []
    environment_in = EnvironmentIn(debugLevel = 4,ASCIIOutputFlag=True,totalNoiseFlag=True,pressureFolderName = '/')    
    
    environment_in = EnvironmentIn(
        		pressureFolderName = acs_params['pressureFolderName'],SPLFolderName =acs_params['SPLFolderName'], sigmaFolderName=acs_params['sigmaFolderName'], debugLevel=acs_params['debugLevel'], ASCIIOutputFlag=acs_params['ASCIIOutputFlag'],
                  OASPLdBFlag=acs_params['OASPLdBFlag'], OASPLdBAFlag=acs_params['OASPLdBAFlag'], spectrumFlag=acs_params['spectrumFlag'],SPLdBFlag= acs_params['SPLdBFlag'], SPLdBAFlag=acs_params['SPLdBAFlag'], pressureGradient1AFlag=acs_params['pressureGradient1AFlag'],
                    acousticPressureFlag=acs_params['acousticPressureFlag'], thicknessNoiseFlag=acs_params['thicknessNoiseFlag'], loadingNoiseFlag=acs_params['loadingNoiseFlag'], totalNoiseFlag=acs_params['totalNoiseFlag'], sigmaFlag=acs_params['sigmaFlag'], loadingNoiseSigmaFlag=acs_params['loadingNoiseSigmaFlag'],
                      thicknessNoiseSigmaFlag=acs_params['thicknessNoiseSigmaFlag'], totalNoiseSigmaFlag=acs_params['totalNoiseSigmaFlag'],normalSigmaFlag= acs_params['normalSigmaFlag'], machSigmaFlag=acs_params['machSigmaFlag'], observerSigmaFlag=acs_params['observerSigmaFlag'], velocitySigmaFlag=acs_params['velocitySigmaFlag'],
                        accelerationSigmaFlag=acs_params['accelerationSigmaFlag'], densitySigmaFlag=acs_params['densitySigmaFlag'], momentumSigmaFlag=acs_params['momentumSigmaFlag'], pressureSigmaFlag=acs_params['pressureSigmaFlag'], loadingSigmaFlag=acs_params['loadingSigmaFlag'],areaSigmaFlag= acs_params['areaSigmaFlag'],MdotrSigmaFlag=acs_params['MdotrSigmaFlag'], iblankSigmaFlag=acs_params['iblankSigmaFlag'])
    
    environment_const = EnvironmentConstants()

    #   Determines whether to write out a observer file or not based on what is provided in the corresponding JSON file
    if not 'nbTheta' in observer_params or 'nbx' in observer_params or 'xLoc' in observer_params:

        radius = observer_params['radius']
        theta = np.array(observer_params['theta'])*np.pi/180
        phi = observer_params['phi']*np.pi/180
        observer_coordinates = np.array([radius*np.cos(theta),radius*np.sin(theta),radius*np.sin(phi)*np.ones(len(theta))]).T
        write_observer_file(os.path.join(saved_params['acs_dir'],f'observer.ascii'),observer=observer_coordinates)
        observer_in = ObserverIn(nt = saved_params['iterations'],Title='mic array',attachedTo = 'aircraft',tMin = saved_params['t'][int(saved_params['iterations']/2)],tMax=saved_params['t'][-1],fileName='observer.ascii',highPassFrequency=observer_params['highPassFrequency'])

    elif observer_params['xLoc'] is not None:
        observer_in = ObserverIn(nt = saved_params['iterations'],Title='mic array',attachedTo = 'aircraft',tMin = saved_params['t'][int(saved_params['iterations']/2)],tMax=saved_params['t'][-1],xLoc=observer_params['xLoc'],yLoc=observer_params['yLoc'],zLoc=observer_params['zLoc'],highPassFrequency = observer_params['highPassFrequency'])

    elif observer_params['nbTheta'] is not None:
        observer_in = ObserverIn(nt = saved_params['iterations'],Title='mic array',attachedTo = 'aircraft',tMin = saved_params['t'][int(saved_params['iterations']/2)],tMax=saved_params['t'][-1],radius=observer_params['radius'],nbTheta=observer_params['nbTheta'],nbPsi=observer_params['nbPsi'],thetaMin=observer_params['thetaMin'],thetaMax=observer_params['thetaMax'],psiMin=observer_params['psiMin'],psiMax = observer_params['psiMax'],highPassFrequency = observer_params['highPassFrequency'])

    elif observer_params['nbx'] is not None:
        observer_in = ObserverIn(nt = saved_params['iterations'],Title='mic array',attachedTo = 'aircraft',tMin = saved_params['t'][int(saved_params['iterations']/2)],tMax=saved_params['t'][-1],nbx=observer_params['nbx'],xMin=observer_params['xMin'],xMax=observer_params['xMax'],nby=observer_params['nby'],yMin=observer_params['yMin'],yMax=observer_params['yMax'],nbz = observer_params['nbz'],zMin = observer_params['zMin'],zMax=observer_params['zMax'],highPassFrequency = observer_params['highPassFrequency'])

    
    aircraft_container = ContainerIn(Title='aircraft',nbContainer = 1)
    rotor_container = ContainerIn(Title='rotor',nbContainer = geom_params['number_of_blades'],nbBase=1)
    rotor_cb = CB(Title='rotation',Rotation = True,AngleType='KnownFunction',Omega=input_params['flight_params']['omega'],AxisValue=[0,0,1])
    nml.extend([environment_in,environment_const,observer_in,aircraft_container,rotor_container,rotor_cb])

    for b_iter in range(geom_params['number_of_blades']):
        nml.append(ContainerIn(Title=f'blade {b_iter} loading',nbBase=1,patchGeometryFile = f'lifting_line_geometry.dat',patchLoadingFile = f'loading_blade_{b_iter}.dat',dtau = saved_params['t'][-1]/(observer_params['nt']-1)))
        nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/geom_params['number_of_blades']*b_iter))
        # nml.append(ContainerIn(Title=f'blade {b_iter} thickness',nbBase=3,patchGeometryFile = f'blade_geometry.dat',dtau = saved_params['t'][-1]/(observer_params['nt']-1)))
        # nml.append(CB(Title=f'blade {b_iter} azimuthal offset',AxisValue=[0,0,1],AngleValue=2*np.pi/geom_params['number_of_blades']*b_iter))
        # nml.append(CB(Title=f'blade {b_iter} th0',AxisValue=[1,0,0],AngleValue=saved_params['th0']))
        # nml.append(CB(Title=f'Aligns blade span with positive x-direction',AxisValue=[0,0,-1],AngleValue=np.pi/2))

    write_nml_file(os.path.join(saved_params['acs_dir'],f"{input_params['case_name']}.nam"),nml)
    case = caseName(globalFolderName = saved_params['acs_dir'],caseNameFile=f"{input_params['case_name']}.nam")
    write_nml_file(os.path.join(saved_params['case_dir'],'cases.nam'),[case])
    constant_compact_geometry_write(os.path.join(saved_params['acs_dir'],f'lifting_line_geometry.dat'),nodes=saved_params['lifting_line_nodes'],norms=saved_params['lifting_line_norms'],ascii=False)

    for b_iter in range(geom_params['number_of_blades']):
        aperiodic_compact_loading_write(os.path.join(saved_params['acs_dir'],f'loading_blade_{b_iter}.dat'),t = saved_params['t'], loads = saved_params['loads'],ascii = False)

