# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

import gmsh
import fmt_scen as scen  
from msh_to_xdmf import msh_to_xdmf 
import numpy as np
from math import sqrt

def create_team_problem_32_msh(modelName: str, sTags, lTags, meshFactor):
    gmsh.initialize()
    #gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.model.add(modelName)

    dim = 2
    meshFine = meshFactor * 2e-3
    meshCoarse  =  meshFactor * 2e-3
    meshCoil    = meshFactor * 2e-3
    meshBoundary   = meshFactor * 1.5e-2
    meshCoilInner  = meshFactor * 3e-3

    ironWidth = 30e-3
    airWidth = 42.25e-3
    coreHeight = 60e-3
    distCoil = 2e-3
    coilWidth = 5e-3
    coilHeightFak = 0.95
    boundDist = 2*coreHeight
    domainRadius = (3*ironWidth + 2*airWidth + 2*boundDist) / 2

    coilArea = 2*coreHeight*coilHeightFak*coilWidth
    turnNo = 90
    currScale  = turnNo / coilArea

    pSym = []
    pBoundCore = []
    # Points for core
    pCore = []
    pCore.append(gmsh.model.geo.add_point(0, 0, 0, meshFine)) # 0
    pBoundCore.append(pCore[len(pCore)-1])
    pCore.append(gmsh.model.geo.add_point(ironWidth               ,0                        , 0, meshFine)) # 1
    pBoundCore.append(pCore[len(pCore)-1])
    pCore.append(gmsh.model.geo.add_point(ironWidth               ,coreHeight               , 0, meshFine)) # 2
    pCore.append(gmsh.model.geo.add_point(ironWidth + airWidth    ,coreHeight               , 0, meshFine)) # 3
    pCore.append(gmsh.model.geo.add_point(ironWidth + airWidth    ,0                        , 0, meshFine)) # 4
    pBoundCore.append(pCore[len(pCore)-1])
    pCore.append(gmsh.model.geo.add_point(2*ironWidth + airWidth  ,0                        , 0, meshFine)) # 5
    pBoundCore.append(pCore[len(pCore)-1])
    pCore.append(gmsh.model.geo.add_point(2*ironWidth + airWidth  ,coreHeight               , 0, meshFine)) # 6
    pCore.append(gmsh.model.geo.add_point(2*ironWidth + 2*airWidth,coreHeight               , 0, meshFine)) # 7
    pCore.append(gmsh.model.geo.add_point(2*ironWidth + 2*airWidth,0                        , 0, meshFine)) # 8
    pBoundCore.append(pCore[len(pCore)-1])
    pCore.append(gmsh.model.geo.add_point(3*ironWidth + 2*airWidth,0                        , 0, meshFine)) # 9
    pBoundCore.append(pCore[len(pCore)-1])
    pCore.append(gmsh.model.geo.add_point(3*ironWidth + 2*airWidth, coreHeight + ironWidth  , 0, meshFine)) # 10
    pCore.append(gmsh.model.geo.add_point(0                       , coreHeight + ironWidth  , 0, meshFine)) # 11

    lCore = []
    lBoundCore = []
    lCore.append(gmsh.model.geo.add_line(pCore[0],pCore[1])) # 0
    lBoundCore.append(lCore[len(lCore)-1])
    lCore.append(gmsh.model.geo.add_line(pCore[1],pCore[2])) # 1 -> start Vac
    lCore.append(gmsh.model.geo.add_line(pCore[2],pCore[3])) # 2
    lCore.append(gmsh.model.geo.add_line(pCore[3],pCore[4])) # 3
    lCore.append(gmsh.model.geo.add_line(pCore[4],pCore[5])) # 4
    lBoundCore.append(lCore[len(lCore)-1])
    lCore.append(gmsh.model.geo.add_line(pCore[5],pCore[6])) # 5
    lCore.append(gmsh.model.geo.add_line(pCore[6],pCore[7])) # 6
    lCore.append(gmsh.model.geo.add_line(pCore[7],pCore[8])) # 7 
    lCore.append(gmsh.model.geo.add_line(pCore[8],pCore[9])) # 8
    lBoundCore.append(lCore[len(lCore)-1])
    lCore.append(gmsh.model.geo.add_line(pCore[9],pCore[10])) # 9
    lCore.append(gmsh.model.geo.add_line(pCore[10],pCore[11])) # 10
    lCore.append(gmsh.model.geo.add_line(pCore[11],pCore[0])) # 11

    coreLoop = gmsh.model.geo.add_curve_loop(lCore)
    core = gmsh.model.geo.add_plane_surface([coreLoop])

    # coil 1
    origCoil = -distCoil-coilWidth
    pBoundCoil1 = []
    
    pCoil = []
    pCoil.append(gmsh.model.geo.add_point(origCoil, 0,0, meshFine))
    pBoundCoil1.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, 0,0, meshFine))
    pBoundCoil1.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, coreHeight *coilHeightFak ,0, meshFine))
    pCoil.append(gmsh.model.geo.add_point(origCoil,coreHeight *coilHeightFak ,0, meshFine))

    lCoil = []
    lBoundCoil1 = []
    for i in range(0, len(pCoil)):
        lCoil.append(gmsh.model.geo.add_line(pCoil[i], pCoil[(i+1)%len(pCoil)]))
        if i == 0:
            lBoundCoil1.append(lCoil[len(lCoil)-1])

    coilLoop1 = gmsh.model.geo.add_curve_loop(lCoil)

    coil1Pos= gmsh.model.geo.add_plane_surface([coilLoop1])

    # coil 2
    origCoil = distCoil + ironWidth
    pCoil = []
    pBoundCoil2 = []
    pCoil.append(gmsh.model.geo.add_point(origCoil, 0,0, meshFine))
    pBoundCoil2.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, 0,0, meshFine))
    pBoundCoil2.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, coreHeight *coilHeightFak ,0, meshFine))
    pCoil.append(gmsh.model.geo.add_point(origCoil,coreHeight *coilHeightFak ,0, meshFine))

    lCoil = []
    lBoundCoil2 = []
    for i in range(0, len(pCoil)):
        lCoil.append(gmsh.model.geo.add_line(pCoil[i], pCoil[(i+1)%len(pCoil)]))
        if i == 0:
            lBoundCoil2.append(lCoil[len(lCoil)-1])

    coilLoop2 = gmsh.model.geo.add_curve_loop(lCoil)

    coil1Neg= gmsh.model.geo.add_plane_surface([coilLoop2])

    # coil 3
    origCoil =  2*ironWidth + 2*airWidth -distCoil-coilWidth
    pCoil = []
    pBoundCoil3 = []
    pCoil.append(gmsh.model.geo.add_point(origCoil, 0,0, meshFine))
    pBoundCoil3.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, 0,0, meshFine))
    pBoundCoil3.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, coreHeight *coilHeightFak ,0, meshFine))
    pCoil.append(gmsh.model.geo.add_point(origCoil,coreHeight *coilHeightFak ,0, meshFine))

    lCoil = []
    lBoundCoil3 = []
    for i in range(0, len(pCoil)):
        lCoil.append(gmsh.model.geo.add_line(pCoil[i], pCoil[(i+1)%len(pCoil)]))
        if i == 0:
            lBoundCoil3.append(lCoil[len(lCoil)-1])

    coilLoop3 = gmsh.model.geo.add_curve_loop(lCoil)

    coil2Pos= gmsh.model.geo.add_plane_surface([coilLoop3])

    # coil 4
    origCoil = 3*ironWidth+ 2*airWidth + distCoil
    pCoil = []
    pBoundCoil4 = []
    pCoil.append(gmsh.model.geo.add_point(origCoil, 0,0, meshFine))
    pBoundCoil4.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, 0,0, meshFine))
    pBoundCoil4.append(pCoil[len(pCoil)-1])
    pCoil.append(gmsh.model.geo.add_point(origCoil + coilWidth, coreHeight *coilHeightFak ,0, meshFine))
    pCoil.append(gmsh.model.geo.add_point(origCoil,coreHeight *coilHeightFak ,0, meshFine))

    lCoil = []
    lBoundCoil4 = []
    for i in range(0, len(pCoil)):
        lCoil.append(gmsh.model.geo.add_line(pCoil[i], pCoil[(i+1)%len(pCoil)]))
        if i == 0:
            lBoundCoil4.append(lCoil[len(lCoil)-1])

    coilLoop4 = gmsh.model.geo.add_curve_loop(lCoil)

    coil2Neg = gmsh.model.geo.add_plane_surface([coilLoop4])

    # points for boundary
    pBoundArc = []
    pBoundArc.append(gmsh.model.geo.add_point(-boundDist,0,0, meshBoundary))
    pBoundArc.append(gmsh.model.geo.add_point(-boundDist + domainRadius, 0, 0, meshBoundary))
    pBoundArc.append(gmsh.model.geo.add_point(-boundDist+domainRadius*2,0,0, meshBoundary))

    lNeuBdry = []
    lDirBdry = []
    lBound = []
    lBound.append(gmsh.model.geo.add_line(pBoundArc[0],pBoundCoil1[0])) # 0
    lBound.append(lBoundCoil1[0]) # 1
    lBound.append(gmsh.model.geo.add_line(pBoundCoil1[1],pBoundCore[0])) # 2
    lBound.append(lBoundCore[0]) # 3
    lBound.append(gmsh.model.geo.add_line(pBoundCore[1],pBoundCoil2[0])) # 4
    lBound.append(lBoundCoil2[0]) # 5
    lBound.append(gmsh.model.geo.add_line(pBoundCoil2[1],pBoundCore[2])) # 6 
    lBound.append(lBoundCore[1]) # 7
    lBound.append(gmsh.model.geo.add_line(pBoundCore[3],pBoundCoil3[0])) # 8
    lBound.append(lBoundCoil3[0]) # 9
    lBound.append(gmsh.model.geo.add_line(pBoundCoil3[1],pBoundCore[4])) # 10
    lBound.append(lBoundCore[2]) # 11
    lBound.append(gmsh.model.geo.add_line(pBoundCore[5],pBoundCoil4[0])) # 12
    lBound.append(lBoundCoil4[0]) # 13
    lBound.append(gmsh.model.geo.add_line(pBoundCoil4[1],pBoundArc[2])) # 14
    lNeuBdry.extend(lBound)
    lBound.append(gmsh.model.geo.add_circle_arc(pBoundArc[2], pBoundArc[1], pBoundArc[0])) # 15
    lDirBdry.append(lBound[len(lBound)-1])

    
    lBound2 = []
    lBound2.append(lCore[1])
    lBound2.append(lCore[2])
    lBound2.append(lCore[3])
    lBound2.append(-lBound[6])
    lBound2.append(-lBound[5])
    lBound2.append(-lBound[4])

    lBound3 = []
    lBound3.append(lCore[5])
    lBound3.append(lCore[6])
    lBound3.append(lCore[7])
    lBound3.append(-lBound[10])
    lBound3.append(-lBound[9])
    lBound3.append(-lBound[8])


    boundaryLoop = gmsh.model.geo.add_curve_loop(lBound)
    vacSurf = gmsh.model.geo.add_plane_surface([boundaryLoop,-coreLoop, -coilLoop1, -coilLoop4 ])

    vacLoop1 = gmsh.model.geo.add_curve_loop(lBound2,reorient=True) 
    vacLoop2 = gmsh.model.geo.add_curve_loop(lBound3,reorient=True) 
    vacSurf2 =  gmsh.model.geo.add_plane_surface([vacLoop1, -coilLoop2])
    vacSurf3 =  gmsh.model.geo.add_plane_surface([vacLoop2, -coilLoop3])


    gmsh.model.geo.synchronize()
    
    gmsh.model.add_physical_group(dim, [core], tag = sTags['core'], name = "core")
    gmsh.model.add_physical_group(dim, [coil1Neg], tag = sTags['coil1Neg'], name = "coil1Neg")
    gmsh.model.add_physical_group(dim, [coil1Pos], tag = sTags['coil1Pos'], name = "coil1Pos")
    gmsh.model.add_physical_group(dim, [coil2Neg], tag = sTags['coil2Neg'], name = "coil2Neg")
    gmsh.model.add_physical_group(dim, [coil2Pos], tag = sTags['coil2Pos'], name = "coil2Pos")
    gmsh.model.add_physical_group(dim, [vacSurf, vacSurf2, vacSurf3], tag = sTags['vacuum'], name = "vacuum")

    gmsh.model.add_physical_group(dim-1, lDirBdry, tag = lTags['dirBound'], name = "DirichletBoundary")
    gmsh.model.add_physical_group(dim-1, lNeuBdry, tag = lTags['neuBound'], name = "NeumannBoundary")

    gmsh.model.mesh.generate(2)
    gmsh.write(modelName+".msh")

    gmsh.finalize()

    return currScale


## folder setup
scenName = "TeamProblem32_case3"
meshName = "Team32_half"

materialLib = "material_data/xml"
inputDir = scen.create_directory("input")
materialDir = scen.create_directory(inputDir + "/"+ "material")
sourcesDir = scen.create_directory(inputDir + "/"+ "sources")
geomDir = scen.create_directory(inputDir+ "/" +"geometry")


vacTag = 2
coreTag = 1
c1NegTag = 3
c1PosTag = 4
c2NegTag = 5
c2PosTag = 6

dirBoundTag = 2
neuBoundTag = 1

coarseFactor = 2

materialName = "TEAM32_cont_60deg"
timeSeriesXML = "TEAM_Problem_32_case_3_two_periods.xml"

Bsat = 1.412229

## define points for post processing
xMid = 0.08725
evalPoints = []
evalPoints.append([xMid, 0.0615])
evalPoints.append([xMid, 0.0640])
evalPoints.append([xMid, 0.0695])
evalPoints.append([xMid, 0.0720])


## create mesh
meshName = geomDir + "/" + meshName
surfTags = {'vacuum': vacTag, 'core':coreTag, 'coil1Pos':c1PosTag, 'coil1Neg':c1NegTag, 'coil2Pos':c2PosTag, 'coil2Neg':c2NegTag}
lineTags = {'dirBound':dirBoundTag, 'neuBound':neuBoundTag}

baseCurrScale = create_team_problem_32_msh(meshName, surfTags,  lineTags, coarseFactor) # returns currents scale <- dependends on coil area

msh_to_xdmf(meshName)

## copy material parameter file
materialLibFile = materialDir + "/" + materialName + ".xml"
materialFile = scen.search_for_path(materialLibFile) # find absolute path for file


## create time series
timeSeriesFile = sourcesDir + "/" + timeSeriesXML # rel. path sufficient
scen.search_for_path(sourcesDir + "/" + timeSeriesXML) # confirm that file exists


## create scenario xml #######
xmlScen = scen.scenario_xml(scenName)

ironDomains = [coreTag]
coil1Domains = [c1NegTag, c1PosTag]
coil2Domains = [c2NegTag, c2PosTag]
dirBound = dirBoundTag


xmlScen.add_mesh_definition(meshName+".xdmf", meshName+"_facets.xdmf")
xmlScen.add_material_hgm(materialFile, Bsat,ironDomains)
xmlScen.add_current_source_file(timeSeriesFile, "i1", coil1Domains, [-baseCurrScale, baseCurrScale])
xmlScen.add_current_source_file(timeSeriesFile, "i2", coil2Domains, [-baseCurrScale, baseCurrScale])

xmlScen.add_zero_boundary(dirBound)


xmlScen.add_solver_parameters("Newton", 1e-6, 50, 1e-4)
xmlScen.add_time_stepping_file(timeSeriesFile)
xmlScen.add_point_evaluation("B", evalPoints)

xmlScen.write_file(scenName)