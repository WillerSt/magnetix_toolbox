# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

from xml.dom import minidom
import numpy
from math import pi
import os
import shutil

def copy_library_file(libFile, targetDir):
    sourcePath = search_for_path(libFile)
    return shutil.copy(sourcePath, targetDir)

def search_for_path(dir_path:str, dispInfo = True):
    maxRec = 10
    resPath = ""
    foundPath = False
    append = ""
    for i in range(0, maxRec):        
        if os.path.exists(append+dir_path):
            foundPath = True
            resPath = os.path.abspath(append+dir_path)
            if dispInfo == True:
                print("Found file or directory: " + dir_path + " in folder herachy: \n\t" + resPath)

        append = append + "../"


    if foundPath == False:
        print("Directory or file " + dir_path + " not found after in " + str(maxRec) + " parent folders")
    
    return resPath


def create_directory(dir_path:str):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print("Directory "+dir_path+ " created successfully!")
    else:
        print("Directory " +dir_path + " already exists!")
    return dir_path

def define_property(rootNode, name, value):
    param = rootNode.createElement('EntityDefintion') 
    param.setAttribute('name', str(name))
    param.setAttribute('value', str(value))
    return param

def transform_to_list(listIn):
    if not isinstance(listIn, list):
        listIn = [listIn]
    
    return listIn

def array_to_string_list(arrayIn):
    stringList = []
    for i in range (0, len(arrayIn)):
        stringList.append(str(arrayIn[i]))
    return stringList

def string_list_to_string(listIn):
    resString = listIn[0]
    for i in range(1, len(listIn)):
        resString = resString +" "+listIn[i]
    return resString

class scenario_xml:
    def __init__(self, name):
        self.name = str(name)
        self.root = minidom.Document()
        self.xml = self.root.createElement('ScenarioDescription')
        self.xml.setAttribute('name', self.name)
        self.root.appendChild(self.xml)

    def write_file(self, filename):
        xml_str = self.root.toprettyxml(indent ="\t")    
        with open(str(filename)+".xml", "w") as f: 
            f.write(xml_str) 

    def add_mesh_definition(self, mshFile, bdryFile):
        meshDef = self.root.createElement('MeshDefinition') 
        meshDef.appendChild(define_property(self.root, 'meshFile', str(mshFile)))
        meshDef.appendChild(define_property(self.root, 'boundaryFile', str(bdryFile)))
        self.xml.appendChild(meshDef)

    def add_material_hgm(self, distFile, Bsat, idxList):

        for idx in transform_to_list(idxList):
            matDef = self.root.createElement('MaterialDefinition') 
            matDef.appendChild(define_property(self.root, 'type', 'HGM'))
            matDef.appendChild(define_property(self.root, 'distFile', str(distFile)))
            matDef.appendChild(define_property(self.root, 'Bsat', float(Bsat)))
            matDef.appendChild(define_property(self.root, 'domainIndex', int(idx)))
            self.xml.appendChild(matDef)

    def add_material_atan(self, Bsat, muInit, muLin, domainList):
        domainList = transform_to_list(domainList)
        for domain in domainList:
            matDef = self.root.createElement('MaterialDefinition') 
            matDef.appendChild(define_property(self.root, 'type', 'atan'))
            matDef.appendChild(define_property(self.root, 'Bsat', float(Bsat)))
            matDef.appendChild(define_property(self.root, 'muInit', float(muInit)))
            matDef.appendChild(define_property(self.root, 'muLin', float(muLin)))
            matDef.appendChild(define_property(self.root, 'domainIndex', domain))
            self.xml.appendChild(matDef)

    def add_material_linear(self, muR, domainList):
            domainList = transform_to_list(domainList)
            for domain in domainList:
                matDef = self.root.createElement('MaterialDefinition') 
                matDef.appendChild(define_property(self.root, 'type', 'linear'))
                matDef.appendChild(define_property(self.root, 'muR', float(muR)))
                matDef.appendChild(define_property(self.root, 'domainIndex', int(domain)))
                self.xml.appendChild(matDef)

    def add_material_spline(self, xmlFile, domainList):
        domainList = transform_to_list(domainList)
        for domain in domainList:
            matDef = self.root.createElement('MaterialDefinition') 
            matDef.appendChild(define_property(self.root, 'type', 'table'))
            matDef.appendChild(define_property(self.root, 'xmlFile', str(xmlFile)))
            matDef.appendChild(define_property(self.root, 'domainIndex', int(domain)))
            self.xml.appendChild(matDef)

    def add_boundary_definition(self, typeDef, value, idxList):
        for idx in transform_to_list(idxList):                        
            bdryDef = self.root.createElement('BoundaryhDefinition') 
            bdryDef.appendChild(define_property(self.root, 'type', str(typeDef)))
            bdryDef.appendChild(define_property(self.root, 'value', float(value)))
            bdryDef.appendChild(define_property(self.root, 'index', int(idx)))
            self.xml.appendChild(bdryDef)

    def add_zero_boundary(self, idxList):
        self.add_boundary_definition('Dirichlet', 0.0, idxList)     

    def add_solver_parameters(self, solverType, rTol, maxIter, minRes):
        solverParameters = self.root.createElement('SolverParameters')
        solverParameters.appendChild(define_property(self.root, 'type', str(solverType)))
        solverParameters.appendChild(define_property(self.root, 'rTol', float(rTol) ))
        solverParameters.appendChild(define_property(self.root, 'maxIter', int(maxIter)))
        solverParameters.appendChild(define_property(self.root, 'minRes', float(minRes)))
        self.xml.appendChild(solverParameters)

    def add_time_stepping_file(self, fileName):
        timeDef = self.root.createElement('TimeStepping')
        timeDef.appendChild(define_property(self.root, 'type', 'xmlFile'))
        timeDef.appendChild(define_property(self.root, 'xmlFile', str(fileName)))
        self.xml.appendChild(timeDef)

    def add_current_source_file(self, fileName, quantityName, domainList, scaleList):
        domainCount = 0
        domainList = transform_to_list(domainList)
        scaleList = transform_to_list(scaleList)
        srcDef = self.root.createElement('SourceDefinition')
        srcDef.appendChild(define_property(self.root, 'type', 'currentDensity'))
        srcDef.appendChild(define_property(self.root, 'distribution', 'const'))
        srcDef.appendChild(define_property(self.root, 'evolution', 'timeSeries'))
        srcDef.appendChild(define_property(self.root, 'xmlFile', str(fileName)))
        srcDef.appendChild(define_property(self.root, 'quantityName', str(quantityName)))

        for i in range(0, len(domainList)):
            domainCount = domainCount + 1
            srcDef.appendChild(define_property(self.root, 'domain_'+str(domainCount), int(domainList[i])))
            srcDef.appendChild(define_property(self.root, 'scale_'+str(domainCount), scaleList[i]))

        self.xml.appendChild(srcDef)

    def add_current_source_sine(self, amp, freq, phase, offset,domainList, scaleList):
        domainCount = 0
        domainList = transform_to_list(domainList)
        scaleList = transform_to_list(scaleList)

        if (len(scaleList)!= len(scaleList)):
            print("Error: domain lists and scaling lists must have the same length")
        
        srcDef = self.root.createElement('SourceDefinition')
        srcDef.appendChild(define_property(self.root, 'type', 'currentDensity'))
        srcDef.appendChild(define_property(self.root, 'distribution', 'const'))
        srcDef.appendChild(define_property(self.root, 'evolution', 'sine'))

        for i in range(0, len(domainList)):
            domainCount = domainCount + 1
            srcDef.appendChild(define_property(self.root, 'domain_'+str(domainCount), int(domainList[i])))
            srcDef.appendChild(define_property(self.root, 'scale_'+str(domainCount), scaleList[i]))

        srcDef.appendChild(define_property(self.root, 'sineAmp', amp))
        srcDef.appendChild(define_property(self.root, 'sinePhase', phase))
        srcDef.appendChild(define_property(self.root, 'sineFreq', freq))
        srcDef.appendChild(define_property(self.root, 'sineOffset', offset))

        self.xml.appendChild(srcDef)

    def add_current_source_const(self, baseVal, domainList, scaleList):
        self.add_current_source_sine(baseVal, 0.0, 0.5*pi, 0.0, domainList, scaleList)

    def add_linear_pemanent_magnet(self, Br, muR, dirVec, domainList):
        domainList = transform_to_list(domainList)
        for domain in domainList:
            srcDef = self.root.createElement('MaterialDefinition')
            srcDef.appendChild(define_property(self.root, 'type', 'permanentMagnet'))
            srcDef.appendChild(define_property(self.root, 'model', 'linear'))
            srcDef.appendChild(define_property(self.root, 'Br', float(Br)))
            srcDef.appendChild(define_property(self.root, 'muR', float(muR)))
            srcDef.appendChild(define_property(self.root, 'xDir', float(dirVec[0])))
            srcDef.appendChild(define_property(self.root, 'yDir', float(dirVec[1])))
            srcDef.appendChild(define_property(self.root, 'domainIndex', int(domain)))
            self.xml.appendChild(srcDef)

            self.add_material_linear(muR, domain)

    def add_force_calculation(self, domainList, domDepth, rotCenter):
        domainCount = 0
        domainList = transform_to_list(domainList)
        forceDef = self.root.createElement('ForceCalculation')
        forceDef.appendChild(define_property(self.root, 'type', 'ForceCalculation'))
        forceDef.appendChild(define_property(self.root, 'domainDepth', str(domDepth)))
        forceDef.appendChild(define_property(self.root, 'rotCenter_x', str(rotCenter[0])))
        forceDef.appendChild(define_property(self.root, 'rotCenter_y', str(rotCenter[1])))

        for i in range(0, len(domainList)):
            domainCount = domainCount + 1
            forceDef.appendChild(define_property(self.root, 'domain_'+str(domainCount), int(domainList[i])))

        self.xml.appendChild(forceDef)


    def add_point_evaluation(self, quantity, queryPoints):
        for p in queryPoints:
            pDef = self.root.createElement('PointEvaluation')
            pDef.appendChild(define_property(self.root, 'quantity', str(quantity)))
            pDef.appendChild(define_property(self.root, 'x', str(p[0])))
            pDef.appendChild(define_property(self.root, 'y', str(p[1])))
            if len(p) > 2:
                pDef.appendChild(define_property(self.root, 'y', str(p[2])))
            self.xml.appendChild(pDef)


    def add_domain_rotation(self, domainList, interfaceList, rotCenter, timeXML ):
        domainCount = 0
        domainList = transform_to_list(domainList)
        interfaceList = transform_to_list(interfaceList)
        if len(interfaceList)>2:
            print('Warning interface list with length ' + str(len(interfaceList)))

        rotDef = self.root.createElement('MovingDomain')
        rotDef.appendChild(define_property(self.root, 'type', 'DomainRotation'))
        rotDef.appendChild(define_property(self.root, 'rotCenter_x', str(rotCenter[0])))
        rotDef.appendChild(define_property(self.root, 'rotCenter_y', str(rotCenter[1])))
        rotDef.appendChild(define_property(self.root, 'masterInterface', str(interfaceList[0])))
        rotDef.appendChild(define_property(self.root, 'slaveInterface', str(interfaceList[1])))
        rotDef.appendChild(define_property(self.root, 'evolution', 'timeSeries'))
        rotDef.appendChild(define_property(self.root, 'xmlFile', str(timeXML)))

        for i in range(0, len(domainList)):
            domainCount = domainCount + 1
            rotDef.appendChild(define_property(self.root, 'domain_'+str(domainCount), int(domainList[i])))

        self.xml.appendChild(rotDef)

class result_xml:
    def __init__(self, fileName):
        self.fileName = fileName
        self.dom = minidom.parse(str(fileName))
        self.elements = self.dom.getElementsByTagName("TimeVariation")

    def gather_point_queries(self, quantity):
        pQuery = []
        for series in self.elements:
            
            if str(series.attributes['name'].value) == 'pointQuery':
                if str(series.attributes['quantity'].value) == str(quantity):
                    noElem = 1
                    if str(series.attributes['format'].value) == 'xy':
                        noElem = 2
                    elif str(series.attributes['format'].value) == 'xy':
                        noElem = 3
                    vals = str.split(series.firstChild.nodeValue)
                    f = []
                    for i in range(0, int(len(vals)/noElem)):
                        tempVal = []
                        for j in range(0, noElem):
                            tempVal.append(float(vals[noElem*i+j]))
                        
                        f.append(numpy.asarray(tempVal))
                    
                    pQuery.append(f)
        return pQuery

    


    def get_force(self, domainIdx):
        breakIt = False
        for series in self.elements:
            if str(series.attributes['name'].value) == 'F':
                if int(series.attributes['domain'].value) == domainIdx:
                    vals = str.split(series.firstChild.nodeValue)
                    F =[]
                    Fx = []
                    Fy = []
                    for i in range(0, int(len(vals)/2)):
                        F.append(numpy.asarray([float(vals[2*i]), float(vals[2*i+1])]))
                        Fx.append(float(vals[2*i]))
                        Fy.append(float(vals[2*i+1]))
                        
                    #F = numpy.asarray([1,1])
            if breakIt == True:
                break
        return F, Fx, Fy
    
    def get_torque(self, domainIdx):
        breakIt = False
        for series in self.elements:
            if str(series.attributes['name'].value) == 'Tz':
                if int(series.attributes['domain'].value) == domainIdx:
                    vals = str.split(series.firstChild.nodeValue)
                    T  = []

                    for i in range(0, len(vals)):
                        T.append(float(vals[i]))

            if breakIt == True:
                break

        T = numpy.asarray(T)
        return T
    
    def get_scalar_quantity(self, name):
        breakIt = False
        for series in self.elements:
            if str(series.attributes['name'].value) == name:
                vals = str.split(series.firstChild.nodeValue)
                T  = []

                for i in range(0, len(vals)):
                    T.append(float(vals[i]))
                        
            if breakIt == True:
                break

        T = numpy.asarray(T)
        return T

class material_xml:
    def __init__(self, name):
        self.name = str(name)
        self.root = minidom.Document()
        self.xml = self.root.createElement('MaterialDefinition')
        self.xml.setAttribute('name', self.name)
        self.root.appendChild(self.xml)

    def define_series(self, name, value, **kwargs):
    

        param = self.root.createElement('MaterialData') 
        param.setAttribute('quantity', str(name))

        if 'unit' in kwargs:
            param.setAttribute('unit', kwargs['unit'])

        textNode =  self.root.createTextNode(string_list_to_string(array_to_string_list(value)))
        param.appendChild(textNode)
        self.xml.appendChild(param)

    def write_file(self, filename):
        xml_str = self.root.toprettyxml(indent ="\t")    
        with open(str(filename)+".xml", "w") as f: 
            f.write(xml_str) 

class time_series_xml:
    def __init__(self, name):
        self.name = str(name)
        self.root = minidom.Document()
        self.xml = self.root.createElement('TimeSeries')
        self.xml.setAttribute('name', self.name)
        self.root.appendChild(self.xml)

    def define_series(self, name, value, **kwargs):
    

        param = self.root.createElement('TimeVariation') 
        param.setAttribute('quantity', str(name))

        if 'unit' in kwargs:
            param.setAttribute('unit', kwargs['unit'])

        textNode =  self.root.createTextNode(string_list_to_string(array_to_string_list(value)))
        param.appendChild(textNode)
        self.xml.appendChild(param)

    def write_file(self, filename):
        xml_str = self.root.toprettyxml(indent ="\t")    
        with open(str(filename)+".xml", "w") as f: 
            f.write(xml_str) 

def split_components(f):
    noElem = numpy.ndim(f[0])+1
    splitRes = []
    for j in range(0, noElem):
        fComp = []
        for i in range(0,len(f)):
            fComp.append(f[i][j])
        splitRes.append(numpy.asarray(fComp))
    return splitRes
            
