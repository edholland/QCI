#!/usr/bin/python
from ConfigParser import SafeConfigParser
from string import Template
import logging
import os
import re
import math
import subprocess
import cclib
import pybel
import openbabel
import numpy

class project(object):
    """Superclass for projects

    Currently holds projectName and makes folder if required
    Attributes:
        projectName: A string holding the project name
    """

    def __init__(self,name):
        """Makes project superclass
        Makes folder is required"""
        self.projectName = name
        try:
            os.makedirs(os.path.join('jobs',self.projectName))
            self.exists = True
        except OSError:
            self.exists = True 
            
class job(project):
    """Superclass for jobs inherits from project superclass

    Currently holds jobName and makes folder if required
    Attributes:
    jobsName: A string holding the job name
    """
    
    def __init__(self,name,projectName):
        """Job Class: subclass of project
        makes job folder if required
        """
        super(job, self).__init__(projectName)
        self.jobName = name
        try:
            os.makedirs(os.path.join('jobs',self.projectName,self.jobName))
            self.exists = True
        except OSError:
            self.exists = True

class step(job):
    """Superclass for step inherits from job superclass
    Holds methods independant of engine choice and calculation
    paramaters. Methods include analyse_job (calls cclib to extract data from
    outputfiles), remove_pbsfile (deletes old pbs file), submit_job (calls qsub
    with pbf file), make_xyz (returns a list of coords given xyz file),
    find_h_spheres (finds H atoms not within 1.2 angstroms from a heavy atom
    for pcm)
    Attributes:
    stepName: A string for holding step name
    ts: A booleanfor holding transition state data
    charge: An int holding charge for calculation
    mult: An int holding mutiplicity for calculation
    xyz: A file handler for an xyz file in read mode
    coords: A list of lists holding coordinates
    h_spheres: A list of ints holding atom numbers for H atoms
               needing H spheres in PCM model
    energy: A float holding SCF energy
    enthalpy: A float holding enthalpy
    freeenergy: A float holding gibbs free energy
    entropy: A float holding entropy
    """

    def __init__(self, projectName, jobName, stepName, **kwargs):
        """Init method for step superclass, many attributes are made here.

        Makes attributes, calling engine specific methods where needed. There
        should definitly be some checking here that subclasses correctly
        implement the methods that will be used here! Also makes step folder if
        needed

        Atrributes:
            projectName: String for project name
                         (passed up to project superclass)
            jobName: String for job name (passed up to job superclass)
            stepName: A string for holding step name
            **ts: A boolean for holding transition state data
            **charge: An int holding charge for calculation
            **mult: An int holding mutiplicity for calculation
            **xyz: A file handler for an xyz file in read mode
            **fragatom: An int holding atom number to be subsituted
            **frag: A string holding fragment name
        """
        super(step, self).__init__(jobName, projectName)
        
        self.stepName = stepName
        
        logging.debug("starting init:- %s-%s-%s"
                      % (self.projectName,self.jobName,self.stepName))

        if 'ts' in kwargs:
            self.ts = kwargs.get('ts')
            logging.debug('ts is %s' % self.ts)
        if 'charge' in kwargs:
            self.charge = int(kwargs.get('charge'))
            logging.debug('charge is %s' % self.charge)
        if 'mult' in kwargs:
            self.mult = int(kwargs.get('mult'))
            logging.debug('mult is %s' % self.mult)
        if 'xyz' in kwargs:
            self.mol = self.make_xyz(kwargs.get('xyz'))
            if 'fragatom' in kwargs and 'frag' in kwargs:
                self.mol = self.replace_atom_by_frag(self.mol, kwargs.get('fragatom'), kwargs.get('frag'))
            self.h_spheres = self.find_h_spheres()


            logging.debug('h_spheres is %s' % self.h_spheres)

        try:
            os.makedirs(os.path.join('jobs', self.projectName,
                                     self.jobName, self.stepName, 'checkpoint'))
            logging.debug('Step does not exisit, creating folder')
            self.exists = True
        except OSError:
            self.exists = True

            logging.debug('finished init: - %s-%s-%s'
                          % (self.projectName,self.jobName,self.stepName))
            
    def __del__(self):
        """Log job name when removing object"""
        logging.debug("Removing step object %s-%s-%s"
                      % (self.projectName, self.jobName, self.stepName))

    def __str__(self):
        """Return basic info as string representation of step"""
        lines = []
        lines.append("Job: %s-%s-%s"
                     % (self.projectName,self.jobName,self.stepName))
        if hasattr(self, 'converged'):
            lines.append('CONVERGED')
        if hasattr(self, 'mol'):
            geom = self.format_coords(self.mol)
            lines.append("Coordinates:")
            lines.append(geom)
        if hasattr(self, 'freeenergy'):
            lines.append("Free energy = %f kcal/mol" % self.freeenergy)
        if hasattr(self, 'enthalpy'):
            lines.append('Enthalpy = %f kcal/mol' % self.enthalpy)
        if hasattr(self, 'entropy'):
            lines.append('Entropy = %f cal/mol-K' % (self.entropy*1000,))
        if hasattr(self, 'temperature'):
            lines.append('Calculated at %f K' % self.temperature)
        lines.append(' ')
        lines.append(' ')
        return "\n".join(lines)
    
    def replace_atom_by_frag(self, mol, atomno, frag):
        """Replaces a given atom with a given fragment"""

        t = cclib.parser.utils.PeriodicTable()        
        atom = [x for x in mol if x.idx == int(atomno)]
        neighbour = [x for x in openbabel.OBAtomAtomIter(atom[0].OBAtom)][0]
        bond = [x for x in openbabel.OBAtomBondIter(atom[0].OBAtom)][0]

        if frag == "CH3":
            atom[0].OBAtom.HtoMethyl()

        elif frag in ('F', 'Cl', 'Br', 'I'):
            atom[0].OBAtom.SetAtomicNum(t.number[frag])
            bond.SetLength(neighbour, bond.GetEquibLength())

        else:
            mol.OBMol.DeleteAtom(atom[0].OBAtom)
            with open('conf/frags/%s.def' % frag, 'r') as def_file:
                new_atoms = []
                count = 0
                for line in def_file:
                    split = line.split()
                    new_atoms.append(mol.OBMol.NewAtom())
                    new_atoms[-1].SetAtomicNum(t.number[split[0]])
                    new_atoms[-1].SetVector(float(split[1]), float(split[2]), float(split[3]))
                    if count > 0:
                        mol.OBMol.AddBond(new_atoms[-2].GetIdx(), new_atoms[-1].GetIdx(), 1)
                    count += 1
            openbabel.OBBuilder.Connect(mol.OBMol, neighbour.GetIdx(),new_atoms[0].GetIdx())  

        
        return mol

        
    def analyse_job(self):
        """Method for extracting useful data from output files

        Uses an external module, cclib, for the parsing of outfile. More info
        about cclib can be found at
        http://sourceforge.net/apps/mediawiki/cclib/index.php?title=Main_Page
        Currently modded version of cclib is required to extract thermochemistry
        data, however hopefully this is be included in the trunk soon. This
        method will extract optimised coordinates from finished jobs, and final
        coordinates from jobs that do not seem to have finished correctly, they
        are stored in coords (an engine specific formatted string). Also
        extracts thermochemistry data from freq calculations for use in
        calculating reaction paramaters.
        Attributes:
            coords: A string of engine specific formatted coordinates
            energy: A float holding SCF energy
            enthalpy: A float holding enthalpy
            freeenergy: A float holding gibbs free energy
            entropy: A float holding entropy
        Returns:
            True if no failues occur
        Throws:
            IOError: Cant open outfile with cclib
        """

        
        try:
            outfile = cclib.parser.ccopen(os.path.join(os.getcwd(), 'jobs',
                                              self.projectName, self.jobName,
                                              self.stepName, '%s-%s-%s.out'
                                              % (self.projectName,
                                                 self.jobName,
                                                 self.stepName)))
            outfile.logger.setLevel(logging.ERROR)
            data=outfile.parse()
        except (IOError, AttributeError):
            e = 'Job: %s %s %s' % (self.projectName, self.jobName, self.stepName)
            raise IOError, e

        final = None
        self.converged = False
        if (hasattr(data,'geovalues') and hasattr(data,'geotargets')
            and hasattr(data,'atomnos')):
            final = 0
            count = 0
            t = cclib.parser.utils.PeriodicTable()
            #look for converged step (final)
            for step in data.geovalues:
                if all(a < b for (a,b) in zip(step,data.geotargets)):
                    self.converged = True
                    final=count
                count += 1
            #gets converged geom or final (should be same for complete jobs)
            if self.converged is True:
                self.mol = pybel.Molecule(cclib.bridge.makeopenbabel(data.atomcoords[final-1], data.atomnos, charge = data.charge, mult = data.mult))
            else:
                self.mol = pybel.Molecule(cclib.bridge.makeopenbabel(data.atomcoords[-1], data.atomnos, charge = data.charge, mult = data.mult))
            
        if hasattr(data,'scfenergies') and self.converged is True:
            self.energy = data.scfenergies[final]
            logging.debug("extracted converged scf: -%s-%s-%s"
                              % (self.projectName,self.jobName,self.stepName))
        elif hasattr(data,'scfenergies'):
            self.energy = data.scfenergies[-1]
            logging.debug("extracted finial scf: -%s-%s-%s"
                          % (self.projectName,self.jobName,self.stepName))
 
        #gets thermochemistry
        if hasattr(data, 'enthalpy'):
            self.enthalpy = data.enthalpy
            logging.debug("extracted enthalpy: -%s-%s-%s"
                              % (self.projectName,self.jobName,self.stepName))
        if hasattr(data, 'freeenergy'):
            self.freeenergy = data.freeenergy
            logging.debug("extracted free energy: -%s-%s-%s"
                              % (self.projectName,self.jobName,self.stepName))
        if hasattr(data, 'entropy'):
            self.entropy = data.entropy
            logging.debug("extracted entropy: -%s-%s-%s"
                          % (self.projectName,self.jobName,self.stepName))
        if hasattr(data, 'temperature'):
            self.temperature = data.temperature
            logging.debug("extracted temp: -%s-%s-%s"
                          % (self.projectName,self.jobName,self.stepName))


        if hasattr(data, 'irccoords'):
            self.irccoords = []
            for geom in data.irccoords:
                self.irccoords.append(pybel.Molecule(cclib.bridge.makeopenbabel(geom, data.atomnos, charge = data.charge, mult = data.mult)))
        if hasattr(data, 'ircenergies'):
            self.ircenergies = data.ircenergies
        if hasattr(data, 'rxcoord'):
            self.rxcoord = data.rxcoord

        
        return True

    def submit_job(self):
        """Submit job using pbs file

        calls qsub to submit job to PBS.Should return success boolean and throw
        revevant exceptions"""

        pbsfile = '%s-%s-%s.pbs' % (
            self.projectName, self.jobName, self.stepName)
        wdir = os.path.join(os.getcwd(), 'jobs', self.projectName, self.jobName, self.stepName)
        subprocess.call(['qsub',pbsfile] , cwd=wdir)
        logging.debug("pbfile: %s \n wdir: %s" % (pbsfile, wdir))
                
    def make_xyz(self, cfile):
        """return a pybel molecule from a filename string

        Takes a string (passed as a command line argument) holding a filename. This file
        has to be either an .xyz or .mdl file (extension matters). The file is parsed by
        openbabel to extract the first molecule in the file.
        Can also take ".job" is a vurtial file that references previous jobs. The respective output
        files are parsed for geometries and were possible an optimised geometery is used, otherwise
        the final geometry is used
        Returns:
            firstmol: A pybel molecule
        Throws:
            IOError: Cant read given file"""
        if ".xyz" in cfile:
            try:
                firstmol = pybel.readfile('xyz',cfile).next()
                return firstmol
            except IOError:
                error = "Cant open .xyz file"
                raise IOError, error
        elif '.mdl' in cfile:
            try:
                firstmol = pybel.readfile('mdl',cfile).next()
                return firstmol
            except IOError:
                error = "Cant open .mdl file"
                raise IOError, error
        elif '.job' in cfile:
            names = cfile[:-4].split('-')
            j = step(names[0],names[1], names[2])
            try:
                j.analyse_job()
                firstmol = j.mol
                return firstmol
            except IOError:
                error = "Cant open .job output file for analysis"
                raise IOError, error
        else:
            error = "File type not supported currently"
            raise IOError, error

    def make_distance_matrix(self, mol):
        """Takes pybel molecule and makes a complete distance matrix"""
        i = 0
        j = 0
        output=numpy.empty((len(mol.atoms),len(mol.atoms)))
        for atomi in mol:
            for atomj in mol:
                distance =  math.sqrt( math.pow(atomi.coords[0]-atomj.coords[0],2) + 
                                       math.pow(atomi.coords[1]-atomj.coords[1],2) + 
                                       math.pow(atomi.coords[2]-atomj.coords[2],2))
                output[i,j] = distance
                output[j,i] = distance
                j += 1
            i += 1
            j = 0
        return output
            

    def find_h_spheres(self):
        """Make list of H atoms needing pcm spheres

        Takes an openbabel molecule and finds H atoms that are over
        the threshold define in the conf file away from a heavy atom
        This is used to produce a list of H atoms that will require PCM
        spheres added in g03"""
        
        spheres = []
        #Loops over H atoms
        for Hatom in self.mol:
            lowest = 1000
            if Hatom.atomicnum == 1:
                for heavyatom in self.mol:
                    if heavyatom.atomicnum != 1:
                        distance =  math.sqrt( math.pow(Hatom.coords[0]-heavyatom.coords[0],2) + 
                                               math.pow(Hatom.coords[1]-heavyatom.coords[1],2) + 
                                               math.pow(Hatom.coords[2]-heavyatom.coords[2],2))
                        lowest = min(lowest, distance)
                if lowest > 1.15:
                    #adds H atoms further than threshold to list for H spheres
                    spheres.append(Hatom.idx)
                        
        return spheres

class gau_step(step):
    """Class for Gaussian steps: subclass of step.

    Specific details on how to write gaussian inputfiles and pbs files

    Methods implemented:
    write_inputfile (generates input file based on definition file)
    write_pbsfile (generates pbs file based on definition file)
    format_coords (formats a pybel molecule to a gaussian formatted string)
    """

    def __init__(self,projectName, jobName, stepName, **kwargs):
        """Makes gaussian step class
        most initation done in superclass"""
        super(gau_step, self).__init__(projectName, jobName, stepName, **kwargs)
        logging.debug("Making gaussian-step object %s" % self.jobName)
        logging.debug("Full gussaian step name %s-%s-%s"
                     % (self.projectName,self.jobName,self.stepName))

    def write_inputfile(self, basis, functional, nodes, cpus, ram, pcm, type):
        """Writes gaussian specific input files.

        Checks for exisitance of input file and fails if exists otherwise writes
        a gaussian specific inputfile to the relevant location based on the
        definition file and args passed
        Args:
            basis: Set basis set for calculation. Default from config file
            functional: Set functional for calculation. Default from config file
            nodes: Set number of nodes to be used. Default from config file
            cpus: Set number of cpus per node. Default from config file
            ram: Set amount of ram per node. Default from config file
            pcm: Set PCM simulation on and select solvent. Deafult = false 
            type: Set type of job, currently allowed OPT and FREQ. Default = OPT
        """
               
        # handles pcm models and any spheres that need adding on H atoms
        if pcm is None:
            pcm_string = ''
            sphere_string = ''
        else:
            pcm_string = 'SCRF=(PCM, SOLVENT=%s)' % pcm
            sphere_string = ''
            if len(self.h_spheres) > 0:
                pcm_string = 'SCRF=(PCM, SOLVENT=%s, READ)' % pcm
                sphere_string ='sphereonh = '
                for atoms in self.h_spheres:
                    sphere_string = sphere_string + "%d ," % atoms
                sphere_string = sphere_string[:-2]
                    
            
        #Implement more job types here
        if type == 'OPT':
            type_string = 'OPT FREQ'
            if self.ts is True:
                type_string = 'OPT=(CalcFC, TS, NOEIGEN) FREQ'
        if type == 'FREQ':
            type_string = 'FREQ'
        if type == 'IRC':
            type_string = 'IRC=(CalcFC, MaxPoints=10)'
        
            
        #open job and definition file
        with open(os.path.join('jobs', self.projectName, self.jobName, self.stepName,
                                      '%s-%s-%s.inp' % (self.projectName, self.jobName, 
                                                        self.stepName)), "w") as inputfile:
            with open('conf/gaussian_input.def', 'r') as def_file:
                content = def_file.read()
            
                #This is where template definition is held
                d = dict(projectName = self.projectName, jobName = self.jobName,
                         stepName = self.stepName, nodes = nodes,  cpus = cpus,
                         ram = ram, functional = functional, basis = basis,
                         PCM = pcm_string, type = type_string, charge = self.charge,
                         mult = self.mult, xyz = self.format_coords(self.mol),
                         sphere_string = sphere_string)
                inputfile.write(Template(content).safe_substitute(d))
                
    def write_pbsfile(self, walltime, nodes, cpus, queue):        
        """Writes gaussian specific pbs files.

        Checks for exisitance of pbs file and fails if exists otherwise writes
        a gaussian specific pbs file to the relevant location based on the
        definition file and args passed
        Args:
            walltime: Set max allowed excution time. Default from config file
            nodes: Set number of nodes to be used. Default from config file
            cpus: Set number of cpus per node. Default from config file
            queue/-q -- Set queue for job. Default from config file
        """        
        
        with open(os.path.join('jobs', self.projectName, self.jobName, self.stepName,
                    '%s-%s-%s.pbs' % ( self.projectName, self.jobName, self.stepName)),
                  "w") as pbsfile: 
            with open('conf/gaussian_pbs.def', 'r') as def_file:
                content = def_file.read()
                
                #template file definition here
                d = dict(projectName = self.projectName, jobName = self.jobName,
                         stepName = self.stepName, nodes = nodes, cpus = cpus,
                         walltime = walltime, queue = queue)
                pbsfile.write(Template(content).safe_substitute(d))

    def format_coords(self, mol):
        """Produces a Gaussian specific coordinate string

        Takes a openbabel molecule and produces a string specific to 
        gaussian input files"""

        formatted = []
        t = cclib.parser.utils.PeriodicTable()
        for atom in mol:
            formatted.append("%s \t %10.5f \t %10.5f \t %10.5f" % (t.element[atom.atomicnum], atom.coords[0], atom.coords[1], atom.coords[2]))
        return "\n".join(formatted)
            

    
class gamess_step(step):
    """Class for Gamess steps: subclass of step.

    Specific details on how to write gaussian inputfiles and pbs files

    Methods implemented:
    write_inputfile (generates input file based on definition file)
    write_pbsfile (generates pbs file based on definition file)
    format_coords (formats a pybel molecule to a gaussian formatted string)
    """

    def __init__(self,projectName, jobName, stepName, **kwargs):
        """Makes Gamess step class
        most initation done in superclass"""
        super(gamess_step, self).__init__(projectName, jobName, stepName, **kwargs)
        logging.debug("Making Gamess-step object %s" % self.jobName)
        logging.debug("Full Gamess step name %s-%s-%s"
                     % (self.projectName,self.jobName,self.stepName))

    def write_inputfile(self, basis, functional, nodes, cpus, ram, pcm, type, sym):
        """Writes Gamess specific input files.

        Checks for exisitance of input file and fails if exists otherwise writes
        a Gamess specific inputfile to the relevant location based on the
        definition file and args passed
        Args:
            basis: Set basis set for calculation. Default from config file
            functional: Set functional for calculation. Default from config file
            nodes: Set number of nodes to be used. Default from config file
            cpus: Set number of cpus per node. Default from config file
            ram: Set amount of ram per node. Default from config file
            pcm: Set PCM simulation on and select solvent. Deafult = false 
            type: Set type of job, currently allowed OPT and FREQ. Default = OPT
        """
        
        extras = []
                                
        #Implement more job types here
        if type == 'OPT':
            type_string = 'OPTIMIZE'
            if self.ts is True:
                type_string = 'SADPOINT'
                extras.append('$STATPT HESS=CALC HSSEND=.TRUE. $END')
            else:
                extras.append('$STATPT HSSEND=.TRUE. $END')
        if type == 'FREQ':
            type_string = 'HESSIAN'
        if type == 'IRC':
            type_string = 'IRC=(CalcFC, MaxPoints=10)'
        
        #BASIS sets
        if basis == '6-31G(d,p)':
            basis = "GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1"
        if basis == '6-31G':
            basis = "GBASIS=N31 NGAUSS=6"
        if basis == '6-31+G(d,p)':
            basis = "GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 DIFFSP=TRUE"
        if basis == 'STO-3G':
            basis = "GBASIS=STO NGAUSS=3"
        
        #SCF type
        if int(self.mult) > 1:
            scf = "ROHF"
        else:
            scf = "RHF"
        
        #Turn ram into mwords
        if ram[-2:] == 'mb':
            ram = int(ram[:-2])/8
        elif ram[-2:] == 'gb':
            ram = (1000*int(ram[:-2]))/8
        
            
        if sym.upper() == "DNH 4":
            sym += '\n'
            if len(self.mol.atoms) == 2:
                lastatom = [x for x in self.mol][-1]
                self.mol.OBMol.DeleteAtom(lastatom.OBAtom)
                

        #open job and definition file
        with open(os.path.join('jobs', self.projectName, self.jobName, self.stepName,
                                      '%s-%s-%s.inp' % (self.projectName, self.jobName, 
                                                        self.stepName)), "w") as inputfile:
            with open('conf/gamess_input.def', 'r') as def_file:
                content = def_file.read()
            
                #This is where template definition is held
                d = dict(projectName = self.projectName, jobName = self.jobName,
                         stepName = self.stepName, ram = ram, functional = functional,
                         basis = basis, scf = scf, type = type_string, 
                         charge = self.charge, extras = "/n".join(extras),
                         mult = self.mult, xyz = self.format_coords(self.mol), sym = sym)
                inputfile.write(Template(content).safe_substitute(d))
                
    def write_pbsfile(self, walltime, nodes, cpus, queue):        
        """Writes Gamess specific pbs files.

        Checks for exisitance of pbs file and fails if exists otherwise writes
        a gaussian specific pbs file to the relevant location based on the
        definition file and args passed
        Args:
            walltime: Set max allowed excution time. Default from config file
            nodes: Set number of nodes to be used. Default from config file
            cpus: Set number of cpus per node. Default from config file
            queue/-q -- Set queue for job. Default from config file
        """        
        
        with open(os.path.join('jobs', self.projectName, self.jobName, self.stepName,
                    '%s-%s-%s.pbs' % ( self.projectName, self.jobName, self.stepName)),
                  "w") as pbsfile: 
            with open('conf/gamess_pbs.def', 'r') as def_file:
                content = def_file.read()
                
                #template file definition here
                d = dict(projectName = self.projectName, jobName = self.jobName,
                         stepName = self.stepName, nodes = nodes, cpus = cpus,
                         walltime = walltime, queue = queue)
                pbsfile.write(Template(content).safe_substitute(d))

    def format_coords(self, mol):
        """Produces a Gamess specific coordinate string

        Takes a openbabel molecule and produces a string specific to 
        gaussian input files"""

        formatted = []
        t = cclib.parser.utils.PeriodicTable()
        for atom in mol:
            formatted.append("%s \t %10.1f \t %10.5f \t %10.5f \t %10.5f" % (t.element[atom.atomicnum], atom.atomicnum, atom.coords[0], atom.coords[1], atom.coords[2]))
        return "\n".join(formatted)
            

    
