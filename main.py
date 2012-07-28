#!/usr/bin/python 
import logging
import classes
import argparse
from ConfigParser import SafeConfigParser

def main():
    """Quantum Chemistry Interface - Fesigned to simplify the use of QC packages

    Quantum Chemistry Interface has been designed to simplify quantum chemistry
    calculations for organic chemists, it provides a common interface to quantum
    chemistry packages, automatically handles creating input file and reading
    output files, furthermore it is capable of automatically calculating
    reaction thermodynamic quantities given basic reaction path info. QCI
    understands the PBS system so is ideal for use on many clusters

    When the script is run you must call it with the 3 positional arguments
    described below, these tell the script how to name your job, and where the
    files relating to each job are found. Both the job and step names can accept
    a list of names to automatically create batch jobs. Warning: Batch modes can
    radpidly create a large number of job, be careful if there is limits to job
    submission on your cluster.
    Furthermore you must provide a flag indicating which mode the script should
    run in. Currently the options are new (creates a new job, no batch mode),
    submit (submits jobs to the queue), analyse (extracts data from jobs),
    reaction, (extracts thermochemistry for reactions)

    It is invisioned that the project name should hold a braod description of
    a series of calulations, the job name should hold a description of the
    substitution pattern, and the step name should hold a description of the
    reaction step.
    For example in a series of calculations studying H atom loss from
    substituted methane moecules you might start with the starting material,
    methane.

    projectName: - Hloss
    jobName: - H
    stepname: - SM

    qci Hloss H SM -n -c 0 -m 1 -x methane.xyz

    The values not specifically passed are taken from the config file
    Next you would calculate the transition state for loss of H atom

    projectName: - Hloss
    jobName: - H
    stepName -HlossTS

    qci Hloss H HlossTS -n -ts -c 0 -m 1 -x HlossTS.xyz

    Notice we have to explictly tell qci we are looking for a TS state
    Next we need to make the product job for methyl radical
    projectName: - Hloss
    jobName: - H
    stepName: - Product

    qci Hloss H Product -n -c 0 -m 2 -x methylradical.xyz

    We have now created the input files for this reaction, we can submit to the
    queue using a batch mode by passing a list of stepnames

    qci Hloss H "SM HlossTS Product" -s -w 4:00:00

    This will submit the jobs we created in the previous commands.
    Once this jobs are complete we can analyse the job to check the geometery
    has converged to the correct structure, again we can use a batch mode here

    qci Hloss H "SM HlossTS Product" -a

    This will print the most important infomaition from the output file
    If all the jobs have finished correctly, and freqency calculations
    were requested it is now possible to calculate reaction paramters.

    qci Hloss H "SM HlossTS Product+H" -re

    Reaction mode requires a list of reaction steps, it is assumed that odd
    numbered steps are starting materials and products, and that even steps are
    transition states connecting them. If a step has more than one molecule you
    ask qci to sum their enthalpies by giving a list seperated by + signs, with
    no spaces. 
    
    Positional Args:
        projectName -- Set name of project
        jobName -- Set name of job. Accepts a list of names for batch modes
        stepName -- Set name of Step. Accepts a list of name for batch modes

    Optional Args:
        Modes:
            new/-n -- Creates inputfile files from options selected
            submit/-s -- Creates pbsfile and submits job to queue. Batch mode
            analyse/-a -- Reads output file with cclib to extract chemical data.
                          Batch mode
            reaction/-re -- Calculates therodynamics quantites for reaction
                            given list in form "SM TS Prod TS2 Prod2 ... TSN
                            ProdN" where each word is a step name
        Chemical Options:
            tstate/-ts -- Mark job as a transition state. Default = False 
            charge/-c -- Set charge for molecule. Default from config file 
            mult/-m -- Set multiplicity for moecule. Default from config file
            basis/-b -- Set basis set for calculation. Default from config file
            functional/-f -- Set functional for system. Default from config file
            xyz/-x -- Set location to read .xyz file for coordinates. No default
            pcm/-p -- Set PCM simulation on and select solvent. Deafult = false 
            type/-t -- Set type of job, currently either OPT and FREQ.
                       Default = OPT
            engine/-e -- Set which QC packake to use for calculation.
                         Default from config file
        PBS Options:
            walltime/-w -- Set allowed excution time. Default from config file
            nodes/-n -- Set number of nodes to be used. Default from config file
            cpus/-cp -- Set number of cpus per node. Default from config file
            ram/-r -- Set amount of ram per node. Default from config file
            queue/-q -- Set queue for job. Default from config file
    """

    logging.basicConfig(filename='qcfoc.log', level=logging.DEBUG)
    logging.info('Started')
    
    confparser = SafeConfigParser()
    confparser.read('conf/qcfoc.conf')

    argparser = argparse.ArgumentParser(
        description='Quantum Chemistry Interface', version='0.1')
    argparser.add_argument('projectName', action='store',
                           help='Set name of project')
    argparser.add_argument('jobName', action='store',
                           help= ('Set name of job. Accepts a list '
                           'of names for batch modes'))
    argparser.add_argument('stepName', action='store', help=
                           ('Set name of Step. Accepts a list'
                           'of name for batch modes'))

    argparser.add_argument('--new', '-n', action='store_true', default=False,
                           help='Creates inputfile files from options selecte')
    argparser.add_argument('--submit', '-s', action='store_true', default=False,
                           help= ('Creates pbsfile and submits'
                                  'job to queue. Batch mode'))
    argparser.add_argument('--analyse', '-a', action='store_true',
                           default=False, help= ('Creates pbsfile and submits '
                                                 'job to queue. Batch mode'))
    argparser.add_argument('--reaction', '-re', action='store_true', default=False,
                           help=('Calcuates reaction thermochem given reaction'
                                  'profile in stepname in form of'
                                  'SM TS Prod1+Prod2'))
    argparser.add_argument('--irc', '-irc', action='store_true', default=False,
                           help='Checks if IRC jobs give correct SM + Prod')

    argparser.add_argument('--fragment', '-fr', action='store',
                           default=None, help= ('Automatically subtitute'
                                                'atom with fragment. pass'
                                                'as "AtomNo Frag"'))
    argparser.add_argument('--tstate', '-ts', action='store_true',
                           default=False, help= ('Flag to mark job as a '
                           'transition state. Default = False'))
    argparser.add_argument('--charge', '-c', action='store',
                           default=confparser.get('default_options','charge'),
                           help= ('Set charge for molecule. Default from config'
                                  'file'), type=int)
    argparser.add_argument('--mult', '-m', action='store',
                           default=confparser.get('default_options', 'mult'),
                           help=('Set multiplicity for moecule. Default'
                                 'from config file'), type=int)
    argparser.add_argument('--basis', '-b', action='store',
                           default=confparser.get('default_options','basis'),
                           help=('Set basis set for calculation. Default'
                                 'from config file'))
    argparser.add_argument('--functional', '-f', action='store',
                           default=confparser.get('default_options','functional'),
                           help=('Set functional for calculation. Default'
                                 'from config file'))
    argparser.add_argument('--xyz', '-x', action='store', help=
                           ('Set location to read .xyz file for coordinates.'
                           'No default'))
    argparser.add_argument('--symmetry', '-sym', action='store',
                           default=confparser.get('default_options','sym'),
                           help=('Molecular point group. Default'
                                 'from config file'))
    argparser.add_argument('--pcm', '-p', action='store', default=None, help=
                           ('Set PCM solvent simulation on and select solvent.'
                            'Deafult = false'))
    argparser.add_argument('--type', '-t', action='store', default='OPT',
                           choices=('OPT', 'FREQ', 'IRC'), help=('Set type of job, '
                                                          'currently supported'
                                                          'OPT and FREQ.'
                                                          'Default = OPT'))
    argparser.add_argument('--engine', '-e', action='store',
                           default=confparser.get('default_options','engine'),
                           choices=('GAU','GAMESS'), help=('Set which quantum '
                                                         'chemistry packake to '
                                                         'use for calculation. '
                                                         'Default from config '
                                                         'file'))
    
    argparser.add_argument('--walltime', '-w', action='store',
                           default=confparser.get('default_options','walltime'),
                           help=('Set max allowed excution time. Default '
                                 'from config file'))
    argparser.add_argument('--nodes', '-no', action='store',
                           default=confparser.get( 'default_options','nodes'),
                           help=('Set number of nodes to be used.'
                                 'Default from config file'))
    argparser.add_argument('--cpus', '-cp', action='store',
                           default=confparser.get('default_options','cpus'),
                           help=('Set number of cpus per node. Default '
                                 'from config file'))
    argparser.add_argument('--ram', '-r', action='store',
                           default=confparser.get('default_options','ram'),
                           help=('Set amount of ram per node. Default'
                           'from config file'))
    argparser.add_argument('--queue', '-q', action='store',
                           default=confparser.get('default_options','queue'),
                           help='Set queue for job. Default from config file')  

    args = argparser.parse_args()
    
    #Call subroutines for making Gaussian new jobs
    if args.new is True:
        if args.engine == 'GAU':
            try:
                if args.fragment is not None:
                    j = classes.gau_step(args.projectName, args.jobName, args.stepName,
                                         ts=args.tstate,charge=args.charge,
                                         mult=args.mult,xyz=args.xyz, 
                                         fragatom=args.fragment.split()[0], 
                                         frag=args.fragment.split()[1] )
                else:
                    j = classes.gau_step(args.projectName, args.jobName, args.stepName,
                                         ts=args.tstate,charge=args.charge,
                                         mult=args.mult,xyz=args.xyz)
            except IOError, error:
                exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, args.jobName, args.stepName))
            try:
                j.write_inputfile(args.basis, args.functional, args.nodes,
                                  args.cpus, args.ram, args.pcm, args.type)
                print 'Wrote inputfile for: %s-%s-%s' % (args.projectName, args.jobName, args.stepName)
            except IOError, error:
                exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, args.jobName, args.stepName))
        elif args.engine == 'GAMESS':
            try:
                if args.fragment is not None:
                    j = classes.gamess_step(args.projectName, args.jobName, args.stepName,
                                         ts=args.tstate,charge=args.charge,
                                         mult=args.mult,xyz=args.xyz, 
                                         fragatom=args.fragment.split()[0], 
                                         frag=args.fragment.split()[1] )
                else:
                    j = classes.gamess_step(args.projectName, args.jobName, args.stepName,
                                         ts=args.tstate,charge=args.charge,
                                         mult=args.mult,xyz=args.xyz)
            except IOError, error:
                exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, args.jobName, args.stepName))
            try:
                j.write_inputfile(args.basis, args.functional, args.nodes,
                                  args.cpus, args.ram, args.pcm, args.type,
                                  args.symmetry)
                print 'Wrote inputfile for: %s-%s-%s' % (args.projectName, args.jobName, args.stepName)
            except IOError, error:
                exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, args.jobName, args.stepName))
    
    #Call subroutines for sumbitting jobs
    if args.submit is True:
        for job in args.jobName.split():
            for step in args.stepName.split():
                if args.engine == 'GAU':
                    j = classes.gau_step(args.projectName, job, step)
                elif args.engine == 'GAMESS':
                    j = classes.gamess_step(args.projectName, job, step)
                try:
                    j.write_pbsfile(args.walltime, args.nodes,
                                    args.cpus, args.queue)
                except IOError, error:
                    exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, job, step))
                j.submit_job()
                
    #Call subroutines for Analysing job
    if args.analyse is True:
        for job in args.jobName.split():
            for step in args.stepName.split():
                if args.engine == 'GAU':
                    j = classes.gau_step(args.projectName, job, step)
                elif args.engine == 'GAMESS':
                    j = classes.gamess_step(args.projectName, job, step)
                try:
                    j.analyse_job()
                    print j
                except IOError, error:
                    exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, job, step))

    #Call subroutines for calculating reaction parmaters
    if args.reaction is True:
        print 'Project: %s' % args.projectName
        print 'Reaction %s' % args.stepName
        print
        for job in args.jobName.split():
            enthalpies = []
            print 'Substitution: %s' % job
            print
            for step in args.stepName.split('=>'): #Makes a list of anthaplies
                step = step.strip()
                jobs = []
                for sub in step.split('+'):
                    sub = sub.strip()
                    if args.engine == 'GAU':
                        jobs.append(classes.gau_step(args.projectName, job, sub))
                    elif args.engine == 'GAMESS':
                        jobs.append(classes.gamess_step(args.projectName, job, sub))
                    try:
                        jobs[-1].analyse_job()
                    except IOError, error:
                        exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, job, sub))
                try:
                    enthalpies.append(sum(a.enthalpy for a in jobs))
                except AttributeError:
                    exit("AttributeError: Missing Enthalpy for job %s-%s-%s" % (a.projectName, a.jobName, a.stepName))

            activation = []
            reaction = []
            total_reaction = 0
            for count in xrange(0, len(enthalpies)-2, 2):#calcautes thermochem
                activation.append(enthalpies[count+1] - enthalpies[count])
                reaction.append(enthalpies[count+2] - enthalpies[count])
            total_reaction = enthalpies[-1] - enthalpies[0]

            for count in xrange(len(activation)):#outputs therochem
                if len(activation) > 1:
                    print "Reaction Number: %d" % (count+1,)
                print "\tActivation energy = %f kcal/mol" % (float(activation[count])*float(627.509),)
                print "\tReaction enthalpy = %f kcal/mol" % (float(reaction[count])*float(627.509),)
                print
            if len(activation) > 1:
                print "Total reaction enthalpy = %f kcal/mol" % (float(total_reaction)*float(627.509),)
                print
    logging.info('Finished')

    if args.irc is True:
        jobs = []
        for job in args.jobName.split():
            for step in args.stepName.split():
                if args.engine == 'GAU':
                    jobs.append(classes.gau_step(args.projectName, job, step))
                elif args.engine == 'GAMESS':
                    jobs.append(classes.gamess_step(args.projectName, job, step))
                try:
                    jobs[-1].analyse_job()
                except IOError, error:
                    exit("IOError: %s for job %s-%s-%s" % (error, args.projectName, job, step))
            print "Optimized Starting Material Coordinates"
            print jobs[0].format_coords(jobs[0].mol)
            print "IRC Starting Material Coordinates"
            print jobs[1].format_coords(jobs[1].irccoords[0])
            print
            print "Optimized Product Coordinates"
            print jobs[2].format_coords(jobs[2].mol)
            print "IRC Product Coordinates"
            print jobs[1].format_coords(jobs[1].irccoords[-1])
    

#run if not being imported
if __name__ == '__main__':
    main()
