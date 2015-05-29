#!/usr/bin/env python

import os
import sys
import subprocess
import datetime

class exostEnv:

    def __init__(self, options):

        self.options = options

        self.legend = '[Environment]:'

        # try to autoconfigure with this release
        self.cmssw_release = None

        # default CMSSW architecture
        self.arch = 'slc5_ia32_gcc434'

        # get ROOTSYS if it is set, None otherwise
        if 'ROOTSYS' in os.environ:
            self.rootsys = os.environ['ROOTSYS']
        else:
            self.rootsys = None
            
        # host could be cmslpc, lxplus or None (default)
        self.host = self.guess_host()

        self.configure()



    # try to guess whether on cmslpc, lxplus
    def guess_host(self):

        # forced standalone
        #if self.options.standalone:
        #    return None

        #if standalone is not forced, guess
        if 'HOSTNAME' in os.environ:
            self.hostname_ = os.environ['HOSTNAME']
            if 'cmslpc' in self.hostname_:
                print self.legend, 'I am guessing that you are on cmslpc...'
                return 'cmslpc'
            elif 'lxplus' in self.hostname_:
                print self.legend, 'I am guessing that you are on lxplus...'
                return 'lxplus'
            else:
                return None
        else:
            return None
        


    def get_cmssw_dict(self):
        # get a dictionary of CMSSW versions using scram
        #   key -   CMSSW version name
        #   value - CMSSW_RELEASE_BASE
        
        cmssw_ = subprocess.Popen(['scram', 'list'],
                                  stdout=subprocess.PIPE).communicate()[0]

        self.cmssw_dict={}
        _next = False # whether the next line contains the path of interest
        for line in cmssw_.splitlines():
            if len(line.strip())!=0:
                split_line_ = line.split()
                #print 'DEBUG:', split_line_
                if split_line_[0] == 'CMSSW':
                    _version = split_line_[1]
                    _next = True
                    
                if _next:
                    self.cmssw_dict[_version] = split_line_[1]

        print self.legend, 'Dictionary of available CMSSW releases:'
        print self.legend, self.cmssw_dict
                
        return self.cmssw_dict


    def configure_cmssw(self):
        # setup the environment using CMSSW
        # need to define:
        #   ROOT_INCLUDE = ${ROOTSYS}/include
        #   ROOFIT_INCLUDE = ROOFIT_INCLUDE ${CMS_PATH}/$arch/lcg/roofit/$roofit_release/include
        # for that, we also need:
        #   arch = slc5_ia32_gcc434  (value as example)
        #   roofit_release = 5.27.04 (value as example)

        print self.legend, 'setting up ExoSt using CMSSW environment'

        # find out CMSSW release
        if 'CMSSW_VERSION' in os.environ.keys():
            self.cmssw_release = os.environ['CMSSW_VERSION']
            print self.legend, 'CMSSW version:', self.cmssw_release

            # find out CMSSW path
            if 'CMS_PATH' in os.environ.keys():
                self.cms_path = os.environ['CMS_PATH']
                print self.legend, 'CMS_PATH:', self.cms_path
            else:
                print self.legend, 'CMS_PATH is not initialized, exiting...'
                sys.exit(-1)

            # check if scram is available,
            # get cmssw_release_base from it or set default
            scram_ = subprocess.Popen(['which', 'scram'],
                                      stdout=subprocess.PIPE).communicate()[0]
            if 'Command not found' in scram_:
                print self.legend, '* Error: SCRAM tool is not available'
                print self.legend, '* Error: Did you set the CMSSW environment properly?'
                print self.legend, '* Error: exiting...'
                sys.exit(-1)
                
            else:
                print self.legend, 'found SCRAM - good'
                if 'SCRAM_ARCH' in os.environ:
                    self.arch = os.environ['SCRAM_ARCH']
                    print self.legend, 'determined architecture as', self.arch
                else:
                    print self.legend, 'using default architecture:', self.arch

            # find RooFit version
            print self.legend, 'querying SCRAM for RooFit version...'
            roofit_ = subprocess.Popen(['scram', 'tool', 'info', 'roofit'],
                                       stdout=subprocess.PIPE).communicate()[0]
            for line in roofit_.splitlines():
                if 'Version' in line:
                    self.roofit_version = line.split(':')[1].strip()
            print self.legend, 'found RooFit version', self.roofit_version
        
            
            # set ROOT_INCLUDE
            os.environ['ROOT_INCLUDE'] = self.rootsys.rstrip('/')+'/include'
            print self.legend, 'ROOT include dir set to', os.environ['ROOT_INCLUDE']

            if self.host=='cmslpc':
                #setenv ROOFIT_INCLUDE ${CMS_PATH}/$arch/lcg/roofit/$roofit_release/include
                os.environ['ROOFIT_INCLUDE'] = self.cms_path.rstrip('/')+'/'+self.arch+'/lcg/roofit/'+self.roofit_version+'/include'
            else:
                #setenv ROOFIT_INCLUDE ${CMS_PATH}/sw/${SCRAM_ARCH}/lcg/roofit/$roofit_release/include
                #os.environ['ROOFIT_INCLUDE'] = self.cms_path.rstrip('/')+'/sw/'+self.arch+'/lcg/roofit/'+self.roofit_version+'/include'
                os.environ['ROOFIT_INCLUDE'] = self.cms_path.rstrip('/')+'/'+self.arch+'/lcg/roofit/'+self.roofit_version+'/include'

            print self.legend, 'RooFit include dir set to', os.environ['ROOFIT_INCLUDE']

            print self.legend, 'ExoSt environment with CMSSW is set'


        elif self.host=='lxplus':
            print self.legend, 'CMSSW_VERSION is not initialized, trying standalone setup...'
            print self.legend, '*** Make sure that you have sourced ExoSt/setup/lxplus_standalone.[c]sh'
            print self.legend, '*** If you have not, we will probably crash about ... now'
            
            print self.legend, 'ExoSt lxplus standalone environment is set'


        elif self.host=='cmslpc':
            print self.legend, 'CMSSW_VERSION is not initialized, trying standalone setup...'
            print self.legend, '*** Make sure that you have sourced ExoSt/setup/cmslpc_standalone_setup.[c]sh'
            print self.legend, '*** If you have not, we will probably crash about ... now'
            
            print self.legend, 'ExoSt lxplus standalone environment is set'


        else:
            print self.legend, 'CMSSW_VERSION is not initialized, exiting...'
            sys.exit(-1)





    # setup the standalone environment
    def configure_standalone(self):
        print self.legend, 'setting up standalone ExoSt environment...'

        #os.environ['LD_LIBRARY_PATH'] = self.rootsys.rstrip('/')+'/lib'
        os.environ['LD_LIBRARY_PATH'] = self.rootsys.rstrip('/')+'/lib:'+os.environ['LD_LIBRARY_PATH']
        print self.legend, 'LD_LIBRARY_PATH set to', os.environ['LD_LIBRARY_PATH']

        print self.legend, 'adding LD_LIBRARY_PATH to the list of module locations'
        sys.path.append(self.rootsys.rstrip('/')+'/lib')

        print self.legend, 'ExoSt standalone environment is set'


    # setup the environment
    def configure(self):

        # check if ROOTSYS is set
        if self.rootsys == None:
            print self.legend, 'ROOTSYS is not set, cannot proceed, exiting...'
            sys.exit(-1)

        # special case of cmslpc and lxplus
        if self.host=='cmslpc' or self.host=='lxplus':

            self.configure_cmssw()

        else:

            self.configure_standalone()

