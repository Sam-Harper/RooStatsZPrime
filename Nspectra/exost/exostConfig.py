#!/usr/bin/env python

import os
import sys
import ConfigParser
import collections

from OrderedDict import *
#from MultiDict import *

from ROOT import TFile
from ROOT import RooArgSet
from ROOT import RooArgList
from ROOT import RooDataSet
from ROOT import RooDataHist
from ROOT import RooStats
from ROOT import RooFit
from ROOT import RooMsgService
from ROOT import RooFormulaVar

from ROOT import gROOT

from ROOT import *

from exostAction import exostAction

legend = '[ExoSt Config]:'

# module banner
print '\033[1m'+legend, 'config parser loaded\033[0m'


class Config:

    def __init__(self, config_file, debug = 0):
        self.legend = '[ExoSt Config]:'        
        # standard section names
        self.sVariables  = 'Variables'
        self.sData  = 'Data'
        self.sModel = 'Model'
        self.sModelConfig = 'Model Config'
        self.sBayes = 'Bayesian Calculator'

        # self-aware sections
        #   - these sections define which action is theirs
        #   - probably should migrate all sections to this mechanism
        #
        # get config file section title [title] by action name:
        # self.section_title[action_name]
        #
        self.section_title = {}
        #
        # get action name by title:
        # self.section_name[title]
        self.section_name = {}
        #
        # key name
        self.section_key = 'action_name'

        # container for the config
        self.items = {}

        if (config_file != None and os.path.isfile(config_file)):
            print self.legend, 'config file ', config_file, 'exists'
            self.parse(config_file, debug)
    
        else:
            print self.legend, 'no valid config file specified, exiting...'
            sys.exit()


    ###########################################
    #
    # Parse the config file
    #
    #   valid section names: [Variables], [Data], [Model]
    # 
    #   returns a dictionary of items, key is section name
    #
    #   Self-aware sections identify their action with
    #     action_name = name
    #
    def parse(self, filename, debug = 0):
        
        if (debug>2):
            with open(filename) as file:
                print 'parse_config: opened config file', filename
                print ' --------------- Contents of the config file:'
                for line in file:
                    print line.strip('\n')
                print ' --------------- end of config file'
                
        with open(filename) as file:
            # config parser
            
            config = ConfigParser.ConfigParser(dict_type=OrderedDict)
            #config = ConfigParser.ConfigParser(dict_type=OrderedMultiDict)
            config.readfp(file)
            
            sections = config.sections()
            
            for section in sections:
                self.items[section] = config.items(section)

                # find self-aware sections
                _key = self.section_key
                for item in self.items[section]:
                    if item[0] == _key:
                        self.section_title[item[1].strip()] = section
                        self.section_name[section] = item[1].strip()

                        if debug>2:
                            print self.legend, 'found action name', \
                                  self.section_name[section], 'in section', \
                                  '['+section+']'

                #pop the action name definition

            if (debug>2):
                print self.legend, 'config sections from the file:', sections
                print self.legend, 'items per section:', self.items

        return self.items



    ###########################################
    #
    # Helper method, checks if the var is present
    # in the workspace
    #
    def is_var_in_ws(self, ws, var_name):

        _set = ws.allVars()

        if _set.getRealValue(var_name, 12345.654321) == 12345.654321:

            return False

        else:

            return True
        


        



                    
    ###########################################
    #
    # Configure self-aware sections
    #
    # Eventually, these should include all sections
    #
    # Self-aware sections specify the action that is associated
    # with them. The action class must implement
    #   - the configuration part that gets saved to the workspace
    #   - the action part that actually executes the action
    #
    # NOTE that this method must not be called before
    # all actions are registered with the action manager
    #
    def self_aware_sections(self, ws, debug):

        print self.legend, 'configuring self-aware sections now'

        # get the action manager singleton
        mActions = exostAction.instance()

        for action in self.section_title:
            print self.legend, 'configuring action', action, \
                  'from the config file section ['+self.section_title[action]+']'

            mActions.configure(ws,
                               self.items,
                               self.section_title[action],
                               debug,
                               action)



    ##################################################
    #
    # helper method to split config item values
    # using semicolon (;) as the separator
    #
    # it also parses the special case of variables
    #
    # input is a pair (key, value)
    #
    # returns a list of values in the order it is
    # specified in the config file
    #
    def split_values(self, item):
        value_list = []
        
        for value in item[1].strip().strip(';').split(';'):
            _v = value.strip()

            # special case if a variable is defined with [] only
            if _v[0] == '[' and _v[len(_v)-1] == ']':
                _v = item[0].strip()+_v
            
            value_list.append(_v)

        return value_list



    ##################################################
    #
    # helper method to process item options
    # that follow the first value (definition)
    #
    # for example, for fixing a variable:
    #   x = [1.0, 0.0, 2.0];
    #       fix; - this option would be processed here
    #
    # input is
    #   - the Workspace (ws)
    #   - item name (key)
    #   - item option to process (value)
    #
    # saves changes to the workspace
    #
    def process_value(self, ws, key, value, debug = 0):
        if value.strip() == 'fix':
            _var = ws.var(key.strip())
            if _var != None:
                _var.setConstant(True)
                if debug > 0:
                    print self.legend, key, 'is set constant'
            else:
                print self.legend, key, 'is not a variable, cannot fix it!'
                
        if value.strip() == 'float':
            _var = ws.var(key.strip())
            if _var != None:
                _var.setConstant(False)
                if debug > 0:
                    print self.legend, key, 'is set floating'
            else:
                print self.legend, key, 'is not a variable, cannot float it!'


                
                
    ##################################################################
    #
    # process [Variables] section
    #
    # Variables, categories and sets are supposed to be defined here
    #
    # We dump the content of the section in a RooStats::HLFactory
    # card file and make the factory process it.
    #
    # One exception is RooArgSet. Sets are identified and parsed
    # "by hand", and added to the Workspace after all variables are
    # created. This is a temporary solution because named sets were
    # not supported in the factory language until June 22, 2010
    # Eventually, we should switch to using the updated factory
    # language for the named sets.
    #
    def variables(self, ws, debug):
        var_items = self.items[self.sVariables]
        
        legend = '['+self.sVariables+']:'
        
        print legend, 'begin processing...'
        print legend, 'define variables, constants, categories, sets and lists'

        # dictionary of value lists for each item
        # this is more detailed than the original items object
        # because the multiple settings for a given item
        # are already split
        # NOTE: sets are special case and are not added here
        item_values = {}

        # using RooStats::HLFactory to parse the model
        # make a temp file with a card file
        card_file_name = 'hlfCardFile.rs'
        with open(card_file_name, 'w') as cardfile:
            vSets={}
            for item in var_items:
                #
                # special case for sets
                # standard RooFactoryWSTool interface exists in newer versions
                # need to switch at some point
                # for now parsing and adding named sets by hand
                #
                if item[1].strip()[0:3]=='set':
                    set_arg_list = item[1].strip()
                    set_arg_list = set_arg_list[3:len(set_arg_list)]
                    set_arg_list = set_arg_list.strip(' ;()').split(',')

                    set_args=[]
                    for arg in set_arg_list:
                        set_args.append(arg.strip())

                    vSets[item[0].strip()]=set_args

                    if (debug>2):
                        print legend, 'Set arg value is ->'+item[1]+'<-'
                        print legend, 'set args are', set_args
                else:
                    # this is a way to handle multiple settings
                    # for each item: they are separated by semicolon
                    # in the value field of the config item
                    item_values[item[0]] = self.split_values(item)

                    #print >> cardfile,  item[0], '=', item[1]
                    print >> cardfile,  item[0], '=', item_values[item[0]][0]

        # add variables using HLFactory cards
        hlfVerbosity = debug > 2
        hlf = RooStats.HLFactory('HLFactoryExost',
                                 ws,
                                 hlfVerbosity)
        hlf.ProcessCard(card_file_name)

        # process the secondary options (fix/float etc.)
        for _key in item_values:
            n_values = len(item_values[_key])
            for _index in range(0,n_values):
                self.process_value(ws, _key, item_values[_key][_index], debug)


        
        # now add sets by hand
        print legend, 'creating sets now...'

        # fixme experimental
        _sets = {}

        for set_args_key in vSets.keys():
            arg_list = ''

            # fixme experimental
            _sets[set_args_key] = RooArgSet(set_args_key)
            
            for arg in vSets[set_args_key]:
                if len(arg_list)>0:
                    arg_list += ','
                arg_list += arg

                # fixme experimental
                _sets[set_args_key].add(ws.var(arg))
            
            #ws.defineSet(set_args_key, arg_list)
            print legend, 'set', set_args_key, 'created'
            
        print legend, 'done creating sets'

        return {'sets':_sets}


    #
    # helper function: creates RooDataSet and adds to Workspace
    #
    def addData(self, ws, ds_name, item_title, var_set_name, var_set_type, ds_object, weight_var_name = None, debug = 0, sets = None):
        legend = '[exostConfig::addData]:'
        
        arg_set = RooArgSet()
        if (var_set_type == 'set'):
            #arg_set.add(ws.set(var_set_name))
            arg_set.add(sets[var_set_name])
        elif (var_set_type == 'var'):
            arg_set.add(ws.var(var_set_name))
        else:
            print legend, 'error: unknown var_set_type, cannot create dataset', ds_name
            return -1

        #create the dataset
        if weight_var_name == None: #no weight
            if (debug>0):
                print legend, 'no weight branch given'
            ds = RooDataSet(ds_name, item_title, ds_object, arg_set)
        else:
            if (debug>0):
                print legend, 'using weight branch', weight_var_name

            # old and incorrect(?) way of applying weights
            #
            #_var_formula = RooFormulaVar('f'+weight_var_name,
            #                             'f'+weight_var_name,
            #                             '@0',
            #                             RooArgList(ws.var(weight_var_name)))
            #ds = RooDataSet(ds_name, item_title, ds_object, arg_set, _var_formula, weight_var_name)
            
            arg_set.add(ws.var(weight_var_name))
            ds = RooDataSet(ds_name, item_title, ds_object, arg_set, "1>0", weight_var_name)
        
        # import the datahist. Note workaround 'import' being a reserved word
        getattr(ws, 'import')(ds)

        #Fixme: debug
        #ds.tree().SaveAs("bTest.root")



    #
    # helper function: creates RooDataSet from ASCII file and adds to Workspace
    #
    # variables in the set correspond to columns in the file
    #
    def addAsciiData(self, ws, ds_name, item_title, var_set_name, var_set_type, ds_file_name, weight_var_name = None, debug = 0, sets = None):
        legend = '[exostConfig::addAsciiData]:'
        
        arg_set = RooArgSet()
        if (var_set_type == 'set'):
            #_var_set = ws.set(var_set_name)
            _var_set = sets[var_set_name]
            arg_set.add(_var_set)

            if _var_set.getSize() != 1:
                print legend, 'Error: too many or too few columns in the input ASCII file'
                print legend, 'Error: Smart program as I am, I can only handle ASCII files'
                print legend, 'Error: with exactly one column at the moment.'
                print legend, 'Error: Support for multiple columns will be implemented'
                print legend, 'Error: eventually, contact the developers.'
                print legend, 'Error: Better yet, switch to ROOT input files,'
                print legend, 'Error: all the cool kids are doing that!'

                return -1
            
        elif (var_set_type == 'var'):
            arg_set.add(ws.var(var_set_name))
        else:
            print legend, 'error: unknown var_set_type, cannot create dataset', ds_name
            return -1

        #create the dataset
        if weight_var_name == None: #no weight
            if (debug>0):
                print legend, 'no weight variable given'

            _arglist = RooArgList(arg_set)
            #ds = RooDataSet()ds_name, item_title)
            ds = RooDataSet.read(ds_file_name, _arglist)
            ds.SetName(ds_name)
            ds.SetTitle(item_title)
        else:
            if (debug>0):
                print legend, 'using variable', weight_var_name, 'as weight'
        
                print legend, 'Error: Smart program as I am, I cannot handle'
                print legend, 'Error: weights when loading from ASCII files yet.'
                print legend, 'Error: Support for weights from ASCII files will'
                print legend, 'Error: be implemented eventually, contact the developers.'
                print legend, 'Error: Better yet, switch to ROOT input files,'
                print legend, 'Error: all the cool kids are doing that!'

                return -1


        # import the datahist. Note workaround 'import' being a reserved word
        getattr(ws, 'import')(ds)

        #Fixme: debug
        #ds.tree().SaveAs("bTest.root")



    #
    # helper function: creates RooDataHist and adds to WS
    #
    def addHistData(self, ws, ds_name, item_title, var_set_name, var_set_type, ds_object, sets = None):
        legend = '[exostConfig::addHistData]:'
        
        arg_list = RooArgList()
        if (var_set_type == 'set'):
            #arg_list.add(ws.set(var_set_name))
            arg_list.add(sets[var_set_name])
        elif (var_set_type == 'var'):
            arg_list.add(ws.var(var_set_name))
        else:
            print legend, 'error: unknown var_set_type, cannot create dataset', ds_name
            return -1

        #create the hist dataset
        ds = RooDataHist(ds_name, item_title, arg_list, ds_object, 1)
        
        # import the datahist. Note workaround 'import' being a reserved word
        getattr(ws, 'import')(ds)


    #
    # helper function: create RooAbsData
    # by combining several existing datasets
    # and add to WS
    #
    def addCombinedData(self, ws, ds_name, item_args, sets = None):
        legend = '[exostConfig::addHistData]:'

        print legend, 'creating a combined dataset'

        if len(item_args) < 8:
            print legend, 'too few arguments specified for a combined dataset, cannot create'
            return -1

        ds_title = item_args[1]

        v_name = item_args[2]

        cat_name = item_args[3]
        _cat = ws.cat(cat_name)

        # this is a workaround for CMSSW
#        if _cat == None:
#            print legend, 'Workspace contains no category named', cat_name, 'cannot create dataset'
#            return -1
        if _cat:
            print legend, 'Workspace contains category named', cat_name
        else:
            print legend, 'Workspace contains no category named', cat_name, 'cannot create dataset'
            return -1
        
        #var_set = ws.set(v_name)
        var_set = sets[v_name]

        # workaround CMSSW issues (commented out part works fine standalone)
        #if var_set == None:
        #    print legend, 'Workspace contains no set named', v_name, 'searching in variables...'
        #    var_ = ws.var(item_args[2])
        #
        #    if var_ == None:
        #        print legend, 'Workspace contains no variable named', v_name, 'cannot create dataset'
        #        return -1
        #
        #    else:
        #        print legend, 'Workspace contains a variable named', v_name, 'will use it to define dataset'
        #        var_set = RooArgSet(var_)

        if var_set:
                print legend, 'Workspace contains a set named', v_name, 'will use it'            
        else:
            print legend, 'Workspace contains no set named', v_name, 'searching in variables...'
            var_ = ws.var(item_args[2])

            if var_:
                print legend, 'Workspace contains a variable named', v_name, 'will use it to define dataset'
                var_set = RooArgSet(var_)
            else:
                print legend, 'Workspace contains no variable named', v_name, 'cannot create dataset'
                return -1

        # get pairs of state name - dataset name
        _state = []
        _dsname = []
        _ds = []

        _state.append(item_args[4])
        _dsname.append(item_args[5])
        _state.append(item_args[6])
        _dsname.append(item_args[7])

        _index = 0
        for _name in _dsname:
            _ds.append(ws.data(_name))
            if _ds[_index] == None:
                print legend, 'Workspace contains no dataset named', _name, 'cannot create dataset'
                return -1
                
            _index += 1

        #create the dataset
        ds = RooDataSet(ds_name, ds_title, var_set,
                        RooFit.Index(_cat),
                        RooFit.Import(_state[0], _ds[0]),
                        RooFit.Import(_state[1], _ds[1])
                        )
                        #,RooFit.WeightVar('n'))
        
        # import the dataset. Note workaround, 'import' being a reserved word
        getattr(ws, 'import')(ds)
        



    ###########################################
    #
    # process [Data] section
    #
    # format of the data string:
    #
    #    dataset_name = source_type name title filename var_set
    # or
    #    dataset_name = source_type(name, title, filename, var_set);
    #
    #      source_type: can be root, ascii or combine
    #      name:        name of tree or histogram that is the data source
    #      title:       must be delimited by "" and can be any string
    #      filename:    name of the file with the data
    #      var_set:     name of a variable or a set of variables defined in
    #                  the [Variables] section
    #
    # example:
    #
    #   drellyanMC = root dyTree "Drell Yan MC" data/dy.root vars
    #
    def data(self, ws, debug = 0, sets = None):

        if self.sData not in self.items:
            if debug>0:
                print self.legend, 'no Data section defined in the config file'
                print self.legend, 'no data to load'
            return

        data_items = self.items[self.sData]
        
        legend = '['+self.sData+']:'
        
        # allowed file types
        file_types = ['root', 'ascii', 'combine']
        
        # put data into the workspace
        if (debug>0):
            print legend, 'Loading data to the workspace', ws.GetName()

        for item in data_items:
            if (debug>2):        
                print legend, item[1]

            # clean up arguments left from the previous config string
            item_args = []

            # identify the delimitor type
            #
            #    can be space-delimited string (legacy)
            # or
            #    python-style with brakets and comma-separated parameters
            if item[1].find('(')>0: # python-style

                # get the source_type
                item_args.append( item[1].strip()[0:item[1].find('(')] )

                # get the insides of the brackets
                item_args_args = item[1].strip()[item[1].find('(')+1:item[1].find(')')].split(',')

                # strip spaces, quotes and append to item_args
                for arg in item_args_args:

                    item_args.append( arg.strip().strip('\"') )

                item_title = item_args[2]



            else:    # legacy space delimiter style

                # cut out the dataset title and split other arguments
                first_quote = item[1].find('\"')
                if (first_quote>0):
                    item_args = item[1].split('\"')
                    item_title = item_args[1]
                    #item_args = (item_args[0]+' '+item_args[2]).split()
                    item_args_new = item_args[0].split()
                    item_args_new.append(item_args[1])
                    for _item in item_args[2].split():
                        item_args_new.append(_item)
                    item_args = item_args_new

                # title is one word without quotes:
                else:
                    item_title = item[1].split()[2]
                    item_args = item[1].split()
                
            
            if (debug>2):        
                print legend, 'dataset title:', item_title
                print legend, 'dataset definition arguments:', item_args
        
            # validate input
            # check if the dataset with this name exists already
            ds_name = item[0].strip()
            if (ws.data(ds_name)!=None):
                print legend, 'Dataset', item[0], 'already exists'
                return -1

            # known file type?
            file_type = item_args[0].strip()
            if (file_type not in file_types):
                print legend, 'Error: unknown data source:', item_args[0]
                return -1


            #
            # combine existing datasets into a categorised dataset
            # this is needed, for example, for simultaneous fits
            #
            if file_type == 'combine':

                self.addCombinedData(ws, ds_name, item_args, sets)



            else: # root or ascii
                
                # remove title from the list
                # (historical reasons, to preserve numbering)
                item_args.pop(2)
                
                # get the object name (hist or tree)
                object_name = item_args[1].strip()

                # input file exists?
                file_name = item_args[2].strip()
                if os.path.isfile(file_name)==False:
                    print legend, 'Error: file ', item_args[2], 'does not exist'
                    return -1

                # check if the var_set is valid
                # if valid, determine type: set or var
                var_set_name = item_args[3].strip()
                #fixme: experimental, sets are outside ws
                #if ws.set(var_set_name)!=None:
                if var_set_name in sets:
                    if (debug>2):
                        print legend, 'Workspace does have a set named', var_set_name
                    var_set_type = 'set'
                
                elif ws.var(var_set_name)!=None:
                    if (debug>2):
                        print legend, 'Workspace does have a variable named', var_set_name
                    var_set_type = 'var'
                else:
                    print legend, 'Error: set or variable', var_set_name, 'is not defined'
                    return -1

                # weights branch (optional)
                weight_var_name = None
                if len(item_args) > 4:
                    weight_var_name = item_args[4].strip()
                    
                    # workaround for checking if variable exists in the workspace
                    #if (ws.var(weight_var_name)==None):
                    if self.is_var_in_ws(ws, weight_var_name) == False:
                        print legend, 'Weight variable', item_args[4], 'does not exist, exiting...'
                        return -1
                
            
                if (debug>0):        
                    print legend, 'dataset parameters seem ok, trying to create dataset now...'


            # load from ROOT file
            if file_type == 'root':

                # open the file
                ds_file = TFile(file_name, "READ")
                ds_object = ds_file.Get(object_name)
                ds_class = ds_object.ClassName()

                if (debug>0):        
                    print legend, 'Dataset object is of type', ds_class

                if (ds_class[0:2]=='TH'):
                    print legend, 'will create RooDataHist'
                    self.addHistData(ws, ds_name, item_title, var_set_name, var_set_type, ds_object, sets)

                elif (ds_class[0:5]=='TTree' or ds_class[0:7]=='TNtuple'):
                    print legend, 'will create RooDataSet'
                    self.addData(ws, ds_name, item_title, var_set_name, var_set_type, ds_object, weight_var_name, debug, sets)
                
                else:
                    print legend, 'Error: object', object_name,
                    'from file', file_name, 'is of unknown type', ds_class
                
            # load from ASCII file
            if file_type == 'ascii':
                print legend, 'will create RooDataSet'
                self.addAsciiData(ws, ds_name, item_title, var_set_name, var_set_type, file_name, weight_var_name, debug, sets)
                
        
        print legend, 'done.\n'




    ###########################################
    #
    # helper function fro the Model section:
    #
    # prints an explanation of a RooFIt warning
    # that is due to a nonstandard CMSSW-like
    # ROOT distribution
    #
    def modelAutoImportWarningExplanation(self, legend, debug = 0):
        from ROOT import gSystem
        #if 'ROOFIT_INCLUDE' in os.environ.keys():
        if 'CMSSW_VERSION' in os.environ.keys():
            if debug > 0:
                print '\033[1m' # bold font on
                print legend, '******************************************'
                print legend, ''
                print legend, 'It seems like you are running in a nonstandard ROOT environment,'
                print legend, 'likely, the one that comes with CMSSW. You may encounter warning'
                print legend, 'messages like the ones given as examples below. They are harmless'
                print legend, 'and have to do with the fact that RooFit is not in the standard'
                print legend, 'ROOT location. RooWorkspace gets confused and attempts to import'
                print legend, 'the standard RooFIt Class code. It fails, which works perfectly'
                print legend, 'since the standard RooFit code does not need to get imported.'
                print legend, ''
                print legend, 'It is safe to disregard the following types of warning messages.'
                print legend, ''
                print legend, '  [#0] WARNING:ObjectHandling -- RooWorkspace::autoImportClass(myWS)'
                print legend, '  WARNING Cannot access code of class RooAbsReal because implementation'
                print legend, '  file roofitcore/src/RooAbsReal.cxx is not found in current directory'
                print legend, '  nor in $ROOTSYS, nor in the search path . To fix this problem add the'
                print legend, '  required directory to the search path using' 
                print legend, '  RooWorkspace::addClassImplDir(const char* dir)'
                print legend, ''
                print legend, '  [#0] WARNING:ObjectHandling -- RooWorkspace::import(myWS)'
                print legend, '  WARNING: problems import class code of object RooRealVar::invmass,'
                print legend, '  reading of workspace will require external definition of class'
                print legend, ''
                print legend, '******************************************'
                print '\033[0m' # bold font off
            else:
                print '\033[1m' # bold font on
                print legend, '******************************************'
                print legend, ''
                print legend, 'It seems like you are running in a nonstandard ROOT environment,'
                print legend, 'likely, the one that comes with CMSSW.'
                print legend, 'Certain harmless RooFit warning messages are turned off. If you wish'
                print legend, 'to see them, turn on the verbosity level to greater than zero.'
                print legend, 'You can do that with the -d command line parameter, e.g.'
                print legend, ''
                print legend, '  ./exost [...] -d 2'
                print legend, ''
                print legend, '******************************************'
                print '\033[0m' # bold font off
                


    
    ###########################################
    #
    # process [Model] section
    #
    # define all PDFs here, including the full likelihood
    # you may also define additional variables here
    #
    # HLFactory language is understood in this section
    #
    def model(self, ws, debug = 0):

        if self.sModel not in self.items:
            if debug>0:
                print self.legend, 'no Model section defined in the config file'
                print self.legend, 'no model will be created'
            return

        legend = '['+self.sModel+']:'

        # FIXME: find a better way to load custom PDFs
        #gROOT.ProcessLine(".L dijetQstarPdf.cxx+");
        #gROOT.ProcessLine(".L Qstar_qg.cxx+");
        #gROOT.ProcessLine(".L Qstar_qg_2.cxx+");
        #gROOT.ProcessLine(".L Jacobian_mt.cxx+");

        model_items = self.items[self.sModel]
        
        if (debug>0):
            print 'Loading model to the workspace', ws.GetName()

        # using RooStats::HLFactory to parse the model
        # make a temp file with a card file
        card_file_name = 'hlfCardFile.rs'
        with open(card_file_name, 'w') as cardfile:
            for item in model_items:
                print >> cardfile,  item[0], '=', item[1]

        hlfVerbosity = debug > 0
        #hlf = ROOT.RooStats.HLFactory('HLFactoryExost', card_file_name, True)
        hlf = RooStats.HLFactory('HLFactoryExost',
                                      ws,
                                      hlfVerbosity)

        self.modelAutoImportWarningExplanation(legend, debug)

        current_messaging_level = RooMsgService.instance().globalKillBelow() ;
        if debug <= 0:
            RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR) ;

        hlf.ProcessCard(card_file_name)

        print legend, 'importing nonstandard class code into the workspace...'
        ws.importClassCode();
        #ws.importClassCode('Qstar_qg', true);
        print legend, '...code import done'

        # return messaging to the original level
        RooMsgService.instance().setGlobalKillBelow(current_messaging_level) ;


    ###########################################
    #
    # process [Model Config] section
    #
    # define here
    #   - model (model = model_pdf)
    #   - parameters of interest (poi = s)
    #   - nuisance parameters: must be a set in a Workspace
    #
    #
    def model_config(self, ws, debug = 0, sets = None):

        if self.sModelConfig not in self.items:
            if debug>0:
                print self.legend, 'no Model Config section in the config file'
                print self.legend, 'cannot create Model Config'
            return

        legend = '['+self.sModelConfig+']:'

        # getting items as a dictionary
        model_config_items = {}
        model_config_items.update(self.items[self.sModelConfig])
        
        print legend, 'Configuring the model... '

        # get the desired name for the model config object
        _model_config_name = self.check_value(self.sModelConfig, 'name')
        if _model_config_name == -1:
            print self.legend, 'using default name for Model Config: exostModelConfig'
            _model_config_name = 'exostModelConfig'
        else:
            print self.legend, 'creating Model Config with name', _model_config_name
            

        # create the ModelConfig object and associate it with the Workspace
        model_config = RooStats.ModelConfig(_model_config_name, _model_config_name)
        model_config.SetWorkspace(ws)

        # check if the model is defined
        if 'model' in model_config_items.keys():
            _model = model_config_items['model']
            print legend, 'The full likelihood model:', _model
            #model_config.SetPdf(ws.pdf(_model))
            model_config.SetPdf(_model)
        else:
            print legend, 'Error: the model PDF is not specified'
            print legend, 'Error: the model config is invalid'
            print legend, 'Error: the model config is not added to the workspace'
            return -1
            

        # check if the observables are specified
        observables_valid = False
        if 'observables' in model_config_items.keys():
            observables_valid = True
            
            _observables = model_config_items['observables']
            print legend, 'Observables:', _observables
            # fixme experimental
            #if ws.set(_observables) != None:
            #    model_config.SetObservables(ws.set(_observables))
            if sets[_observables] != None:
                model_config.SetObservables(sets[_observables])
            else:
                observables_valid = False
                print legend, 'Warning: observables set', _observables, 'is not defined'
        else:
            observables_valid = False
            print legend, 'Warning: no observables specified'

        # observables is an optional field in ModelConfig, no need to invalidate
        #if observables_valid == False:
        #    print legend, 'Error: the model config is invalid'
        #    print legend, 'Error: the model config is not added to the workspace'
        #    return -1
            

        # check if the global global_observables are specified
        global_observables_valid = False
        if 'global_observables' in model_config_items.keys():
            global_observables_valid = True
            
            _global_observables = model_config_items['global_observables']
            print legend, 'Global observables:', _global_observables
            # fixme experimental
            #if ws.set(_global_observables) != None:
            #    model_config.SetGlobalObservables(ws.set(_global_observables))
            if sets[_global_observables] != None:
                model_config.SetGlobalObservables(sets[_global_observables])
            else:
                global_observables_valid = False
                print legend, 'Warning: global observables set', _global_observables, 'is not defined'
        else:
            global_observables_valid = False
            print legend, 'Warning: no global observables specified'

        #if global_observables_valid == False:
        #    print legend, 'Error: the model config is invalid'
        #    print legend, 'Error: the model config is not added to the workspace'
        #    return -1
            

        # check if the parameters of interest are specified
        poi_valid = False
        if 'poi' in model_config_items.keys():
            poi_valid = True
            
            _poi = model_config_items['poi']
            print legend, 'Parameters of interest:', _poi
            # fixme experimental
            #if ws.set(_poi) != None:
            #    model_config.SetParametersOfInterest(ws.set(_poi))
            if sets[_poi] != None:
                model_config.SetParametersOfInterest(sets[_poi])
            else:
                poi_valid = False
                print legend, 'Error: POI set', _poi, 'is not defined'
        else:
            poi_valid = False
            print legend, 'Error: no parameters of interest specified'

        if poi_valid == False:
            print legend, 'Error: the model config is invalid'
            print legend, 'Error: the model config is not added to the workspace'
            return -1
            

        # check if the nuisance parameters are specified
        nuis_valid = False
        if 'nuisance_parameters' in model_config_items.keys():
            nuis_valid = True
            
            _nuis = model_config_items['nuisance_parameters']
            print legend, 'Nuisance parameters:', _nuis
            #fixme experimental
            #if ws.set(_nuis) != None:
            #    model_config.SetNuisanceParameters(ws.set(_nuis))
            if sets[_nuis] != None:
                model_config.SetNuisanceParameters(sets[_nuis])
            else:
                nuis_valid = False
                print legend, 'Error: Nuisance parameter set', _nuis, 'is not defined'
        else:
            nuis_valid = False
            print legend, 'Error: no nuisance parameters specified'

        if nuis_valid == False:
            print legend, 'Error: the model config is invalid'
            print legend, 'Error: the model config is not added to the workspace'
            return -1
            

        # check if a prior PDF is specified
        prior_valid = False
        if 'prior' in model_config_items.keys():
            prior_valid = True
            
            _prior = model_config_items['prior']
            print legend, 'Prior:', _prior
            if ws.pdf(_prior) != None:
                model_config.SetPriorPdf(ws.pdf(_prior))
            else:
                prior_valid = False
                print legend, 'Error: prior', _prior, 'is not defined'
        else:
            prior_valid = False
            print legend, 'Error: no prior PDF specified'

        if prior_valid == False:
            print legend, 'Error: the model config is invalid'
            print legend, 'Error: the model config is not added to the workspace'
            return -1
            

        # import the model config
        getattr(ws, 'import')(model_config, _model_config_name)

        print legend, 'Model config is added to workspace ', ws.GetName()
        print legend, 'done'



    #
    # check and report a config value
    #
    def check_value(self, section, key):

        legend = '['+section+']:'

        # getting items as a dictionary
        _items = {}
        _items.update(self.items[section])

        if key in _items.keys():
            _value = _items[key]
            print legend, key+':', _value
            return _value
        else:
            print legend, 'Warning: key', key, 'is not specified'
            return -1


    ###########################################
    #
    # process [Bayesian Calculator] section
    #
    # define here
    #   - confidence_level = 0.90
    #   - posterior_plot = True
    #
    #
    def bayesian(self, ws, debug = 0):

        if self.sBayes not in self.items:
            if debug>0:
                print self.legend, 'no Bayesian Calculator section in the config file'
                print self.legend, 'cannot configure Bayesian calculator'
            return {'status':'fail'}

        legend = '['+self.sBayes+']:'

        _name = 'exostBayes'

        # getting items as a dictionary
        bayes_items = {}
        bayes_items.update(self.items[self.sBayes])
        
        print legend, 'Configuring Bayesian calculator... '

        # check if a CL is specified
        m_conf_name = self.check_value(self.sBayes, 'model_config')
        if m_conf_name == -1:
            print legend, 'Error: model config is not specified'
            print legend, 'Error:', self.sBayes, 'cannot be configured'
            return {'status':'fail'}
            

        # check if data is specified
        data_name = self.check_value(self.sBayes, 'data')

        data_valid = False
        if data_name != -1:
            data_valid = True
            if ws.data(data_name) != None:
                data = ws.data(data_name)
            else:
                data_valid = False
                print legend, 'Error: dataset', data_name, 'is not defined'

        if data_valid == False:
            print legend, 'Error:', self.sBayes, 'cannot be configured'
            return {'status':'fail'}


        # check if a CL is specified
        conf_level = self.check_value(self.sBayes, 'confidence_level')
        if conf_level == -1:
            print legend, 'Warning: desired confidence level is not specified, setting to 0.90'
            conf_level = 0.90


        # check if a posterior plot is requested
        post_plot = self.check_value(self.sBayes, 'posterior_plot')
        b_post_plot = False
        if post_plot == -1:
            print legend, 'will generate a posterior plot'
            b_post_plot = True
        elif post_plot == 'True':
            b_post_plot = True
        else:
            b_post_plot = False

        # plot format
        if post_plot:
            plot_format = self.check_value(self.sBayes, 'plot_format')
            if plot_format == -1:
                plot_format = 'png'


        # to suppress messgaes when pdf goes to zero
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

        #mconf = ws.obj('exostModelConfig')
        mconf = ws.obj(m_conf_name)

        if mconf == None:
            if debug > 0:
                print self.legend, 'fail to get Model Config'
                print self.legend, 'unable to configure Bayesian calculator'
                print self.legend, 'calculator will not be imported into the workspace'
            return {'status':'fail'}
                    
        bCalc = RooStats.BayesianCalculator(data, mconf)
        bCalc.SetConfidenceLevel(float(conf_level))


        ############# FIXME: debug test
        #bInt = bCalc.GetInterval()
        #print legend, 'DEBUG3*******************************'
        # 
        #cl = bCalc.ConfidenceLevel()
        #print legend,  str(cl)+'% CL central interval: [', bInt.LowerLimit(), ' - ', bInt.UpperLimit(), ']'
        #print legend, 'or', str(cl+(1.0-cl)/2), '% CL limits'
        ##############################


        #size = 1.0 - float(conf_level)
        #bCalc.SetTestSize(size)
        
        # import the calculator into the Workspace
        getattr(ws, 'import')(bCalc, _name)

        print legend, 'Bayesian calculator', _name, 'is configured and added to workspace', ws.GetName()
        print legend, 'done'

        return {'do_posterior_plot':b_post_plot,
                'status':'success',
                'model_config_name':m_conf_name,
                'plot_format':plot_format}
