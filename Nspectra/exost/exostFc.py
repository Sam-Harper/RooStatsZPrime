#!/usr/bin/env python
#
# ExoSt: CMS EXOtica Statistics Toolkit
#
# Module: Feldman-Cousins calculator
#
# Author: Gena Kukartsev
#
# January 2011
#
############################################

from ROOT import RooMsgService
from ROOT import RooFit

#from ROOT import RooLinkedList
#from ROOT import RooCmdArg
#from ROOT import RooDataHist

from ROOT import RooStats
from ROOT import TCanvas
from ROOT import TStopwatch
from ROOT import gROOT

from exostAction import exostAction

legend = '[Feldman-Cousins]:'

# module banner
print '\033[1m'+legend, 'Feldman-Cousins calculator module loaded\033[0m'


class exostFc(exostAction):

    def __init__(self):

        # call the base class default constructor
        exostAction.__init__(self)
        
        self.legend = '[Feldman-Cousins]:'
        self.name = 'fc'
        self.config_name = 'Feldman-Cousins Calculator'

        self.register()
        
    #########################################################
    #
    # Configure the action
    #
    #  - 'items' is a dictionary of lists of pairs
    #     from the config file as in exostConfig.items
    #
    def configure(self, ws, items, section, debug = 0):
        print self.legend, 'configuring Feldman-Cousins calculator...'

        _items = items[section]

        # check mandatory parameters

        # dataset name in the workspace
        _dataset_name = self.check_value(items, section, 'data')
        if _dataset_name == None:
            print self.legend, 'ERROR: dataset is not specified, cannot configure'
            return
        else:
            if ws.data(_dataset_name) != None:
                _data = ws.data(_dataset_name)
            else:
                print legend, 'Error: dataset', _dataset_name, 'is not defined, cannot configure'
                return
    
        # model config name in the workspace
        self._model_config_name = self.check_value(items, section, 'model_config')
        if self._model_config_name == None:
            print self.legend, 'ERROR: Model Config name is not specified, cannot configure'
            return
        else:
            if self.check_generic_object(ws, self._model_config_name, debug):
                _mconf = ws.obj(self._model_config_name.strip())
            else:
                print self.legend, 'Error: model config object is not found, cannot configure'
                return
            

        # check optional parameters, set defaults if not specified

        # desired confidence level
        _conf_level = self.check_value(items, section, 'confidence_level')
        if _conf_level == None:
            print self.legend, 'no confidence level specified, setting to 0.95'
            _conf_level = 0.95
                
        # fluctuate number of data entries
        _fluctuate_entries = self.check_value(items, section, 'fluctuate_entries')
        if _fluctuate_entries == None:
            print self.legend, 'number of data entries will be fluctuated'
            _fluctuate_entries = False
        elif _fluctuate_entries == 'True':
            _fluctuate_entries = True
        else:
            _fluctuate_entries = False

        # adaptive sampling
        _adaptive_sampling = self.check_value(items, section, 'adaptive_sampling')
        if _adaptive_sampling == None:
            print self.legend, 'use adaptive sampling'
            _adaptive_sampling = True
        elif _adaptive_sampling == 'True':
            _adaptive_sampling = True
        else:
            _adaptive_sampling = False
            
        # number of bins
        _n_bins = self.check_value(items, section, 'n_bins')
        if _n_bins == None:
            _n_bins = 40
            print self.legend, 'number of bins not specified, using',self._n_bins
                

        # configuring now...

        self._fc = RooStats.FeldmanCousins(_data,_mconf)
        self._fc.SetConfidenceLevel( float(_conf_level) )
        self._fc.FluctuateNumDataEntries(_fluctuate_entries) 
        self._fc.UseAdaptiveSampling(_adaptive_sampling);
        self._fc.SetNBins(int(_n_bins));

        print legend, 'Feldman-Cousins calculator is configured'

        # import the calculator into the Workspace
        #getattr(ws, 'import')(self._fc, 'exostMcmcCalc')
        #print legend, 'Markov chain MC calculator', 'exostMcmcCalc', 'is added to workspace', ws.GetName()
        #print legend, 'done'

        ws.Print()

        # flag that the action is configured successfully
        self.configured = True

        return {'do_plot':False}

    def execute(self, ws, debug = 0):
        print self.legend, 'Feldman-Cousins calculation started...'

        # time the action
        t = TStopwatch()
        t.Start()

        # model config is guaranteed to be here
        # the action won't configure without it
        # and this code will never run unless valid model config is found
        mconf = ws.obj(self._model_config_name)
            
        _poi = mconf.GetParametersOfInterest()

        fcInt = self._fc.GetInterval()

        # stop watch and print how long it took to get the interval
        t.Print()

        # iterate over all parameters of interest and print out the intervals
        # (could there be more than one?)
        _iter = _poi.createIterator()
        while True:
            _poi_name = _iter().GetTitle()
            
            lower_limit = fcInt.LowerLimit( _poi[_poi_name] )
            upper_limit = fcInt.UpperLimit( _poi[_poi_name] )
            
            print self.legend, 'Feldman-Cousins interval for', _poi_name, 'is ['+ \
                  str(lower_limit) + ', ' + \
                  str(upper_limit) + ']'

            if _iter.Next() == None:
                break
            

            
        return (lower_limit, upper_limit)
