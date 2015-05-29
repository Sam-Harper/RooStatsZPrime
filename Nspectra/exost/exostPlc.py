#!/usr/bin/env python
#
# ExoSt: CMS EXOtica Statistics Toolkit
#
# Module: Profile Likelihood calculator
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

legend = '[Profile Likelihood]:'

# module banner
print '\033[1m'+legend, 'Profile Likelihood calculator module loaded\033[0m'


class exostPlc(exostAction):

    def __init__(self):

        # call the base class default constructor
        exostAction.__init__(self)
        
        self.legend = '[Profile Likelihood]:'
        self.name = 'plc'
        self.config_name = 'Profile Likelihood Calculator'

        self.register()
        
    #########################################################
    #
    # Configure the action
    #
    #  - 'items' is a dictionary of lists of pairs
    #     from the config file as in exostConfig.items
    #
    def configure(self, ws, items, section, debug = 0):
        print self.legend, 'configuring Profile Likelihood calculator...'

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
                
        # whether to make the posterior plot or not
        self._scan_plot = self.check_value(items, section, 'make_scan_plot')
        if self._scan_plot == None:
            print self.legend, 'likelihood intervalscan plot is not requested'
            self._scan_plot = False
        elif self._scan_plot == 'True':
            self._scan_plot = True
        else:
            print self.legend, 'invalid assignment make_scan_plot =', self._scan_plot
            print self.legend, 'likelihood interval scan plot will not be made'
            self._scan_plot = False

        # plot format
        if self._scan_plot:
            self._plot_format = self.check_value(items, section, 'plot_format')
            if self._plot_format == None:
                print self.legend, 'no plot format specified, setting to PNG'
                self._plot_format = png
            
                

        # configuring now...

        self._plc = RooStats.ProfileLikelihoodCalculator(_data,_mconf)
        self._plc.SetConfidenceLevel( float(_conf_level) )

        print legend, 'Profile Likelihood calculator is configured'

        # import the calculator into the Workspace
        #getattr(ws, 'import')(self._plc, 'exostMcmcCalc')
        #print legend, 'Markov chain MC calculator', 'exostMcmcCalc', 'is added to workspace', ws.GetName()
        #print legend, 'done'

        ws.Print()

        # flag that the action is configured successfully
        self.configured = True

        return {'do_scan_plot':self._scan_plot}

    def execute(self, ws, debug = 0):
        print self.legend, 'Profile Likelihood calculation started...'

        # time the action
        t = TStopwatch()
        t.Start()

        # model config is guaranteed to be here
        # the action won't configure without it
        # and this code will never run unless valid model config is found
        mconf = ws.obj(self._model_config_name)
            
        _poi = mconf.GetParametersOfInterest()

        plcInt = self._plc.GetInterval()

        # stop watch and print how long it took to get the interval
        t.Print()

        # iterate over all parameters of interest and print out the intervals
        # (could there be more than one?)
        _iter = _poi.createIterator()
        while True:
            _poi_name = _iter().GetTitle()
            
            lower_limit = plcInt.LowerLimit( _poi[_poi_name] )
            upper_limit = plcInt.UpperLimit( _poi[_poi_name] )
            
            print self.legend, 'Profile Likelihood interval for', _poi_name, 'is ['+ \
                  str(lower_limit) + ', ' + \
                  str(upper_limit) + ']'

            if _iter.Next() == None:
                break
            

        if self._scan_plot:
            # draw scan plot

            print self.legend, 'making the likelihood scan plot'
            _plot_name = _poi_name+'_plc_scan_exost.'+self._plot_format
            c1 = TCanvas("c1", "c1", 600, 600)
            plcPlot = RooStats.LikelihoodIntervalPlot(plcInt)
            gROOT.SetStyle("Plain")
            plcPlot.Draw()
            c1.SaveAs(_plot_name)
            
        return (lower_limit, upper_limit)
