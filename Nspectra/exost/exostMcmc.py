#!/usr/bin/env python
#
# Exost: CMS EXOtica Statistics Toolkit
#
# Module: Markov Chain Monte Carlo
#
# Author: Gena Kukartsev
#
# September 2010
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

legend = '[MCMC]:'

# module banner
print '\033[1m'+legend, 'Markov chain Monte Carlo module loaded\033[0m'


class exostMcmc(exostAction):

    def __init__(self):

        # call the base class default constructor
        exostAction.__init__(self)
        
        self.legend = '[MCMC]:'
        self.name = 'mcmc'
        self.config_name = 'MCMC Calculator'

        self.register()
        
    #########################################################
    #
    # Configure the action
    #
    #  - 'items' is a dictionary of lists of pairs
    #     from the config file as in exostConfig.items
    #
    def configure(self, ws, items, section, debug = 0):
        print self.legend, 'configuring Markov chain MC calculator...'

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

        # proposal helper cache size
        _ph_cache_size = self.check_value(items, section, 'ph_cache_size')
        if _ph_cache_size == None:
            print self.legend, 'using default Proposal Helper cache size: 100'
            _ph_cache_size = 100 

        # desired confidence level
        _conf_level = self.check_value(items, section, 'confidence_level')
        if _conf_level == None:
            print self.legend, 'no confidence level specified, setting to 0.95'
            _conf_level = 0.95
                
        # left side tail fraction (0 for upper limit)
        _left_side = self.check_value(items, section, 'left_side_tail_fraction')
        if _left_side == None:
            print self.legend, 'no left side tail fraction specified, setting to 0 (upper limit)'
            _left_side = 0.0
                
        # number of burn-in steps
        _burn_in = self.check_value(items, section, 'burn_in_steps')
        if _burn_in == None:
            print self.legend, 'number of burn-in steps set to default 500'
            _burn_in = 500
                
        # number of bins in the posterior pdf hist
        _nbins = self.check_value(items, section, 'n_bins')
        if _nbins == None:
            print self.legend, 'number of bins in the posterior PDF histogram set to default 50'
            _nbins = 50
                
        # number of iterations
        _n_iter = self.check_value(items, section, 'number_of_iterations')
        if _n_iter == None:
            print self.legend, 'number of iterations steps set to default 100000'
            _n_iter = 100000
                
        # whether to make the posterior plot or not
        self._posterior = self.check_value(items, section, 'make_posterior_plot')
        if self._posterior == None:
            print self.legend, 'posterior plot is not requested'
            self._posterior = False
        elif self._posterior == 'True':
            self._posterior = True
        else:
            print self.legend, 'invalid assignment make_posterior_plot =', self._posterior
            print self.legend, 'posterior plot will not be made'
            self._posterior = False
            
        # plot format
        if self._posterior:
            self._plot_format = self.check_value(items, section, 'plot_format')
            if self._plot_format == None:
                print self.legend, 'no plot format specified, setting to PNG'
                self._plot_format = 'png'
                

        # configuring now...


        # to suppress messages when pdf goes to zero
        if debug > 0:
            print self.legend, 'will now fit in order to build a proposal function'
        #else:
        #    RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

        _model = _mconf.GetPdf()

        _fit = _model.fitTo(_data, RooFit.Save())

        _ph = RooStats.ProposalHelper()
        _ph.SetVariables(_fit.floatParsFinal())
        _ph.SetCovMatrix(_fit.covarianceMatrix())
        _ph.SetUpdateProposalParameters(True)
        _ph.SetCacheSize(int(_ph_cache_size))
        _pf = _ph.GetProposalFunction()
        
        self._mc = RooStats.MCMCCalculator(_data,_mconf)
        self._mc.SetConfidenceLevel( float(_conf_level) )
        self._mc.SetProposalFunction(_pf)
        self._mc.SetNumBurnInSteps( int(_burn_in) )
        self._mc.SetNumBins( int(_nbins) )
        self._mc.SetNumIters( int(_n_iter) )
        self._mc.SetLeftSideTailFraction( float(_left_side) )

        print legend, 'Markov chain MC calculator', 'exostMcmcCalc', 'is configured'

        # import the calculator into the Workspace
        #getattr(ws, 'import')(self._mc, 'exostMcmcCalc')
        #print legend, 'Markov chain MC calculator', 'exostMcmcCalc', 'is added to workspace', ws.GetName()
        #print legend, 'done'

        ws.Print()

        # flag that the action is configured successfully
        self.configured = True

        return {'do_posterior_plot':self._posterior}

    def execute(self, ws, debug = 0):
        print self.legend, 'Markov chain MC calculation started...'

        # time the action
        t = TStopwatch()
        t.Start()

        # model config is guaranteed to be here
        # the action won't configure without it
        # and this code will never run unless valid model config is found
        mconf = ws.obj(self._model_config_name)
            
        _poi = mconf.GetParametersOfInterest()

        mcInt = self._mc.GetInterval()

        # stop watch and print how long it took to get the interval
        t.Print()

        # iterate over all parameters of interest and print out the intervals
        # (could there be more than one?)
        _iter = _poi.createIterator()
        while True:
            _poi_name = _iter().GetTitle()
            
            lower_limit = mcInt.LowerLimit( _poi[_poi_name] )
            upper_limit = mcInt.UpperLimit( _poi[_poi_name] )
            
            print self.legend, 'MCMC interval for', _poi_name, 'is ['+ \
                  str(lower_limit) + ', ' + \
                  str(upper_limit) + ']'

            if _iter.Next() == None:
                break
            

        if self._posterior:
            # draw posterior plot

            print self.legend, 'making the posterior plot'
            _plot_name = _poi_name+'_mcmc_posterior_exost.'+self._plot_format
            c1 = TCanvas("c1", "c1", 600, 600)
            mcPlot = RooStats.MCMCIntervalPlot(mcInt)
            gROOT.SetStyle("Plain")
            mcPlot.Draw()
            c1.SaveAs(_plot_name)
            
        return (lower_limit, upper_limit)
