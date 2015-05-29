#!/usr/bin/env python

# globals
exost_vtag='00.01.14'
debug = 0
legend = '[ExoStCore]:'

# global banner sequence
def banner():
    print ''
    print '\033[1mExoSt v', exost_vtag, '-- Developed by Gena Kukartsev\033[0m'
    print '                    Copyright (C) 2010 Brown University'
    print '                    All rights reserved, please read http://brux.hep.brown.edu/~kukarzev/exost/license.txt'
    print ''

    #print 'uses'    
    print ''

banner()
from ROOT import RooFit
roofit = RooFit()
# end of the banner sequence

# module banner
print '\033[1m'+legend, 'loaded\033[0m'


from ROOT import TCanvas
from ROOT import RooMsgService
#import exostMcmc as mcmc


###################################################
#
# actions
#
def action(ws, action_name, settings = None):
    legend = '[Action]:'
    print legend, 'beginning the ', action_name, 'action'

    if action_name == 'bayes':

        if settings['status'] == 'fail':
            print legend, 'Bayesian calculator failed to configure'
            print legend, 'action stopped.'
            return

        legend = '[Bayesian calculator]:'

        #mconf = ws.obj('exostModelConfig')
        mconf = ws.obj(settings['model_config_name'])

        bCalc = ws.obj('exostBayes')

        # to suppress messages when pdf goes to zero
        RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
  
        if mconf.GetParametersOfInterest().getSize() == 1:

            bInt = bCalc.GetInterval()


            cl = bCalc.ConfidenceLevel()
            print legend,  str(cl)+'% CL central interval: [', bInt.LowerLimit(), ' - ', bInt.UpperLimit(), ']'
            print legend, 'or', str(cl+(1.0-cl)/2), '% CL limits'

            # make a posterior plot
            if settings['do_posterior_plot']:

                plot = bCalc.GetPosteriorPlot()
                c1 = TCanvas('c1','Bayesian Calculator Result');
                c1.cd(1)
                plot.Draw()
                c1.SaveAs('bayesian_example_plot.'+settings['plot_format'])
  
        else:
            print legend, 'Error: Bayesian Calc. only supports one parameter of interest'

    print legend, 'the', action_name, 'action is done'

