#!/usr/bin/env python

from optparse import OptionParser

legend = '[Command line parser]:'

# module banner
print '\033[1m'+legend, 'loaded\033[0m'

#############################################
#
#   command line argument parser
#
class exostCLParser:

    def __init__(self):
        add_help_option = "exost -a ACTION -c CONFIG_FILE"
    
        self.parser = OptionParser(add_help_option)
    
        self.parser.add_option("-a", "--action", dest="action", default=None,
                          help="Action to perform", metavar="ACTION")
    
        self.parser.add_option("-c", "--config-file", dest="config_file", default=None,
                          help="Configuration file", metavar="CONFIGFILE")
    
        self.parser.add_option("-d", "--debug-level", dest="debug", default=0,
                          help="Verbosity of the output", metavar="DEBUG_LEVEL")
    
        self.parser.add_option("--force-standalone", dest="standalone", default=False,
                         help="assume standalone ROOT installation",
                         action="store_true", metavar="STANDALONE")
        
        (self.options, args) = self.parser.parse_args()

        print legend, 'Config file:', self.options.config_file
        print legend, 'Action to perform:', self.options.action
        if self.options.debug>0:
            print legend, 'Output verbosity level set to', self.options.debug

    def get_options(self):
        return self.options
    
    
    def print_help(self):
        self.parser.print_help()
#
##### end of command line parser
