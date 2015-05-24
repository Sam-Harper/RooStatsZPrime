#!/usr/bin/env python
#
# Exost: CMS EXOtica Statistics Toolkit
#
# Module: Implements base class for actions
#
# Author: Gena Kukartsev
#
# September 2010
#
############################################

# define special exception class to be raised
# if a singleton class is attempted to be
# recreated via constructor
# just because I don't program enough
class singletonException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



#
# singleton base class for all actions
# serves as the action manager
#
class exostAction:
    theInstance = None

    @classmethod
    def instance(cls):
        if cls == exostAction:
            if exostAction.theInstance == None:
                exostAction.theInstance = exostAction()

        return exostAction.theInstance

    def __init__(self):
        if self.__class__ == exostAction:
            if exostAction.theInstance:
                raise singletonException('The exostAction class is a singleton. Use exostAction.instance() to get its instance')

            exostAction.theInstance = self

        self.name = 'exostAction'
        self.config_name = 'Base Action'
        self.legend = '[ExoSt action manager]:'

        if self.__class__ == exostAction:
            self.action_dict = {}
        else:
            #need to have this flag set to be able to execute
            self.configured = False
        
    def register(self):
        print self.legend, 'registering', self.__class__.__name__
        exostAction.instance().action_dict[self.get_name()]=self

    def action(self, ws = None):
        print self.legend, 'action'

    def get_legend(self):
        return self.legend

    def get_name(self):
        return self.name

    def get_config_name(self):
        return self.config_name

    def print_actions(self):
        print self.legend, 'Available actions:'
        for act in self.action_dict:
            print self.legend, '  ', self.action_dict[act].get_config_name()

    def action_names(self):
        return self.action_dict.keys()


    #########################################################
    #
    # check and report a config value
    #
    def check_value(self, items, section, key, verbose = False):

        #legend = self.legend + '['+section+']:'
        legend = self.legend

        # getting items as a dictionary
        _items = {}
        _items.update(items[section])

        if key in _items.keys():
            _value = _items[key]
            print legend, key+' =', _value
            return _value
        else:
            if verbose:
                print legend, 'Error:', key, 'is not specified'
            return None



    #########################################################
    #
    # Configure an action
    #
    #  - 'items' is a dictiionary of lists of pairs
    #     from the config file as in exostConfig.items
    #
    def configure(self, ws, items, section, debug = 0, theAction = None):

        if self.__class__ != exostAction:
            print self.legend, 'configuration for this action is not implemented'
            return
        
        print self.legend, 'searching for action', theAction

        if theAction in self.action_dict:
            print self.legend, 'action', theAction, 'is found, configuring...'

            self.action_dict[theAction].configure(ws, items, section, debug)

        else:
            print self.legend, 'action', theAction, \
                  'is not registered, cannot configure'



    def execute(self, ws, debug = 0, theAction = None):

        if self.__class__ != exostAction:
            print self.legend, 'execute() method for this action is not implemented'
            return
        
        print self.legend, 'searching for action', theAction

        if theAction in self.action_dict:
            print self.legend, 'action', theAction, 'is found, executing...'

            if self.action_dict[theAction].configured:
                self.action_dict[theAction].execute(ws, debug)
            else:
                print self.legend, 'action', theAction, 'was not configured, cannot execute'

        else:
            print self.legend, 'action', theAction, \
                  'is not registered, cannot execute'
            

    # check if a generic object is in the workspace
    # FIXME: THIS IS CURRENTLY BROKEN! MADE to always return True
    def check_generic_object(self, ws, name, debug = 0):
        return True

        # workaround the crash when a nonexistent object is accessed
        gen_list = ws.allGenericObjects()
        _list = []
        print 'DEBUG', gen_list.size()
        #print 'DEBUG', getattr(gen_list, '')()
        cpp_iter = gen_list.begin()
        for ind in range(gen_list.size()):
            _list.append( gen_list.front().GetTitle() )
            print 'DEBUG', gen_list.size()
            print 'DEBUG front', gen_list.front().GetTitle()
            print 'DEBUG back', gen_list.back().GetTitle()
            #print 'DEBUG', gen_list.front().GetTitle()
            

        if name in _list:
            return True
        else:
            return False
        
