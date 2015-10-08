import abc

class MainView(object):

    #__metaclass__ = abc.ABCMeta

    #@abc.abstractmethod
    def get_processing_view(self):
        pass

    #@abc.abstractmethod
    def select_processing_view(self):
        pass