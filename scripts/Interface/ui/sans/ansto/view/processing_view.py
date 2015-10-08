import abc

class ProcessingView(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def show_empty(self):
        pass