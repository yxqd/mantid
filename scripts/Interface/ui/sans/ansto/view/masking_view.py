import abc

class MaskingView(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def show_empty(self):
        pass

