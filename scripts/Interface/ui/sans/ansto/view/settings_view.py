import abc

class SettingsView(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def show_empty(self):
        pass
    