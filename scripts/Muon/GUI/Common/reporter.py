from __future__ import (absolute_import, division, print_function)

from mantid.kernel import ConfigService
import os
import getpass

class Reporter(object):

    """
    A class for recording user interaction
    """

    def __init__(self,file_name):
        config = ConfigService.Instance()
        path = os.path.join(config.getAppDataDirectory(),"MuonRecovery")
        user = getpass.getuser()
        recover_path = os.path.join(path, user)
        self.makePath(recover_path)
        self.file_name = os.path.join(recover_path, file_name)
 
    def makePath(self,path):
        config = ConfigService.Instance()
        if os.path.exists(path):
            return
        os.mkdir(path)
        return

    def exists(self):
        return os.path.exists(self.file_name)

    def clear(self):
        store = open(self.file_name,"w")
        store.close()

    def report(self, command):
        store = open(self.file_name,"a")
        store.write(command)
        store.close()
