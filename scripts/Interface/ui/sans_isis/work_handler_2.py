from PyQt4.QtCore import pyqtSlot, QThreadPool
from abc import ABCMeta, abstractmethod
from six import with_metaclass
from worker import Worker


class WorkHandler2(object):
    class WorkListener(with_metaclass(ABCMeta, object)):
        def __init__(self):
            pass

        @abstractmethod
        def on_processing_finished(self, result):
            pass

        @abstractmethod
        def on_processing_error(self, error):
            pass

    def __init__(self):
        self.thread_pool = None
        self._listener = None
        self._worker = None

    @pyqtSlot()
    def on_finished(self):
        print('returned on finished')
        result = self._worker.result
        self._worker = None
        self._finished_callback(result)

    @pyqtSlot()
    def on_error(self, error):
        self._worker = None
        self._error_callback(error)

    def process(self, finished_callback, error_callback, func, *args, **kwargs):
        # Add the caller
        self._error_callback = error_callback
        self._finished_callback = finished_callback

        # Generate worker
        self._worker = Worker(func, *args, **kwargs)
        self._worker.signals.finished.connect(self.on_finished)
        self._worker.signals.error.connect(self.on_error)

        QThreadPool.globalInstance().start(self._worker)