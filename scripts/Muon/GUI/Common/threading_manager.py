from __future__ import (absolute_import, division, print_function)

from qtpy import QtCore, QtWidgets
from qtpy.QtCore import Signal

import traceback
import sys


def split_list_into_n_parts(lst, n):
    return [lst[i::n] for i in xrange(n)]


def column(matrix, i):
    return [row[i] for row in matrix]


def split_kwarg_list(kwarg_list, n):
    split = lambda x: split_list_into_n_parts(x, n)
    chunks = list(map(split, kwarg_list.values()))
    return [dict(zip(kwarg_list.keys(), column(chunks, i))) for i in range(n)]


def threading_decorator(fn):
    def wrapper(progress_callback, **kwargs):
        results = []
        failed_results = []
        num_evals = len(kwargs.values()[0])
        input_list = [{key: value[i] for key, value in kwargs.items()} for i in range(num_evals)]
        for inputs in input_list:
            try:
                result = fn(**inputs)
            except Exception as e:
                failed_results += [e]
                continue
            results += [result]
            progress_callback.emit(1 / num_evals)
        return results, failed_results

    return wrapper


class WorkerSignals(QtCore.QObject):
    """
    Defines the signals available from a running worker thread.
    """
    started = Signal()
    finished = Signal()
    error = Signal(dict)
    result = Signal(tuple)
    progress = Signal(int)
    cancelled = Signal()


class Worker(QtCore.QThread):
    """
    Worker thread. Inherits from QThread to handle worker thread setup, signals and wrap-up. If we use QRunnable
    instead we would not be able to cancel the thread.

    fn : A function to evaluate
    kwargs : list of keyword arguments to pass to fn

    returns : A list whose first element is the kwargs dict and second argument the result
    """

    def __init__(self, fn, **kwargs):
        super(Worker, self).__init__()
        self.fn = threading_decorator(fn)
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        # Add callback so fn can update on its progress
        self.kwargs['progress_callback'] = self.signals.progress

    def run(self):
        self.signals.started.emit()
        try:
            result, fails = self.fn(**self.kwargs)
            if 'progress_callback' in self.kwargs: del self.kwargs['progress_callback']
            for fail in fails:
                exctype, value = type(fail), next(iter(fail.args), "")
                self.signals.error.emit({"inputs": self.kwargs,
                                         "excetype": exctype,
                                         "value": value,
                                         "traceback": traceback.format_exc()})
        except:
            # traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            if 'progress_callback' in self.kwargs: del self.kwargs['progress_callback']
            self.signals.error.emit({"inputs": self.kwargs,
                                     "excetype": exctype,
                                     "value": value,
                                     "traceback": traceback.format_exc()})
        else:
            if 'progress_callback' in self.kwargs: del self.kwargs['progress_callback']
            self.signals.result.emit((self.kwargs, result))
        finally:
            self.signals.finished.emit()

    def cancel(self):
        self.signals.cancelled.emit()
        # must use terminate (not quit/exit) as there is no event loop for the thread
        self.terminate()


class WorkerManager(QtWidgets.QWidget):
    """
    The WorkerManager class handles multi-threading of a function, by splitting an arbitrary number of list arguments
    equally across a specified number of threads.

    fn : A function to be evaluated, it must take all arguments as lists.
    arg_list : A list of (a single) input parameter to be given to fn
    num_threads : The max number of threads to be executed across the arg_list


    A list of callbacks can be passed into the constructor to respond to events (all threads finishing,
    threads updating on progress, threads catching exceptions from the function evaluation)
    """

    finished = Signal()  # used for unit tests
    cancelled = Signal()  # used for unit tests

    def __init__(self, fn, num_threads,
                 callback_on_threads_complete=lambda: 0,
                 callback_on_progress_update=lambda x: 0,
                 callback_on_thread_exception=lambda x: 0,
                 **kwarg_list):
        super(WorkerManager, self).__init__()

        # Callback called when all threads have executed
        self.thread_complete_callback = callback_on_threads_complete
        # Callback called when a single thread updates its progress
        self.progress_callback = callback_on_progress_update
        # Callback called when an error is thrown from a thread
        self.error_callback = callback_on_thread_exception

        self.fn = fn
        self.kwarg_list = kwarg_list

        # maintain a list of running threads to allow for cancelling
        self._threads = []

        self._results = {key: [] for key in self.kwarg_list.keys() + ['results']}
        self._failed_results = {key: [] for key in self.kwarg_list.keys()}

        self._progress = 0.0

        self._num_threads = num_threads  # total number of threads started
        self._thread_count = 0  # current threads running

        self.mutex = QtCore.QMutex()

    def _cancel_threads(self):
        for thread in self._threads:
            if thread.isRunning():
                thread.cancel()

    def _clear_threads(self):
        self._cancel_threads()
        self._threads = []
        self._thread_count = 0

    def _clear_results(self):
        self._progress = 0.0
        self._results = {key: [] for key in self.kwarg_list.keys() + ['results']}
        self._failed_results = {key: [] for key in self.kwarg_list.keys()}

    def cancel(self):
        self._clear_threads()
        self.cancelled.emit()

    def clear(self):
        self._clear_threads()
        self._clear_results()

    def is_running(self):
        return len(self._threads) > 0

    def start(self):
        if not self.is_running():
            self._clear_results()
            self._thread_count = self._num_threads
            self.on_thread_progress(0.0)

            # split each argument into num_threads equally sized lists
            thread_list = split_kwarg_list(self.kwarg_list, self._num_threads)
            for i, arg in enumerate(thread_list):
                worker = Worker(self.fn, **arg)
                self.connect_worker_signals(worker)
                print("\tStarting thread " + str(i + 1) + " with arg : ", arg)
                self._threads += [worker]
                worker.start()
        else:
            raise RuntimeError("Cannot start threads, "
                               "WorkerManager has active threads (call cancel() or clear() first)")

    def connect_worker_signals(self, worker):
        worker.signals.started.connect(self.on_thread_start)
        worker.signals.result.connect(self.on_thread_result)
        worker.signals.finished.connect(self.on_thread_complete)
        worker.signals.progress.connect(self.on_thread_progress)
        worker.signals.error.connect(self.on_thread_exception)
        worker.signals.cancelled.connect(self.on_thread_cancelled)

    def on_thread_start(self):
        pass

    def on_thread_cancelled(self):
        pass

    def on_thread_exception(self, kwargs):
        print("THREAD EXCEPTION", kwargs)
        self.error_callback(**kwargs)

        self.mutex.lock()
        kwargs = kwargs["inputs"]
        for key, item in kwargs.items():
            self._failed_results[key] += item
        self.mutex.unlock()

    def on_thread_result(self, s):
        self.mutex.lock()
        kwargs, results = s
        for key, item in kwargs.items():
            self._results[key] += item
        self._results['results'] += results
        self.mutex.unlock()

    def on_thread_complete(self):
        print("\tTHREAD " + str(self._thread_count) + " COMPLETE!\n")
        self.mutex.lock()
        self._thread_count -= 1
        self.mutex.unlock()
        if self._thread_count == 0:
            self._cancel_threads()
            # in case of exceptions in threads, set the progress to 100%
            if self._progress < 100.0:
                self._update_progress(100.0 - self._progress)
            self.thread_complete_callback()
            self.finished.emit()

    def _update_progress(self, percentage):
        self._progress += percentage
        self.progress_callback(self._progress)

    def on_thread_progress(self, percentage):
        self.mutex.lock()
        self._update_progress(percentage * 1 / self._num_threads)
        self.mutex.unlock()

    @property
    def results(self):
        return self._results
