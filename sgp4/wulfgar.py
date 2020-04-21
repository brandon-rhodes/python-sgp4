"""Support test functions, in a module small enough to carry inline."""

from importlib import import_module
from unittest import TestCase
__unittest = 1  # Tell unittest not to include run() in test tracebacks.

def add_test_functions(loader, tests, module_name):
    """Run our main documentation as a test and test functions in this file."""

    def wrap_test_function(test):
        def run(self):
            return test()
        return run

    module = import_module(module_name)
    test_methods = dict((name, wrap_test_function(getattr(module, name)))
                        for name in dir(module) if name.startswith('test_'))
    TestFunctions = type('TestFunctions', (TestCase,), test_methods)
    TestFunctions.__module__ = module_name
    tests.addTest(loader.loadTestsFromTestCase(TestFunctions))
