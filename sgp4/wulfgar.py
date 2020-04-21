"""Support test functions, in a module small enough to carry inline."""

from unittest import TestCase
__unittest = 1  # Tell unittest not to include run() in test tracebacks.

def add_test_functions(loader, tests, module_name):
    """Run our main documentation as a test and test functions in this file."""

    module = __import__(module_name)

    def wrap_test_function(test):
        def run(self):
            return test()
        return run

    test_functions = [getattr(module, name) for name in dir(module)
                      if name.startswith('test_')]

    class TestFunctions(TestCase): pass

    for f in test_functions:
        setattr(TestFunctions, f.__name__, wrap_test_function(f))

    tests.addTest(loader.loadTestsFromTestCase(TestFunctions))
