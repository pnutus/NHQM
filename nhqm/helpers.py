from __future__ import division
import scipy as sp

def matrix_from_function(function, order, dtype=complex, 
                         hermitian=False, symmetric=False):
    matrix = sp.empty((order, order), dtype)
    for i in xrange(order):
        limit = (i + 1) if hermitian or symmetric else order
        for j in xrange(limit):
            matrix[i, j] = function(i, j)
            if hermitian:
                matrix[j, i] = sp.conj(matrix[i, j])
            elif symmetric:
                matrix[j, i] = matrix[i, j]
    return matrix
    
#
#   Decorators
#
    
import collections
import functools

class memoize(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)