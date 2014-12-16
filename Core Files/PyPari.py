# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:24:20 2013

@author: Benjamin Manifold
"""

import ctypes
pari = ctypes.CDLL("/usr/lib/libpari.so")

#Pari function return type declarations. This makes Ctypes happy. Avoids Segfaults.
from restypes import *


initialized = False
#Must initialize the Pari stack before we do anything else.
def pari_init(stackSize, maxPrime):
  global initialized
  if initialized:
    raise RuntimeError, "Pari has already been initialized."
  else:
    pari.pari_init(stackSize, maxPrime)
    initialized = True

def cgiv(Gen):
  return pari.cgiv(Gen)

# Returns current stack pointer, which also points to most recently created Gen.
def getAVMA():
  return ctypes.c_long.in_dll(pari, "avma").value

def setAVMA(arg): #arg should be a recorded stack pointer
  ctypes.c_long.in_dll(pari, "avma").value = arg

#A class Wrapper for Pari Gens. Hopefully we will hide this behind a shell implementation.
class Gen:
  #"Generic Pari Gen object"
  def __init__(self, arg):
    if isinstance(arg, int) or isinstance(arg, long):
      av = getAVMA() # record stack pointer location
      self.ref = ctypes.cast(pari.gclone(pari.stoi(arg)), #clonse result to heap
                 ctypes.POINTER(ctypes.c_long)) 
      setAVMA(av) # restore stack pointer location
    if isinstance(arg, float):
      av = getAVMA()
      self.ref = ctypes.cast(pari.gclone(pari.dbltor(ctypes.c_double(arg))),
                 ctypes.POINTER(ctypes.c_long))
      setAVMA(av) 
#Evaluates a string passed in as arg as it would in GP and makes the correct Pari Object.
    if isinstance(arg, str):
      av = getAVMA() 
      self.ref = ctypes.cast(pari.gclone(pari.gp_read_str(arg)), 
                 ctypes.POINTER(ctypes.c_long))
      setAVMA(av)
        
    if isinstance(arg, list):
      av = getAVMA()
      self.ref = ctypes.cast(pari.gclone(pari.listcreate(len(arg))),
                 ctypes.POINTER(ctypes.c_long))
      for element in range(0, len(arg)): 
        pari.listput(self.ref, Gen(arg[element]).ref, element + 1)
        # self.ref = ctypes.cast(pari.gclone(pari.gp_read_str(arg.__repr__())),
        #           ctypes.POINTER(ctypes.c_long))
        setAVMA(av)
        #Treat any ctypes pointer to a c_long as a gen in need of a class wrapper
    if isinstance(arg, ctypes.POINTER(ctypes.c_long)):
      av = getAVMA()
      self.ref = pari.gclone(arg)
      setAVMA(av)
     
  def __add__(self, arg):
    av = getAVMA()
    if arg.__class__.__name__ == 'Gen':
      result = Gen(ctypes.cast(pari.gclone(pari.gadd(self.ref, arg.ref)), 
                   ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result
    else:
      result = Gen(ctypes.cast(pari.gclone(pari.gadd(self.ref, 
               pari.gp_read_str(repr(arg)))), ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result            
    
  def __sub__(self, arg):
    av = getAVMA()
    if arg.__class__.__name__ == 'Gen':
      result = Gen(ctypes.cast(pari.gclone(pari.gsub(self.ref, arg.ref)),
               ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result   
    else:
      result = Gen(ctypes.cast(pari.gclone(pari.gsub(self.ref, 
               pari.gp_read_str(repr(arg)))), ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result

  def __mul__(self, arg):
    av = getAVMA()
    if arg.__class__.__name__ == 'Gen':
      result = Gen(ctypes.cast(pari.gclone(pari.gmul(self.ref, arg.ref)),
               ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result   
    else:
      result = Gen(ctypes.cast(pari.gclone(pari.gmul(self.ref, 
               pari.gp_read_str(repr(arg)))), ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result
    
  def __div__(self, arg):
    av = getAVMA()
    if  arg.__class__.__name__ == 'Gen':
      result = Gen(ctypes.cast(pari.gclone(pari.gdiv(self.ref, arg.ref)), 
               ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result   
    else:
      result = Gen(ctypes.cast(pari.gclone(pari.gdiv(self.ref, 
               pari.gp_read_str(repr(arg)))), ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result
       
  def __floordiv__(self, arg):
    av = getAVMA()
    if  arg.__class__.__name__ == 'Gen':
      result = Gen(ctypes.cast(pari.gclone(pari.gdivent(self.ref, arg.ref)), 
               ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result 
    else:
      result = Gen(ctypes.cast(pari.gclone(pari.gdivent(self.ref, 
               pari.gp_read_str(repr(arg)))), ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result

  def __mod__(self, arg):
    av = getAVMA()
    if arg.__class__.__name__ == 'Gen':
      result = Gen(ctypes.cast(pari.gclone(pari.gmod(self.ref, arg.ref)), 
               ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result
    else:
      result = Gen(ctypes.cast(pari.gclone(pari.gmod(self.ref, 
               pari.stoi(arg))), ctypes.POINTER(ctypes.c_long)))
      setAVMA(av)
      return result

  def __repr__(self):
    return ctypes.string_at(pari.GENtostr(self.ref))

  def __del__(self):
    pari.gunclone(self.ref)
  
  #Returns the Gen type of some pari Gen object
  def gen_type(self):
    return ctypes.string_at(PyPari.pari.GENtostr(PyPari.pari.type0(self.ref)))
    
#Takes an integer, returns the address of a c_long, representing a GEN, that is a T_INT. Can go anywhere Pari expects a Gen.
def t_int(integer):
#Type checking/truncating code goes here. We don't want reals or lists of something that looks like a c_long.    
  return pari.stoi(integer)
    
#Takes a python real, returns a t_real type Gen Pari object.
def t_real(real):
  return pari.dbltor(ctypes.c_double(real))

#Takes Gen class objects, performs modular exponentiation through repeated squares method               
def pow_mod(gen_m, gen_n, gen_p):
  av = getAVMA()
  reduced = pari.stoi(1)
  bit_array = eval(ctypes.string_at(pari.GENtostr(pari.binaire(gen_n.ref))))
  for bit_index in range(len(bit_array)):
    reduced = pari.gmod(pari.gmul(reduced, reduced), gen_p.ref)
    if bit_array[bit_index] == 1:
      reduced = pari.gmod(pari.gmul(reduced, gen_m.ref), gen_p.ref)
  return Gen(reduced)
  
  
  setAVMA(av)
  return 
##Pari convertion functions

#Gen must be a T_INT. Returns a python integer converted from a c_long.
def genToInt(Gen):
  return pari.itos(Gen)
    
#Gen must be a T_REAL. Returns a python real converted from a c_double.
def genToReal(Gen):
  return pari.rtodbl(Gen)    

def fibo(n): # n can be entered in as an integer. ctypes will convet to a long int. Needs more work if function expects a Gen
   global initialized
   if not initialized:
     print("Pari has not yet been initialized.\nInitializing with default parameters (500000,500000).")
     pyPari_init(500000, 500000)
   return pari.fibo(n)




## This should convert an arbitrary python list pyarr into a C array. Or at least a list of integers.
#   arr = (ctypes.c_int * len(pyarr))(*pyarr)

from functionMap import *