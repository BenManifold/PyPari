# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 15:12:35 2014

@author: ben
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May 25 00:17:25 2014

@author: ben
"""

#Pari public function return declarations from paridecl.h
import keyword
decl_lines = open("paridecl.h", "r").read().splitlines()
desc_functions = open("pari.desc", "r").read().split("Function: ")[1:]
split_desc_functions = []
specials = ('_','!','@','#','$','%','O(','[','-')
cnames = {}
help_dict = {}
param_dict = {}

output_lines = []
#output_lines.append('import ctypes')
#output_lines.append('pari = ctypes.CDLL("/usr/lib/libpari.so")')
output_lines.append('import PyPari')
output_lines.append('import ctypes')
output_lines.append('PyPari.pari_init(400000,400000)\n')
#Utility Function for filtering lines that are not function declarations
def is_clean(line):
  decl_line_indicators = ('GEN ','ulong ', 'void ','long',
  'int ','double ', 'char* ','INLINE GEN ', 'INLINE long ', 
  'INLINE double ', 'INLINE void ', 'INLINE int ', 'INLINE ulong ',
  'INLINE char* ')  
  if len(line.split(' ',1)) > 1:
    if not line.split(' ',1)[1].strip().startswith('('):
      if not (" print(" in line or 'hash' in line):
        if line.startswith(decl_line_indicators):
          return True
        else:
          return False
      else:
        return False
    else:
      return False
  else:
    return False

for funct in desc_functions:
  split_desc_functions.append(funct.splitlines())

reserved_pyIds = keyword.kwlist + dir(__builtins__)
clean_functs = [entry for entry in split_desc_functions if 
  (entry[0][0].isalpha() and entry[0].isalnum() and 
  (not entry[0] in reserved_pyIds))]
  
  
gp_params_Map = {}
gp_params_list = []  
for funct in clean_functs:
    
    funct_params_list = []
    for i in range(1, len(funct)):
        if funct[i].startswith("Description: "):
            for j in range(i, len(funct)):
                if funct[j].startswith(" "):
                    #Do things to things here
                    gp_params = funct[j].strip().split('(')[1].split(')')[0].split(",")
                    gp_params = [p.replace(" ", "") for p in gp_params]
                    c_call_prototype = funct[j].split('Doc:')[0].split('Variant')[0].split(')',1)[1].split(' ', 1)[1].strip()
                    funct_params_list.append([gp_params, c_call_prototype])
                    #print(funct[0] + ": " + c_call_prototype)                    
                    #print(gp_params)
                    #if c_call_prototype.startswith('$'):
                    #  print(funct[0] + ": " + c_call_prototype)  #For investigating problematic functions
                      #print(c_params)
                if funct[j].startswith("Doc") or funct[j].startswith("Variant"):
                    break
            gp_params_Map[funct[0]] = funct_params_list
            

#Building GP - Parilib function name dictionary
for funct in clean_functs:
  for line in funct:
    if line.startswith("C-Name: ") and not funct[0].startswith(specials):
      cnames[ funct[0].strip() ] = line.split(":")[1].strip()
    if line.startswith("Help: ") and not funct[0].startswith(specials):
      help_dict[ funct[0] ] = '\n'.join(funct)
      
#Inverting Function Name Dictionary
gpnames = dict([v,k] for k,v in cnames.items())

#Building Parameter-list dictionaries tied to c-names from paridecl.h      
cleanlines = [line for line in decl_lines if is_clean(line)]


  
  
  
for line in cleanlines:
  declaration = line.split(' ',1)[1].strip().split('(',1)
  c_lib_name = declaration[0]
  c_params = declaration[1].split(')',1)[0].split(',')
  param_dict[ c_lib_name ] = c_params

#Build rest of output file
for line in cleanlines:
  c_lib_name = line.split(' ',1)[1].strip().split('(',1)[0]
  c_restype = line.split(' ',1)[0].strip()
  params_amount = len(param_dict[ c_lib_name ])
  A = 65
  params = param_dict[ c_lib_name ]
  if c_lib_name in cnames.values() and gpnames[ c_lib_name ] not in reserved_pyIds:
    output_lines.append("def " + gpnames[ c_lib_name ] + "(*argv):")
    #Introduce valid docstring derived from help_dict
    output_lines.append("  '''")
    helplines = []
    for helpline in help_dict[ gpnames[ c_lib_name]].splitlines():
      helplines.append("  " + helpline)
      docstring = '\n'.join(helplines)
    output_lines.append("  " + docstring.strip())
    output_lines.append("  '''")
    #parse params list from dictionary, add '.ref' to appropriately indexed arg if gen is expected.
    output_lines.append("  c_params = []")
    for arg_num in range(len(params)):
      if 'GEN' in params[arg_num]:
        output_lines.append("  c_params.append(argv["+repr(arg_num)+"].ref)")
      else:
        output_lines.append("  c_params.append(argv["+repr(arg_num)+"])")
    output_lines.append("  c_arg_tuple = tuple(c_params)")
    if not line.strip().startswith('void') or line.strip().startswith('INLINE void'):
      #Resetting pari stack after returning Gen, slow and safe      
      if line.strip().startswith('GEN') or line.strip().startswith('INLINE GEN'):
        output_lines.append("  av = PyPari.getAVMA()")
        output_lines.append("  result = PyPari.Gen(ctypes.cast(PyPari.pari."
        + "gclone(PyPari.pari."+ c_lib_name + "(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))")
        output_lines.append("  PyPari.setAVMA(av)")
        output_lines.append("  return result")
      else:
        output_lines.append("  return PyPari.pari."+c_lib_name+"(*c_arg_tuple)")
    else:
      output_lines.append("  PyPari.pari."+c_lib_name+"(*c_arg_tuple)")
    output_lines.append("")
  
declFile = open('functionMap.py', "w")

for line in output_lines:
    declFile.write(line + '\n')
declFile.close()    




