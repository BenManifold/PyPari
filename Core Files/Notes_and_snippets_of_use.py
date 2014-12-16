# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 23:45:34 2014

@author: ben
"""

# Gp_param_to_c_param Logic Prototype
output_lines.append("gen_desc_types = ['int', 'real', 'mp', 'vec', 'list', 'smallvec', 'pol', 'var', 'gen', '?gen', '?int', '?var', 'genstr']")
output_lines.append("mp_list = [ 'mp', '?mp' ]")
output_lines.append("int_list = ['int','?int']")
output_lines.append("real_list = ['real','?real']")
output_lines.append("pol_list = ['pol','?pol']")
output_lines.append("vec_list = ['vec','?vec']")
output_lines.append("list_list = ['list','?list']")
output_lines.append("smallvec_list = ['smallvec','?smallvec']")

for gp_function_name in gp_params_Map:
  output_lines.append("def " + gp_function_name + "(*argv):")
  #adding valid python docstrings
  output_lines.append("  '''")
  helplines = []
  for helpline in help_dict[ gp_function_name ].splitlines():
    helplines.append("  " + helpline)
    docstring = '\n'.join(helplines)
    output_lines.append("  " + docstring.strip())
    output_lines.append("  '''")
    
  output_lines.append("  c_params = []")
  output_lines.append("  py_param_types = []")
  #collect type information at runtime to match up with correct gp param combo
  ## 
  output_lines.append("  for arg in argv:")
  output_lines.append("    if arg.__class__.__name__ == Gen:")
  output_lines.append("      py_param_types.append(arg.gen_type())")
  output_lines.append("    else:")
  output_lines.append("      py_param_types.append('small)") #no non-small, non-gen gp params as of pari 2.5.4, may need to add doubles
  output_lines.append("  param_pairs = " + repr(gp_params_Map[gp_function_name]))
  output_lines.append("  for pair in param_pairs:")
  output_lines.append("    if len(argv) > len(pair[0]):")
  output_lines.append("      continue")
  output_lines.append("    for gp_param_index in range(len(pair[0][:len(argv)])):") #We can't test more parameters than the user supplied!
  output_lines.append("      if (pair[0][gp_param_index] in gen_desc_types) and (arg.__class__.__name__ is not 'Gen'):")
  output_lines.append("        break")
  output_lines.append("      if (pair[0][gp_param_index] in  mp_list) and (arg.gen_type() not in ['t_INT', 't_REAL']):")
  output_lines.append("        break")
  output_lines.append("      if (pair[0][gp_param_index] in int_list) and (arg.gen_type() is not 't_INT'):")
  output_lines.append("        break")
  output_lines.append("      if pair[0][gp_param_index] in real_list and arg.gen_type() is not 't_REAL':")
  output_lines.append("        break")
  output_lines.append("      if pair[0][gp_param_index] in pol_list and arg.gen_type() is not 't_POL':")
  output_lines.append("        break")
  output_lines.append("      if pair[0][gp_param_index] in vec_list and arg.gen_type() is not 't_VEC':")
  output_lines.append("        break")
  output_lines.append("      if pair[0][gp_param_index] in list_list and arg.gen_type() is not 't_LIST':")
  output_lines.append("        break")
  output_lines.append("      if pair[0][gp_param_index] in smallvec_list, and arg.gen_type() is not 't_SMALLVEC':")
  output_lines.append("        break")
  output_lines.append("      ")
    
  # and so on
  #output_lines.append("        py_param_types.append('int')")
  #output_lines.append("      if arg.gen_type() == 't_REAL':")
  #output_lines.append("        py_param_types.append('real')")
  #output_lines.append("      py_param_types.append(arg.gen_type())")
  #output_lines.append("    else:")
  #assuming non-gen args to gp functions are just long ints, safe right now as per pari.desc
  #output_lines.append("      py_param_types.append('small')") 
  #output_lines.append("  param_pairs = " + repr(gp_params_Map[gp_function_name]))
  #output_lines.append("  for param in py_param_types:")
  #output_lines.append("    param.replace('Gen', 'gen')")
  #output_lines.append("    param.replace('t_INT', 'int')")
  #output_lines.append("    param.replace('t_REAL)")
  for pair in gp_params_Map[gp_function_name]:
    
    output_lines.append("  if len(argv) > " + repr(len(pair[0])) + ":")
    output_lines.append("    continue")
    cname = pair[1].split('(',1)[0]
    c_params = pair[1].split('(',1)[1].split(')',1)[0]
        
    for gp_param_index in range(len(pair[0])):
      
      if pair[0][gp_param_index] in ['gen', 'pol', 'vec', 'vecsmall', 'list', 'var', 'genstr']: ##any manner of gens described in pari.desc
        c_call_actual = c_call_actual.replace("$" + repr(gp_param_index+1), "argv[" + repr(gp_param_index) + "].ref")
        #output_lines.append("  c_params.append(argv[" + repr(gp_param_index) + "].ref")
      if pair[0][gp_param_index] in ['small','int', 'real', 'mp','str', 'bool', 'negbool', 'float']:    
        c_call_actual = c_call_actual.replace("$" + repr(gp_param_index+1), "argv[" + repr(gp_param_index) + "]")
        #output_lines.append("  c_params.append(argv[" + repr(gp_param_index) + "]")
      #Handling optional parameters.   
      if pair[0][gp_param_index] in ['?small','?int', '?real', '?mp','?str', '?bool', '?negbool', 'float']:
        c_call_actual = c_call_actual.replace()
        #output_lines.append("  try:")
        #output_lines.append("    c_params.append(argv[" + repr(gp_param_index) + "]")
        #output_lines.append("  except IndexError:")
        #output_lines.append("    continue")
      #Do some things to allow optional python parameter
      if pair[0][gp_param_index] in ['?gen', '?pol', '?vec', '?vecsmall', '?list', '?var', '?genstr']:
      #Do some things to allow optional python parameter respecting the need to pass the .ref attribute of the gen class to ctypes.
        #output_lines.append("  try:")
        #output_lines.append("    c_params.append(argv[" + repr(gp_param_index) + "].ref")
        #output_lines.append("  except IndexError:")
        output_lines.append("    continue")   
        
        
for pair in gp_params_Map['bnfinit']:
  c_call_actual = pair[1]
  for gp_param_index in range(len(pair[0])):
    #Gp_parameters in the description file not preceeded by '?' are mandatory. 
    if pair[0][gp_param_index] in ['gen', 'pol', 'vec', 'vecsmall', 'list', 'var', 'genstr']:
      c_call_actual = c_call_actual.replace("$" + repr(gp_param_index+1), "argv[" + repr(gp_param_index) + "].ref")
    if pair[0][gp_param_index] in ['small','int', 'real', 'mp','str', 'bool', 'negbool']:    
      c_call_actual = c_call_actual.replace("$" + repr(gp_param_index+1), "argv[" + repr(gp_param_index) + "]")
    #Handling optional parameters.   
    if pair[0][gp_param_index] in ['?small','?int', '?real', '?mp','?str', '?bool', '?negbool']:
      #Do some things to allow optional python parameter
    if pair[0][gp_param_index] in ['?gen', '?pol', '?vec', '?vecsmall', '?list', '?var', '?genstr']:
      #Do some things to allow optional python parameter respecting the need to pass the .ref attribute of the gen class to ctypes.
      

        
for function in gp_params_Map:
  for param_pair in gp_params_Map[function]:
    for gp_param_index in range(len(param_pair[0])):
      try:
        param_pair[1] = param_pair[1].replace("$" + repr(gp_param_index+1), "argv["+repr(gp_param_index) + "]")
      except IndexError:
        print param_pair
      
for function in gp_params_Map:
  for param_pair in gp_params_Map[function]:
    for gp_param_index in range(len(param_pair[0])):
      if pair[0][gp_param_index] in ['gen', 'pol', 'vec', 'vecsmall', 'list', 'var', 'genstr']:
        param_pair[1] = param_pair[1].replace("$" + repr(gp_param_index+1), "argv["+repr(gp_param_index) + "].ref")
        param_pair[1] = param_pair[1].replace("&", "ctypes.addressof(")
        param_pair[1] = param_pair[1].replace(repr(gp_param_index)+"].ref", repr(gp_param_index) + "].ref)")
      if pair[0][gp_param_index] in ['small','int', 'real', 'mp','str', 'bool', 'negbool']:
        param_pair[1] = param_pair[1].replace("$" + repr(gp_param_index+1), "argv["+repr(gp_param_index) + "]")
        param_pair[1] = param_pair[1].replace("&", "ctypes.addressof(")
        param_pair[1] = param_pair[1].replace(repr(gp_param_index)+"]", repr(gp_param_index) + "])")
      if pair[0][gp_param_index] in ['?gen', '?pol', '?vec', '?vecsmall', '?list', '?var', '?genstr']:
        param_pair[1] = param_pair[1].replace("$" + repr(gp_param_index+1), "argv["+repr(gp_param_index) + "].ref")
        param_pair[1] = param_pair[1].replace("&", "ctypes.addressof(")
        param_pair[1] = param_pair[1].replace(repr(gp_param_index)+"].ref", repr(gp_param_index) + "].ref)")  
      if pair[0][gp_param_index] in ['?small','?int', '?real', '?mp','?str', '?bool', '?negbool']:  
        param_pair[1] = param_pair[1].replace("$" + repr(gp_param_index+1), "argv["+repr(gp_param_index) + "]")
        param_pair[1] = param_pair[1].replace("&", "ctypes.addressof(")
        param_pair[1] = param_pair[1].replace(repr(gp_param_index)+"]", repr(gp_param_index) + "])")  
      

