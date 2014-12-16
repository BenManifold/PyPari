#Pari public function return declarations from paridecl.h
import copy
filelines = open("paridecl.h", "r").read().splitlines()
restypelist = []
restypelist.append('import ctypes')
restypelist.append('pari = ctypes.CDLL("/usr/lib/libpari.so")')

for line in filelines:
    if line.startswith('GEN') and (not line[8:].startswith('(')):
        restypelist.append('pari.' + line[8:line.index('(')] + '.restype = ctypes.POINTER(ctypes.c_long)')
    if line.startswith('INLINE GEN')  and (not line[8:].startswith('(')):
        restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.POINTER(ctypes.c_long)')
    #if line.startswith('ulong'):
    #    restypelist.append('pari.'+ line[8:line.index('(')] + '.restype = ctypes.c_ulong')
    #if line.startswith('void') and (not line[8:].startswith('(')):
    #    restypelist.append('pari.'+ line[8:line.index('(')] + '.restype = ctypes.c_void_p')
    #if line.startswith('long'):
    #    restypelist.append('pari.'+ line[8:line.index('(')] + '.restype = ctypes.c_long')
    #if line.startswith('int'):
    #    restypelist.append('pari.'+ line[8:line.index('(')] + '.restype = ctypes.c_int')
    #if line.startswith('double'):
    #    restypelist.append('pari.'+ line[8:line.index('(')] + '.restype = ctypes.c_double')
    #if line.startswith('char*'):
    #    restypelist.append('pari.'+ line[8:line.index('(')] + '.restype = ctypes.c_char_p')
    #if line.startswith('INLINE long'):
    #    restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.c_long')
    #if line.startswith('INLINE double'):
    #    restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.c_double')
    #if line.startswith('INLINE void'):
    #    restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.c_void_p')
    #if line.startswith('INLINE int'):
    #    restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.c_int')
    #if line.startswith('INLINE ulong'):
    #    restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.c_ulong')
    #if line.startswith('INLINE char*'):
    #    restypelist.append('pari.' + line[14:line.index('(')] + '.restype = ctypes.c_char_p')
    
#remove declarations with sensitive python syntax, or from hash.c
#clean_restypelist = []
#clean_restypelist[:] = [line for line in restypelist if ((line.count("hash") == 0) and
#			line.count(".print.")==0)]

declFile = open('restypes.py', "w")

for line in restypelist:
    declFile.write(line + '\n')
declFile.close()


 


