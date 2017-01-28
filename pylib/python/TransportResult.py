import numpy as np


class TransportResult (object):
    def __init__ (self, num_energy=0):
        self.NE = 0
        self.E = np.zeros(0)
        self._set_num_energy (num_energy)
        self.TE_op = {}
        self.I_op = {}
        self.DOS_op = {}
        self.n_op = {}
        self.neq_op = {}

        self.operators = {
                'TE_op': self.TE_op,
                'I_op': self.I_op,
                'DOS_op': self.DOS_op,
                'n_op': self.n_op,
                'neq_op': self.neq_op,
        }

    def get_num_energy (self):
        return self.NE

    def get_energy (self):
        return self.E

    def get_energy (self, iE):
        return self.E[iE]

    def get_matrix (self, type, iblock, jblock, iE):
        return operators [type][(iblock, jblock)][iE]

    def get_block_ids (self, type):
        return operators [type].keys ()

    def is_close (self, other, tolerance = 1E-1):
        if self.NE != other.NE:
            return False, "Number of energy point are not equal"
        for iE in range (self.NE):
            diff = abs ((self.E[iE] - other.E[iE])/self.E[iE])
            if diff > tolerance:
                return False, "Energey point # " + str (iE) + " are not equal"

        for type in self.operators.keys ():
            #print ("DBG: " + str (type))
            if type not in other.operators.keys ():
                return False, "Operator " + str (type) + " does not exist."
            for ijblocks in self.operators[type].keys ():
                #print ("DBG: " + str (ijblocks))
                if ijblocks not in other.operators[type].keys ():
                    return False, "Operator " + str (type) + " does not have blocks " + str (ijblocks) + "."

                this_matrix_arr = self.operators[type][ijblocks]
                other_matrix_arr = other.operators[type][ijblocks]
                if len (this_matrix_arr) != len (other_matrix_arr):
                    return False, "Number of energy points for Operator " + str (type) + " with blocks " + str (ijblocks) + " are not equal."

                for iE in range (len (this_matrix_arr)):
                    this_matrix = this_matrix_arr [iE]
                    other_matrix = other_matrix_arr [iE]
                    error = np.abs (this_matrix - other_matrix)
                    max_error = np.max (error)
                    if max_error > tolerance:
                        return False, "Error exeeds tolerance for Operator: " + str (type) + " blocks: " + str (ijblocks) + " energy # " + str (iE) + "."
                    #print ("DBG:" + str (error))

        return True
                

    def __eq__ (self, other):
        return self.is_close (other)
    
    def _set_num_energy (self, NE):
        self.NE = NE
        self.E = np.zeros (self.NE, dtype=float)

    def _set_energy (self, iE, E):
        self.E[iE] = E

    def _append_matrix (self, type, iblock, jblock, matrix):
        if type not in self.operators.keys ():
            raise RuntimeError ("Unknown matrix operator " + type)

        operator = self.operators [type]
        if not (iblock, jblock) in operator.keys ():
            operator [(iblock, jblock)] = []
        operator [(iblock, jblock)].append (matrix)

    @classmethod
    def read (cls, fileName): 
        """
        Imports output of transport calculation into numpy data structures.
        """
        out =  cls ()
    
        with open(fileName) as fid:
            for line in fid:
                line = line.strip()
    
                if line == 'ENERGY':
                    cls._scan_energy (out, fid)
                elif line == 'TRANSMISSION':
                    cls._scan_matrix (out, fid, 'TE_op')
                elif line == 'CURRENT':
                    cls._scan_matrix (out, fid, 'I_op')
                elif line == 'DOS':
                    cls._scan_matrix (out, fid, 'DOS_op')
                elif line == 'n':
                    cls._scan_matrix (out, fid, 'n_op')
                elif line == 'neq':
                    cls._scan_matrix (out, fid, 'neq_op')
    
        return out


    @classmethod
    def _scan_matrix (cls, out, fid, type):
        line = fid.readline ()
        NE = int (line)
    
        line = fid.readline ().strip ()
        blocks = line.split ()
        iblock = int (blocks[0])
        jblock = int (blocks [1]) 
    
        line = fid.readline ().strip()
    
        matrix_size = int (line)
        for iE in range (NE):
            matrix = np.zeros ((matrix_size, matrix_size), dtype=complex)
            row = 0
            while (row < matrix_size):
                line = fid.readline ().strip ()
                if len (line) == 0:
                    continue
                col_data = line.split ()
                if len (col_data) != matrix_size:
                    raise RuntimeError ("WW> Parse error at the line containing: '" + line + "': matrix does not have " + str (matrix_size) + " of columns")
    
                for col in range (matrix_size):
                    data = col_data [col].strip ("()").split(",")
    
                    if len (data) != 2:
                        raise RuntimeError ("WW> Parse error at the line containing: '" + line + "': expected complex number, got " + col_data[col])
                    data = complex (float (data[0]), float (data[1]))
                    matrix [row,col] = data
                row += 1
            out._append_matrix (type, iblock, jblock, matrix)
    
    @classmethod
    def _scan_energy (cls, out, fid):
        line = fid.readline ()
    
        NE = int (line)
        out._set_num_energy (NE)
    
        for iE in range (NE):
            line = fid.readline ()
            EE = float (line)
            out._set_energy (iE, EE)

