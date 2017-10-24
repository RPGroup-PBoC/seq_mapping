def seq_to_array(seq):
    '''
    Converts DNA sequence to array for use in other functions

    INPUTS
    ------
    seq: string: DNA sequence

    OUTPUTS
    -------
    Array of numerical values corresponding with sequence
    '''
    mat = []
    seq_dict = {'A':0, 'C':1, 'G':2, 'T':3}
    for base in seq:
        mat.append(seq_dict[base])
    return mat

def binding_energy(seq, mat, wt_energy):
    '''
    Returns the binding energy of an operator sequence

    INPUTS
    ------
    seq : string : sequence whose energy you would like to calculate
    mat : array : raw energy matrix, 4 columns of length L where L is the operator length.
    wt_energy : float : binding energy of wild-type operator sequence in k_BT

    OUTPUTS
    -------
    Binding energy predicted for seq by mat
    '''

    # Convert seq to array if string
    if type(seq) is str:
        seq = seq_to_array(seq)

    # Determine binding energy for given sequence and convert to kBT
    Onew_energy = 0
    for i, val in enumerate(seq):
        Onew_energy += mat[i][val]
    Onew_energy = Onew_energy + wt_energy

    return Onew_energy

def fix_wt(matrix, seq):
    '''
    Alters an energy matrix so that the binding energy of the wild-type
    sequence is fixed at 0.

    INPUTS
    ------
    matrix: array: 4 X N energy matrix
    seq: string: wild-type sequence used to generate mutant library

    OUTPUTS
    -------
    Fixed energy matrix
    '''
    seq = seq_to_array(seq)
    mat_fixed = []
    for i, row in enumerate(matrix):
        mat_fixed.append([val - row[seq[i]] for val in row])

    return mat_fixed

def find_mult(seq, seq_energy, mat, wt_energy):
    '''
    Finds scaling factor for a matrix given known binding energy for a sequence

    INPUTS
    ------
    seq: string: sequence of interest
    seq_energy: float: measured binding energy for seq
    mat: array: raw 4 X N energy matrix
    wt_energy: float: binding energy of wild-type operator sequence

    OUTPUTS
    -------
    Scaling factor
    '''
    return (seq_energy - wt_energy)/binding_energy(seq, mat, 0)

def two_point_binding_energy(seq, df_mat, wt_en, mult):
    '''
    Predicts binding energy from two-point energy matrix model

    INPUTS
    ------
    seq: string: sequence of interest
    df_mat: DataFrame: two-point energy matrix model formatted as a tidy
                       DataFrame with headers "pos_1", "pos_2", "base_1",
                       "base_2", and "value"
    wt_en: float: binding energy of wild-type sequence in k_BT
    mult: float: matrix scaling factor

    OUTPUTS
    -------
    Predicted binding energy for seq
    '''
    import numpy as np
    
    energy = 0
    for i in range(len(seq)):
        for j in df_mat.pos_2[(df_mat.pos_1==i)].unique():
            energy += df_mat.value[(df_mat.pos_1==i) &\
                                   (df_mat.base_1==seq[int(i)]) &\
                                   (df_mat.pos_2==j) &\
                                   (df_mat.base_2==seq[int(np.floor(j))])].unique()[0]
    return(energy * mult + wt_en)
