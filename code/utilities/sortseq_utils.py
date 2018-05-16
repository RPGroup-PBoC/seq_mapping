import numpy as np

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

# Functions for calculating leakiness, saturation, and dynamic range

def leakiness(K_A, K_I, e_AI, R, Op):
    '''
    Computes the leakiness of a simple repression construct
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    R : float
        Number of repressors per cell
    Op : float
        Operator binding energy
    Returns
    -------
    leakiness
    '''
    return 1 / (1 + 1 / (1 + np.exp(-e_AI)) * R / 5E6 * np.exp(-Op))

def saturation(K_A, K_I, e_AI, R, Op):
    '''
    Computes the saturation of a simple repression construct
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    R : float
        Number of repressors per cell
    Op : float
        Operator binding energy
    Returns
    -------
    saturation
    '''
    return 1 / (1 + 1 / (1 + np.exp(-e_AI) * (K_A / K_I)**2) * R / 5E6 * np.exp(-Op))

def dynamic_range(K_A, K_I, e_AI, R, Op):
    '''
    Computes the dynamic range of a simple repression construct
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    R : float
        Number of repressors per cell
    Op : float
        Operator binding energy
    Returns
    -------
    dynamic range
    '''
    return saturation(K_A, K_I, e_AI, R, Op) - leakiness(K_A, K_I, e_AI, R, Op)

def EC50(K_A, K_I, e_AI, R, Op):
    '''
    Computes the concentration at which half of the repressors are in the active state
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    Concentration at which half of repressors are active (EC50)
    '''
    t = 1 + (R / 4.6E6) * np.exp(-Op) + (K_A / K_I)**2 * (2 * np.exp(-e_AI) + 1 + (R / 4.6E6) * np.exp(-Op))
    b = 2 * (1 + (R / 4.6E6) * np.exp(-Op)) + np.exp(-e_AI) + (K_A / K_I)**2 * np.exp(-e_AI)
    return K_A * ((K_A / K_I - 1)/(K_A / K_I - (t/b)**(1/2)) -1)

def pact(iptg, Ka, Ki, epsilon):
    '''returns probability that repressor is in active state'''
    return (1 + iptg/Ka)**2/((1 + iptg/Ka)**2 + np.exp(-epsilon) * (1 + iptg/Ki)**2)

def fold_change(IPTG, K_A, K_I, e_AI, R, Op):
    '''
    Computes fold-change for simple repression
    Parameters
    ----------
    IPTG : array-like
        Array of IPTG concentrations in uM
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    R : float
        Number of repressors per cell
    Op : float
        Operator binding energy
    Returns
    -------
    probability that repressor is active
    '''
    return 1 / (1 + R / 4.6E6 * pact(IPTG, K_A, K_I, e_AI) * np.exp(-Op))
