import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def draw_A(x=0, y=0, s=1, color='#bf0087', ax=None):
    points_A = [[x, y], [0.2 + x, y], [0.275 + x, s * 0.3 + y],
                [0.525 + x, s * 0.3 + y], [0.475 + x, s * 0.5 + y],
                [0.325 + x, s * 0.5 + y], [0.4 + x, s * 0.8 + y],
                [0.6 + x, y], [0.8 + x, y],
                [0.55 + x, s + y], [0.25 + x, s + y]]
    polygon_A = Polygon(points_A, color=color)

    if ax == None:
        ax = plt.gca()
    ax.add_patch(polygon_A)

def draw_C(x=0, y=0, s=1, color='#e5b600', ax=None):
    theta_outer = np.linspace(np.radians(53), np.radians(307), 50)
    theta_inner = np.linspace(np.radians(-53), np.radians(-307), 50)
    points_C_outer = [[0.5*np.cos(theta) + 0.5 + x,
                       s*0.5*np.sin(theta) + y + s * 0.5] for theta in theta_outer]
    points_C_inner = [[0.3*np.cos(theta) + 0.5 + x,
                       s*0.3*np.sin(theta) + y + s * 0.5] for theta in theta_inner]
    polygon_C = Polygon(points_C_outer + points_C_inner, color=color)

    if ax == None:
        ax = plt.gca()
    ax.add_patch(polygon_C)

def draw_G(x=0, y=0, s=1, color='#40ac00', ax=None):
    theta_outer = np.linspace(np.radians(53), np.radians(307), 50)
    theta_inner = np.linspace(np.radians(-70), np.radians(-307), 50)
    points_G_outer = [[0.5*np.cos(theta) + 0.5 + x,
                       s*0.5*np.sin(theta) + y + s * 0.5] for theta in theta_outer]
    points_G_inner = [[0.3*np.cos(theta) + 0.5 + x,
                       s*0.3*np.sin(theta) + y + s * 0.5] for theta in theta_inner]
    points_G_stem = [[points_G_outer[-1][0], s * 0.5 + y], [points_G_inner[0][0], s * 0.5 + y]]
    polygon_G = Polygon(points_G_outer + points_G_stem + points_G_inner, color=color)

    if ax == None:
        ax = plt.gca()
    ax.add_patch(polygon_G)

def draw_T(x=0, y=0, s=1, color='#5233ea', ax=None):
    points_T = [[x, s * 0.8 + y], [0.3 + x, s * 0.8 + y], [0.3 + x, y],
                [0.5 + x, y], [0.5 + x, s * 0.8 + y], [0.8 + x, s * 0.8 + y],
                [0.8 + x, s + y], [x, s + y]]
    polygon_T = Polygon(points_T, color=color)

    if ax == None:
        ax = plt.gca()
    ax.add_patch(polygon_T)

def prob_array(matrix):
    '''Returns an array of probabilities associated with each base at each position.
    Should be a vertical matrix (i.e. 4 columns of length L, where L is the length of the
    binding site.)'''
    probs = []
    for row in matrix:
        probs.append([np.exp(-val) / sum(np.exp(-row)) for val in row])
    return probs

def info_array(probs):
    '''Returns an array of information content for each sequence position.'''
    info_list = []
    for row in probs:
        info_list.append(2 + sum(val * np.log2(val) for val in row))
    return info_list

def seq_logo(matrix, ax=None, PWM=False, info_scaling=True, colors=['#bf0087', '#e5b600', '#40ac00', '#5233ea']):
    # If the input is an energy-weight matrix, create a matrix of probabilities
    if PWM==False:
        matrix = np.array(prob_array(matrix))

    if info_scaling==True:
        array = []
        for i, val in enumerate(info_array(matrix)):
            array.append(val * matrix[i])
    else:
        array = np.copy(matrix)

    # Map matrix values to base identity
    bases = ['A', 'C', 'G', 'T']
    bases_dict = {}
    for i, row in enumerate(array):
        bases_dict[i] = dict(zip(row, bases))

    # Sort matrix values in each row from smallest to largest
    sorted_array = []
    for row in array:
        sorted_array.append(sorted(row))

    for i, row in enumerate(sorted_array):
        a = 0
        for val in row:
            if bases_dict[i][val] == 'A': draw_A(x = i + 0.6, y = a, s = val, color=colors[0], ax=ax)
            elif bases_dict[i][val] == 'C': draw_C(x = i + 0.6, y = a, s = val, color=colors[1], ax=ax)
            elif bases_dict[i][val] == 'G': draw_G(x = i + 0.6, y = a, s = val, color=colors[2], ax=ax)
            else: draw_T(x = i + 0.6, y = a, s = val, color=colors[3], ax=ax)
            a += val
