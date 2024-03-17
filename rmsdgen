import numpy as np
from scipy.spatial.transform import Rotation as R
import glob
import matplotlib.pyplot as plt

def load_vectors(file_path):
    with open(file_path, 'r') as file:
        file.readline()
        file.readline()
    # Read vectors representing the atoms
        coordinates = []
        for line in file:
            data = line.strip().split()
            if len(data) == 4:
                x, y, z = float(data[1]), float(data[2]), float(data[3])
                coordinates.append([x, y, z])
    return np.array(coordinates)


def calculate_rmsd(vector_set1, vector_set2, rotation):
    squared_diff = np.sum((vector_set1 - rotation.apply(vector_set2)) ** 2, axis=1)
    mean_squared_diff = np.mean(squared_diff)
    rmsd = np.sqrt(mean_squared_diff)
    return rmsd

def normalize_vectors(vectors):
    magnitudes = np.linalg.norm(vectors, axis=1)
    normalized_vectors = vectors / magnitudes[:, np.newaxis]
    return normalized_vectors

def process_structures(isomer_list, destination_list, iteration, threshold):
    for i in range(iteration, len(isomer_list)):
        vectors_a = load_vectors(isomer_list[iteration-1])
        vectors_b = load_vectors(isomer_list[i])
        rot, rssd = R.align_vectors(vectors_a, vectors_b, return_sensitivity=False)
        normalized_vectors_a = normalize_vectors(vectors_a)
        normalized_vectors_b = normalize_vectors(vectors_b)
        norm_rmsd = calculate_rmsd(normalized_vectors_a, normalized_vectors_b, rot)
        if norm_rmsd > threshold:
            destination_list.append(isomer_list[i])

xyz_files = sorted(glob.glob('*.xyz'))

num_files = len(xyz_files)
rmsd_matrix = np.zeros((num_files, num_files))

for i in range(num_files):
    for j in range(i+1, num_files):
        file_a = xyz_files[i]
        file_b = xyz_files[j]

        vectors_a = load_vectors(file_a)
        vectors_b = load_vectors(file_b)

        rot, rssd = R.align_vectors(vectors_a, vectors_b, return_sensitivity=False)

        normalized_vectors_a = normalize_vectors(vectors_a)
        normalized_vectors_b = normalize_vectors(vectors_b)

        norm_rmsd = calculate_rmsd(normalized_vectors_a, normalized_vectors_b, rot)
        rmsd_matrix[i, j] = rmsd_matrix[j, i] = norm_rmsd

unique_isomers_1 = [xyz_files[0]]

unique_isomers_2 = []
unique_isomers_3 = []
unique_isomers_4 = []
unique_isomers_5 = []
unique_isomers_6 = []
unique_isomers_7 = []
unique_isomers_8 = []
unique_isomers_9 = []

for i in range(1, num_files):
        file_a = xyz_files[0]
        file_b = xyz_files[i]

        vectors_a = load_vectors(file_a)
        vectors_b = load_vectors(file_b)

        rot, rssd = R.align_vectors(vectors_a, vectors_b, return_sensitivity=False)

        normalized_vectors_a = normalize_vectors(vectors_a)
        normalized_vectors_b = normalize_vectors(vectors_b)

        norm_rmsd = calculate_rmsd(normalized_vectors_a, normalized_vectors_b, rot)
        if norm_rmsd > 0.4:
            unique_isomers_1.append(xyz_files[i])

print(f'Unique isomers1: {unique_isomers_1}\n\n')

unique_isomers_2.extend(unique_isomers_1[:2])

process_structures(unique_isomers_1,unique_isomers_2,2,0.4)

print(f'\n\nnew unique isomers 2:\n\n{unique_isomers_2}\n')

unique_isomers_3.extend(unique_isomers_2[:3])

process_structures(unique_isomers_2,unique_isomers_3,3,0.4)

print(f'\n\nnew unique isomers 3:\n\n{unique_isomers_3}\n')

unique_isomers_4.extend(unique_isomers_3[:4])

process_structures(unique_isomers_3,unique_isomers_4,4,0.4)

print(f'\n\nnew unique isomers 4:\n\n{unique_isomers_4}\n')

unique_isomers_5.extend(unique_isomers_4[:5])

process_structures(unique_isomers_4,unique_isomers_5,5,0.4)

print(f'\n\nnew unique isomers 5:\n\n{unique_isomers_5}\n')

unique_isomers_6.extend(unique_isomers_5[:6])

process_structures(unique_isomers_5,unique_isomers_6,6,0.4)

print(f'\n\nnew unique isomers 6:\n\n{unique_isomers_6}\n')

unique_isomers_7.extend(unique_isomers_6[:7])

process_structures(unique_isomers_6,unique_isomers_7,7,0.4)

print(f'\n\nnew unique isomers 7:\n\n{unique_isomers_7}\n')

unique_isomers_8.extend(unique_isomers_7[:8])

process_structures(unique_isomers_7,unique_isomers_8,8,0.4)

print(f'\n\nnew unique isomers 8:\n\n{unique_isomers_8}\n')

unique_isomers_9.extend(unique_isomers_8[:9])

process_structures(unique_isomers_8,unique_isomers_9,9,0.4)

print(f'\n\nnew unique isomers 9:\n\n{unique_isomers_9}\n')

print(len(unique_isomers_1))
print(len(unique_isomers_2))
print(len(unique_isomers_3))
print(len(unique_isomers_4))
print(len(unique_isomers_5))
print(len(unique_isomers_6))
print(len(unique_isomers_7))
print(len(unique_isomers_8))
print(len(unique_isomers_9))
