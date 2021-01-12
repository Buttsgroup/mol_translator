
import numpy as np

def get_bond_angles(aemol):
    bond_angles = np.zeros((aemol.structure['size'], aemol.structure['size'], aemol.structure['size']))
    
    for i in range(aemol.structure['size']):
        for j in range(aemol.structure['size']):
            ba_vec = aemol.structure['xyz'][i] - aemol.structure['xyz'][j]
            
            for k in range(aemol.structure['size']):
                if len(set([i, j, k])) != 3 or bond_angles[i][j][k] != 0:
                    continue
                
                bc_vec = aemol.structure['xyz'][k] - aemol.structure['xyz'][j]
                
                cosine_angle = np.dot(ba_vec, bc_vec) / (np.linalg.norm(ba_vec) * np.linalg.norm(bc_vec))
                angle = np.degrees(np.arccos(cosine_angle))
                
                bond_angles[i][j][k] = angle
                bond_angles[k][j][i] = angle
                
    return bond_angles