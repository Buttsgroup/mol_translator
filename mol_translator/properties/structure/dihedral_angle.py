import numpy as np

def get_dihedral_angle(aemol):
    dihedral_angle = np.zeros((aemol.structure['size'], aemol.structure['size'], aemol.structure['size'], aemol.structure['size']))
    
    for x in range(aemol.structure['size']):
        for y in range(aemol.structure['size']):
            ab_vec = aemol.structure['xyz'][y] - aemol.structure['xyz'][x]
            for z in range(aemol.structure['size']):
                ac_vec = aemol.structure['xyz'][z] - aemol.structure['xyz'][y]
                ac_vec /= np.linalg.norm(ac_vec)
                for q in range(aemol.structure['size']):
                    if len(set([x, y, z, q])) != 4:
                        continue
    
                    cd_vec = aemol.structure['xyz'][q] - aemol.structure['xyz'][z]
                    
                    # v: projection of b0 onto plane perpendicular to b1
                    v = ab_vec - np.dot(ab_vec, ac_vec)*ac_vec
                    # w: projection of b2 onto plane perpendicular to b1
                    w = cd_vec - np.dot(cd_vec, ac_vec)*ac_vec
                    
                    t = np.dot(v, w)
                    b = np.dot(np.cross(ac_vec, v), w)
                    
                    angle = np.degrees(np.arctan2(t, b))
                    
                    dihedral_angle[x][y][z][q] = angle
                    
    return dihedral_angle
                    
    
    
    
    
    
    