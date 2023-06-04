import numpy as np
from tqdm import tqdm



def getDist(mat1,mat2):
    n, m = np.meshgrid(mat1,mat2)
    dist = m - n
    dist -= np.rint(dist)
    return dist


def D(frac1,frac2,M):
    
    da, db, dc = (
        getDist(frac1[:,0],frac2[:,0]),
        getDist(frac1[:,1],frac2[:,1]),
        getDist(frac1[:,2],frac2[:,2]),
    )  # Fracrtional difference matrix

    # ---------cartesian differences------------

    DX = M[0][0] * da + M[1][0] * db + M[2][0] * dc
    DY = M[0][1] * da + M[1][1] * db + M[2][1] * dc
    DZ = M[0][2] * da + M[1][2] * db + M[2][2] * dc

    # -----------distance matrix--------------

    D = np.sqrt(np.square(DX) + np.square(DY) + np.square(DZ))
    
    return D


def get_molecules_ids(structure,rcut_OH=1.1): 
    
    molecule = {}
    H_dists = []

    count = 0
    
    print("getting molecules....")
    
    for site in tqdm(structure.sites):
        if site.specie.symbol == "O":
            b = []
            b.append(site.frac_coords.tolist()) 
            nbrs = structure.get_neighbors(site,rcut_OH)
            
            Hd = [] 
            
            if len(nbrs) !=2:
                print("OH cutoffs not working")
                break
            for s in nbrs:
                if  s.specie.symbol == "H":
                    b.append(s.frac_coords.tolist())
                    Hd.append(s[1])
                    
            H_dists.append(Hd)      
            molecule[count] = b
            count+=1
    return molecule,np.array(H_dists)



def get_adjacency(
                 structure,
                 rcut_OH=1.1,
                 rcut_OO=3.5,
                 phi = 30,
                 consider_Hbond=False 
                 
                    ):
    molecules,Hd = get_molecules_ids(structure)
    O_mat = np.array([item[0] for item in molecules.values()])
    M =  structure.lattice.matrix
    
    #-----------oridinary A calculation---------------
    
    A =  D(O_mat,O_mat,M)
    A[A<rcut_OO] = 1
    A[A>rcut_OO] = 0
    
    if not consider_Hbond:
        return A
    
    indices = np.tril_indices(A.shape[0], k=-1)
    connection = np.where(A[indices]==1)[0]
    connection_indices = np.array(indices)[:,connection].T
    
    
    #-----finding H bond network----------------
    

    not_connected = []
    
    print("getting H bonds....")

    for index in tqdm(connection_indices):
        mol1 = np.array(molecules[index[0]])
        mol2 = np.array(molecules[index[1]])
        dist = D(mol1,mol2,M)

        O1_O2 = dist[0,0]
        (O1H1,O1H2) = Hd[index[0],:]
        (O2H1,O2H2) = Hd[index[1],:]
        (O1H1_2,O1H2_2) = (dist[1,0],dist[2,0])
        (O2H1_1,O2H2_1) = (dist[0,1],dist[0,2])
        theta_1 = np.arccos((O1_O2**2+O1H1**2-O2H1_1**2)/(2*O1_O2*O1H1))*(180/np.pi)
        theta_2 = np.arccos((O1_O2**2+O1H2**2-O2H2_1**2)/(2*O1_O2*O1H2))*(180/np.pi)
        theta_3 = np.arccos((O1_O2**2+O2H1**2-O1H1_2**2)/(2*O1_O2*O2H1))*(180/np.pi)
        theta_4 = np.arccos((O1_O2**2+O2H2**2-O1H2_2**2)/(2*O1_O2*O2H2))*(180/np.pi)

        angles = [theta_1,theta_2,theta_3,theta_4]

        if all(i > phi for i in angles ):
            not_connected.append(index)
            
            
    # adjacancy matrix with H bond connection

    not_connected = np.array(not_connected).T

    i = not_connected[0,:].tolist()+not_connected[1,:].tolist()
    j = not_connected[1,:].tolist()+not_connected[0,:].tolist()

    A[i,j] = 0

    return A


