
import find_existing_matrix
import He5
import HO_basis
import succesive_matrix_generation

matris_path = '/Users/Spill/Documents/kandidatarbete/matrisdata'
matris_format = '.txt'
matris_namn = 'He5_matris_'

n = 8
toggle = 0

def H_element_wrapper(n,nprim):
    #lement(n, n_prim, l, s, omega, V):
    return HO_basis.H_element(n,nprim,0,0.5,1,He5.V)

if toggle:
    exist, existingFilePath = succesive_matrix_generation.findExistingMatrix( n, matris_namn, matris_path, matris_format )
    if len( exist ) != 0:
        print 'existing matrix: '
        print exist
    else:
        print 'there is no applicable existing matrix'    

    hamilton = succesive_matrix_generation.generateSuccMatrixFromMatrix( n, H_element_wrapper, exist,matris_namn, matris_path, matris_format )
    print 'generated matrix: ' 
    print hamilton
else:
    hamilton = succesive_matrix_generation.generateSuccMatrix( n, H_element_wrapper, matris_namn, matris_path, matris_format, 1 )
    print 'generated matrix: ' 
    print hamilton
