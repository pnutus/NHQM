import He5
import HO_basis
import successive_matrix_generation as smg

matris_path = '/Users/Spill/Documents/kandidatarbete/matrisdata'
matris_format = '.txt'
matris_namn = 'He5_matris_'

n = 14
toggle = False

#not working

#element generating function
def H_element_wrapper(n,nprim):
    #lement(n, n_prim, l, s, omega, V):
    return HO_basis.H_element(n,nprim,0,0.5,1,He5.V)

#a quick switch to alternate between the two generation methods. generateSuccMatrix is the preferred 
#method as is calls findExistingMatrix internally and deletes uneeded matrix files
if toggle:
    exist, existing_file_path = smg.find_existing_matrix( matris_namn, matris_path, matris_format, n )
    if len( exist ) != 0:
        print 'existing matrix: '
        print exist
    else:
        print 'there is no applicable existing matrix'    

    hamilton = smg.generate_succ_matrix_from_matrix( n, H_element_wrapper, exist,matris_namn, matris_path, matris_format )
    print 'generated matrix: ' 
    print hamilton
else:
    hamilton = smg.generate_succ_matrix( n, H_element_wrapper, matris_namn, matris_path, matris_format, True )
    print 'generated matrix: ' 
    print hamilton
