# two-seg_3ry_search

An example implemenation of a two-segment tertiary epitope geometric match. 

The goal of this routine is to execute a fast search in a protein structure file (PDB format) for similar discontinuous peptide fragments separated by any number of residues. The search routine is completely sequence-agnostic and relies only on matching backbine dihedral angles profiles, and the relative internal orientation of the two segments. This routine can be run mutliple times sequentially to successively narrow a PDB library for searching three or more segments of discontinuous peptides. 

This routine was used to neofunctionalise simple and small proteins into granulopoietic agents by measuring their geometric compatibility as scaffolds for a migrated guest-epitope from human G-CSF.

## Dependencies 

  - NumPy
  
## Usage

```sh
$ cd dillinger
$ ./geo_mtch_fex.py -q qry_fn.pdb -s sbj_fn.pdb -d 0.2 -a 0.2 -1 10 -2 10 -b 1 -e 1 -n 2 -x 50
```
This attempts to find matches between an input 2-segment fragment defined by a backbone dihedrals vector <img src="https://render.githubusercontent.com/render/math?math=(\phi_1, \psi_1, .., \phi_n, \psi_n)">, and four centroid position vectors\; <img src="https://render.githubusercontent.com/render/math?math=\textbf{p_1}, \textbf{p_2},\textbf{p_3}, \textbf{p_4}">. These four centroids represent the centres of mass of the first residue in the first segment, the last residue in the first segment, the first residue on the second segment, and the last residue on the second segment. A sliding window with these segment sizes and a variable gap size is used to scan the subject structure for the segment pair that matches the backbone dihedrals and end-to-end orientation of the segments residues to that of the query. The search cutoffs thus define the maximum average dihedral deviation permitted, and the maximum average <img src="https://render.githubusercontent.com/render/math?math=\textbf{p_{1,2}} \leftrightarrow \textbf{p_{4,3}}"> distance permitted. 


input parameters: 

- -q query pdb file name\; /path/to/query_file.pdb
- -s subject pdb file name; /path/to/subject_file.pdb
- -d the maximum spacing permitted between points 1<->4 or 2<->3 in Angstroms
- -a average (phi, psi) absolute dihedral deviations permitted in radians
- -1 segment 1 sequence length; segment between points 1->2 <int>
- -2 segment 2 sequence length; segment between points 3->4 <int>
- -b nth order of the first amino acid order in segment 1 to be considered <int>
- -e reverse nth order of the last amino acid order in segment 2 to be considered <int>
- -n minimum sequence gap length <int>
- -x maximum sequence gab length <int>

## Example usage

```
$ ./geo_mtch_fex.py -q frg.pdb -s 2qup.pdb -d .4 -a 0.3 -1 10 -2 10 -b 1 -e 1 -n 2 -x 50
--backbone geometric search with params: 
query: frg.pdb 
subject: 2qup.pdb
spacing cutoff: 0.400000
dihedral deviation cutoff: 0.300000
segment 1 length: 10 
segment 2 length: 10 
segment 1 start position: 1 
segment end reverse position: 1
minimum sequence gap length: 2
maximum sequence gap length: 50

subject chain length: 119, abs dihedral dev: 0.254350, match start index 11, seg_1 length 10, gap length 31, seg_2 length 10
```
Where the last output line represents a signle hit found along a subject chain of length 119 residues. The match starts at the 11th amino acid and the optimal gap length is 31 amino acids, with 10 upstream and 10 downstream residues matched with an average dihedral deviation between query and subject of 0.25 radians, and the match fullfiled the interfragment average distance cutoff. 
