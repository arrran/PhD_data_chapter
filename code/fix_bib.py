#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:55:25 2020

@author: whitefar
"""


def fix_bib(bib_path):
            
        """
        Put curly brackets around all capitals in the titles of a .bib file e.g.
 
        
@article{smith2019land,
  title={Land ice height-retrieval algorithm for NASA's ICESat-2 photon-counting laser altimeter},
  author={Smith, Benjamin and Fricker, Helen A and Holschuh, Nicholas and Gardner, Alex S and Adusumilli, Susheel and Brunt, Kelly M and Csatho, Beata and Harbeck, Kaitlin and Huth, Alex and Neumann, Thomas and others},
  journal={Remote Sensing of Environment},
  volume={233},
  pages={111352},
  year={2019},
  publisher={Elsevier}
}
  
becomes:      
  
@article{smith2019land,
  title={{L}and ice height-retrieval algorithm for {N}{A}{S}{A}'s {I}{C}{E}{S}at-2 photon-counting laser altimeter},
  author={Smith, Benjamin and Fricker, Helen A and Holschuh, Nicholas and Gardner, Alex S and Adusumilli, Susheel and Brunt, Kelly M and Csatho, Beata and Harbeck, Kaitlin and Huth, Alex and Neumann, Thomas and others},
  journal={Remote Sensing of Environment},
  volume={233},
  pages={111352},
  year={2019},
  publisher={Elsevier}
}

        
         INPUT: bib_path
         a string. e.g. /Users/home/sam/sams_bib.bib
                
        
         OUTPUT: /Users/home/sam/sams_bib_caps.bib the .bib file with caps all in curly brackets
        """

        with open(bib_path,'r') as of:
            data = of.readlines()
            
            find_string = "title={"
                    
                #find line with this string
            line_i = [i for i, line in enumerate(data) if line.find(find_string) != -1 ]
                        
            for i in line_i:
                line = data[i]
                new_line = ""
                for letter in line:
                    if letter.isupper():
                        new_line += '{'+letter+'}'
                    else:
                        new_line += letter
                        
                data[i] = new_line
                
        with open(bib_path[:-4]+"_caps.bib",'w') as of:
            for line in data:
                of.write(line)
            print(f'written to {bib_path[:-4]}_caps.bib')
            