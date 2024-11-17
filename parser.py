import sys
from typing import List

class Parser:
    def __init__(self) -> None:
        pass
    
    def parse_file(self, filename: str) -> List[str]:
        with open(filename, "r") as file:
            nodes_seq = {}
            node_name = ""
            
            lines = file.readlines()
            
            for i in range(len(lines)):
                line = lines[i].strip()
                
                if line.startswith(">"):
                    node_name = line.split(' ')[0][1:]
                    
                    if node_name not in nodes_seq:
                        nodes_seq[node_name] = ""        
                else:
                    nodes_seq[node_name] += line
            
            nodes_split = {}
                
            for node in nodes_seq.keys():
                sequence = nodes_seq[node]
                
                seq_list = [sequence[i:i+50] for i in range(0, len(sequence), 50)]
                nodes_split[node] = seq_list
        return nodes_split  
    
    def parse_nosplit(self, filename):
        with open(filename, "r") as file:
            nodes_seq = {}
            node_name = ""
            
            lines = file.readlines()
            
            for i in range(len(lines)):
                line = lines[i].strip()
                
                if line.startswith(">"):
                    node_name = line.split(' ')[0][1:]
                    
                    if node_name not in nodes_seq:
                        nodes_seq[node_name] = ""       
                else:
                    nodes_seq[node_name] += line
            
        return nodes_seq 