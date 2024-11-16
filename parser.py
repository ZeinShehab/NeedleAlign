import sys
from typing import List

class Parser:
    def __init__(self) -> None:
        pass
    
    def parse_file(self, filename: str) -> List[str]:
        strings: List[str] = []
        file = open(filename, 'r')
        
        while True:
            line = file.readline()
            
            if line.startswith(">"):
                continue
            
            if len(line) <= 0:
                return strings
            
            strings.append(line.strip())