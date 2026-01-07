'''
A class for building hpc submission scripts.
'''


class ScriptBuilder():
    
    def __init__(self) -> None:
        self.script = ''
        self.header_bool = False
        return None
    
    def set_header(self, header_str: str) -> None:
        
        self.header_bool = True
        return None
    
    def add_cmd(self):
        return None
    
    def write_script(self, output_filename: str):
        return None
    
    def submit_script(self) -> None:
        return None 