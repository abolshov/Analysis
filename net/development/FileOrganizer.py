import pathlib
import re
import os
from typing import Dict

class FileOrganizer:
    """
    Class for grouping files corresponding to each masspoint
    into a dictionary of the form
        masspoint: [list of paths to files for this masspoint]
    """
    def __init__(self,
                 start_dir: str | os.PathLike | pathlib.Path,
                 pattern: re.Pattern[str],
                 verbose: bool=False):

        self.verbose = verbose
        if self.verbose:
            print(f'Initilizing FileOrganizer with pattern={pattern} and start_dir={start_dir}.')
            
        self.pattern = re.compile(pattern)
        self.start_dir = pathlib.Path(start_dir)
        self.file_map = self._search(self.start_dir)

    def __getitem__(self, 
                    mp: int) -> list[pathlib.Path]:
        try:
            return self.file_map[mp]
        except KeyError:
            raise KeyError(f'Files corresponding to mass M={mp} not found in {self.start_dir}.')

    def __iter__(self):
        return iter(self.file_map)

    def _search(self, 
                start_dir: pathlib.Path) -> Dict[int, list[pathlib.Path]]:
        
        if self.verbose:
            print('Start search.')

        file_map = {}
        for file_path in start_dir.rglob('*.root'):
            file_path = file_path.resolve()
            match = self.pattern.search(str(file_path))
            if match:
                mass = int(match.group(1))
                if mass in file_map:
                    file_map[mass].append(file_path)
                else:
                    file_map[mass] = [file_path]
                
                if self.verbose:
                    print(f'M={mass}, file={file_path}')

        if self.verbose:
            print('Search done.')

        return file_map