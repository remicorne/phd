from dataclasses import dataclass
from module.core.ProjectMember import ProjectMember
from matplotlib import pyplot as plt
from typing import ClassVar
from threading import Thread
from IPython.display import Image, display

@dataclass
class Figure(ProjectMember):
    
    from_scratch = None
    _subfolder: ClassVar[str] = None
    
    def __post_init__(self):
        super().__post_init__()
        self.figure_parameters = None
        if self.from_scratch:
            self.get()
            
    def generate(self):
        fig = self.generate_figure()
        return self.plot(fig)
    
    def generate_figure(self):
        return plt.subplots(self.figure_parameters)

    def plot(self):
        pass
            
    def initialize(self):
        fig = super().initialize()
        fig.show()
    
    def save(self, fig):
        def target():
            fig.savefig(f"{self.filepath}.svg")  # dpi also?
            print(f"SAVED {self.filepath}.svg")
            fig.savefig(f"{self.filepath}.png", dpi=fig.dpi)
            print(f"SAVED {self.filepath}.png")
        Thread(target=target).start()
        
    def load(self):
        display(Image(filename=f'{self.filepath}.png'))
        
    @property
    def identifier(self):
        pass
    
    @property
    def filepath(self):
        return f"{self.location}/{self.identifier}"
    
    
class Correlogram(Figure):
    