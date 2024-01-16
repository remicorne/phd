import pandas as pd
from abc import abstractmethod
from dataclasses import dataclass
from module.utils import (
    subselectDf,
    save_figure,
    parallel_process,
    generate_figure,
    dictToFilename,
)
from module.Matrices import Matrices
from module.Matrix import Matrix


@dataclass
class Figure:
    
    data: pd.DataFrame
    
    @abstractmethod
    def generate_fig(self):
        pass
    
    @abstractmethod
    def generate_fig_data(self):
        pass
    
    @abstractmethod
    def plot_fig_data(self):
        pass
    
    @abstractmethod
    def generate_save_filepath(self):
        pass
    
    @abstractmethod
    def save(self):
        

@dataclass
class SingleCorrelogram:
    ax: object
    matrix: Matrix
    
    
    
    

@dataclass
class Correlogram:
    matrices: Matrices = None
    filename: str = None
    group_by: str = None
    between: str = None
    variables: str | list = None
    accross: str = None
    columns: list[str] = None
    sub_selector: str = None
    n_minimum: int = 5
    method: str = "pearson"
    pvalue_threshold: float = 0.05
    from_scratch: bool = False

    def __post_init__(self):
        if self.matrices is None:
            self.matrices = Matrices(
                filename=self.filename,
                group_by=self.group_by,
                between=self.between,
                variables=self.variables,
                accross=self.accross,
                columns=self.columns,
                sub_selector=self.sub_selector,
                n_minimum=self.n_minimum,
                method=self.method,
                pvalue_threshold=self.pvalue_threshold
            )

    def plot(self, from_scratch=True):
        # [plotting code here, similar to plot_correlogram in Matrices]
        # Use self.matrices.matrices for data
    
    