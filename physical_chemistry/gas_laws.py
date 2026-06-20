import numpy as np
import pandas as pd
from sympy import symbols, Eq, solve, fraction

# ideal gas equations

class ideal_gas_equation():
    
    def __init__(self):
        self.P, self.V, self.n, self.R, self.T = symbols("P V n R T")
        self.ideal_gas_eq = Eq(self.P * self.V, self.n * self.R * self.T)

        self.variable_map = {
            "P": self.P,
            "V": self.V,
            "n": self.n,
            "R": self.R,
            "T": self.T
        }

        print(str(self.P) + str(self.V) + " = " + str(self.n) \
              + str(self.R) + str(self.T))
        
    def solve_for(self, target_variable, values = None):
        target_str = target_variable.upper()
        target_symbol = self.variable_map.get(target_str)
        
        if target_symbol is None:
            raise ValueError(f"Invalid variable. Choose from: {list(self.variable_map.keys())}")
        
        solution = solve(self.ideal_gas_eq, target_symbol)
        solution_expr = solution[0]
        
        num, denom = fraction(solution_expr)
        
        if denom != 1:
            equation_str = f"{target_str} = ({num})/({denom})"
        else:
            equation_str = f"{target_str} = {solution_expr}"
        
        if values:
            sympy_subs = {self.variable_map[k]: v for k, v in values.items() if k \
                          in self.variable_map}
            numeric_solution = solution_expr.subs(sympy_subs)
            return numeric_solution, equation_str
        

        return solution, equation_str


class van_der_waals():

    def __init__(self):
        self.P, self.a, self.n, self.V, self.b, self.n, self.R, self.T = \
        symbols("P a n V b n R T")
        self.vdw_eq = Eq(self.P + ((self.a * self.n**2)/self.V**2) * (self.V - self.n * self.b), 
                    self.n * self.R * self.T)
        
        self.variable_map = {
            "P": self.P,
            "a": self.a,
            "b": self.b,
            "V": self.V,
            "n": self.n,
            "R": self.R,
            "T": self.T
        }

        print("(" + str(self.P) + " + ((" + str(self.a) + str(self.n) + "^2))" + "/" + str(self.V) + "^2)" \
              + "(" + str(self.V) + " - " + str(self.n) + str(self.b) + ") = " \
                + str(self.n) + str(self.R) + str(self.T))
        
    def solve_for(self, target_variable, values = None):
        target_str = target_variable.upper()
        target_symbol = self.variable_map.get(target_str)
        
        if target_symbol is None:
            raise ValueError(f"Invalid variable. Choose from: {list(self.variable_map.keys())}")
        
        solution = solve(self.vdw_eq, target_symbol)
        solution_expr = solution[0]
        
        num, denom = fraction(solution_expr)
        
        if denom != 1:
            equation_str = f"{target_str} = ({num})/({denom})"
        else:
            equation_str = f"{target_str} = {solution_expr}"
        
        if values:
            sympy_subs = {self.variable_map[k]: v for k, v in values.items() if k \
                          in self.variable_map}
            numeric_solution = solution_expr.subs(sympy_subs)
            return numeric_solution, equation_str
        

        return solution, equation_str

ige_dict = {
    "P" : 6,
    "V" : 1,
    "n" : 2,
    "R" : 3,
    "T" : 4
}


vdw_dict = {
    "a" : 1.35,
    "b" : 0.0386,
    "V" : 800,
    "n" : 200,
    "R" : 0.08206,
    "T" : 898
}

# ig = ideal_gas_equation()
# p = ig.solve_for("T", ige_dict)
# print(p)

vdw = van_der_waals()
p = vdw.solve_for("P", vdw_dict)
print(p)