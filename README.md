# Reactive-Distillation-with-Pyomo
Hi there, welcome!

In this project, our goal is to build a **Fischer–Tropsch reactive distillation simulator** using **_Pyomo_**.

## Model Navigation
**Code Map:**
* **archive:** Model archives, completed with all necessary files to guarantee compatibility. Might not always reflect the latest changes.
* **global_sets:** A place for global sets (e.g `COMP_TOTAL`). Each sub-module in this project will be referencing the same sets of components and universal parameters.
* **saved_solutions:** A place for .json or .pickle files, usually saved after certain pyomo solve, will be used either as a quick access archive or simply a temporary file.
* **utility:** A place for helper functions we wrote to assist data gathering, model construction, debugging, analysis and data visulization.
* **data:** This module contains all the `physical` parameters for the model. For each of the pyomo blocks, we recommend writing only `one python script`. This way each pyomo blocks only need to import once from the data module. For large sets of data, use of excel files are recommended by using `xlrd` inside corresponding data python script.

* **physics:** This module is often refered to as `1st level physics block` in this project. Here resides all the fundemental physcis equations (or `block construction rules`) that will be used in higher level model structures. Each module will be using certain `parent_block variables` and `local variables` depending on its nature. Currently we have the following physics modules:
  1. **Environment Variables:** Each physics block is responsible for defining its own sets, parameters and variables, however, there are certain `environment variables` that represent shared environment values. **Note:** These `environment variables` will be defined during construction of `2nd level stage block`.
      * Example Usage:  

      ```python
      # In 2nd level stage block
      block.T = pe.Var(within=pe.NonNegativeReals,bounds=(200+273.15,300+273.15)) # K

      # In 1st level physics block
      def k_FT_rule(block):
        return block.k_FT == k.k0_FT * pe.exp(-k.E_FT*1e3/(k.R*block.parent_block().T))
      block.k_KT_con = pe.Constraint(rule=k_FT_rule)
      ```

      * List of `shared environment variables`:

      | Description | Name | Description | Name |
      | --- | :---: | --- | :---: |
      | Temperature | T | Feed Temperature | T_F |
      | Pressure | P | Feed Flow Rate| F |
      | Catalyst Loading | cat| Feed Composition | z |
      | Heat Duty | Q_main | Feed Enthalpy | H_F |
      | Vapor Outflow Rate | V | Liquid Outflow Rate | L |
      | Vapor Molar % | y | Liquid Molar % | x |
      | Vapor Enthalpy | H_V | Liquid Enthalpy | H_L |
      | Vapor Fugacity | f_V | Liquid Fugacity | f_L |
      | Reaction Rates | r_total_comp |


  2. **kinetics:** Given catalyst loading, gas phase fugacity, temperature and pressure, calculate each component's reaction rate.
  3. **enthalpy:** Given gas and liquid composition, temperature and pressure, return molar enthalpy of gas and liquid phase.
  4. **VLE:** Given gas and liquid composition, temperature and pressure, return vapor and liquid fugacity.
  5. **VLLE:** Similar to VLE, water is treated as a seperate components, with its own set of rules.


* **stages:** This module is often refered as `2nd level stage block` in this project. Here resides different types of single trays. Based on tray type (e.g reactive / non-reactive), physics equations in `1st level physics block` are selected and assmbled using `MESH` equations. To retain the ability to customize column configuration, **flow variables in this module is entirely local**, and will be linked manually in a higher level structure.
  1. **Assemble Physics:**
      * Example Usage:  

      ```python
      # In 1st level physics block
      def kinetic_block_rule(block):
        # ... ... ... ... (relevant sets, parameters, variables and equations)

      # In 2nd level stage block
      block.kinetics_block = pe.Block(rule=kinetic_block_rule)
      ```
  2. **Reactive Stage:** Equivalent to a reactive flash set-up.
  3. **Reboiler:** Equivalent to a non-reactive flash set-up.
  4. **Condenser Stage:** Currently obsolete, treat water as Henry components.
  5. **Condenser Stage2:** Treat water as a separate component and a seperate phase.

* **columns:** This module is refered as `3rd level block` in this project. Here `2nd level stage blocks` are assembled and linked to simulate a reactive distillation column. For better support of presentations and data analysis, files in this section use  `Jupyter Notebook` format.

* **Showcase:** Any completed model will be transfered to this section.


```python
# mass balance around each tray
def mass_balance(model,j,i):
    return model.F[j]*model.z[j,i] + model.r_total_comp[j,i] \
           + model.L[j-1]*model.x[j-1,i] + model.V[j+1]*model.y[j+1,i]
           == model.L[j]*model.x[j,i] + model.V[j]*model.y[j,i]
# For all trays and all components
model.mass_balance = pe.Constraint(model.TRAY,model.COMP_TOTAL,rule=mass_balance)
```

```python
def gamma_rule(model,j,i):
    return pe.log(model.gamma[j,i])*(model.n_ave[j] - u.cal_cnumber('C6H14')) == \
            pe.log(model.gamma_ref[j])*(model.n_ave[j] - u.cal_cnumber(i))
model.gamma_con = pe.Constraint(model.TRAY | model.COND | model.REBO, model.COMP_NONHENRY,rule=gamma_rule)
```


```python
# deactivate henry components constraints
model.f_L_HENRY_con.deavtivate()
model.Hen_con.deactivate()
model.Hen_ref_con.deactivate()
# ...
# activate activity coefficient constraints
model.f_L_NONHENRY_con.activate()
model.gamma_con.activate()
model.P_sat_con.activate()
# ...
```
```python
# In 1st level physics block
def VLE_block_rule(block):
    # declare sets
    block.COMP_HENRY = pe.Set(initialize=['H2','CO','CO2','H2O','C1H4','C2H6','C3H8','C2H4','C3H6'])
    block.COMP_NONHENRY = m.COMP_TOTAL - block.COMP_HENRY
    # ...    
    # declare variables
    block.Hen = pe.Var(block.COMP_HENRY,within=pe.NonNegativeReals,bounds=Hen_bounds)  # Bar
    block.gamma = pe.Var(block.COMP_NONHENRY,within=pe.PositiveReals,initialize=0.1,bounds=gamma_bounds)
    block.P_sat = pe.Var(block.COMP_NONHENRY,within=pe.PositiveReals,initialize=1e-10,bounds=P_sat_bounds)  # Bar
    block.poynting = pe.Var(m.COMP_TOTAL,within=pe.Reals,bounds=poynting_bounds)
    # ...
    # declare equations
    # henry component
    def f_L_HENRY_rule(block,i):
        return block.parent_block().f_L[i] == block.Hen[i]*block.parent_block().x[i]*block.poynting[i]
    block.f_L_HENRY_con = pe.Constraint(block.COMP_HENRY,rule=f_L_HENRY_rule)

    # henry's constant
    def Hen_rule(block,i):
        return block.Hen[i]*pe.exp(block.n_ave*e.henry.dHen[i]) == pe.exp(block.Hen0[i])
    block.Hen_con = pe.Constraint(block.COMP_HENRY,rule=Hen_rule)


# Kinetics Block
block.kinetics_block = pe.Block(rule=kinetic_block_rule)
# Energy Block
block.energy_block = pe.Block(rule=energy_block_rule)
# VLE block
block.VLE_block = pe.Block(rule=VLE_block_rule)
```
