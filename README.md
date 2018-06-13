# Reactive-Distillation-with-Pyomo
Hi there, welcome!

In this project, our goal is to build a **Fischer–Tropsch reactive distillation simulator** using **_Pyomo_**.

## Model Navigation
**Code Map:**
* **global_sets:** A place for global sets (e.g `COMP_TOTAL`). Each sub-module in this project will be referencing the same sets of components and universal parameters.
* **saved_states:** A place for .json or .pickle files, usually saved after certain pyomo solve, will be used either as a quick access archive or simply a temporary file.
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
