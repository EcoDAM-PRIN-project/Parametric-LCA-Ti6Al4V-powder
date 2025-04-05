# Ti6Al4V Powder Production LCA Optimizer
Comprehensive parametric Life Cycle Assessment (LCA) model for Ti6Al4V powder production aiming to define the process parameters minimizing the environmental impact of the process under specific conditions and constraints

<b>Christian Spreafico (1*), Baris Ördek‬ (1)</b><br>
(1) Department of Management Information and Production Engineering, University of Bergamo, Viale Marconi 5, 24044 Dalmine (Bg), Italy.
(*) christian.spreafico@unibg.it

## 📌 Features

- Multi-objective optimization using `fmincon`
- Support for 18 different environmental impact categories
- Country-dependent LCA coefficients (EU and China)
- Full breakdown of material and energy usage
- Calculation of dependent variables and process flows
- Modular structure for future extension

## 🛠️ Technologies Used

- **Language**: MATLAB
- **Optimization**: Sequential Quadratic Programming (`fmincon`)
- **License**: GNU General Public License v3.0

## 🧮 Input Parameters

The tool requests several user inputs upon running:

1. `m_finalPowder`: Final powder mass to be produced `[kg]`
2. `impact_category`: Impact category to minimize  
   (choose from: `TA`, `GW`, `FET`, `MET`, `TET`, `FF`, `ME`, `HTPc`, `HTPnc`, `IR`, `LO`, `SO`, `OD`, `PMF`, `HOF`, `EOF`, `WC`)
3. `country`: Region for country-specific LCA coefficients (`EU` or `CN`)
4. `d_target`: Target particle diameter `[µm]`

## 🧠 Optimization Variables

The tool optimizes the following variables within pre-defined constraints:

- `φ_ATelectrode` (Electrode diameter) `[m]`
- `p_ATargon` (Atomization pressure) `[MPa]`
- `β_TiO2` (TiO₂ content in Ti slag) `[-]`

The objective is to minimize total environmental impact across the selected impact category.

## 📤 Output

Upon completion, the tool displays:

- Optimized values for `φ`, `p_ATargon`, `β_TiO2`
- Mass flows of intermediate and final products
- Energy consumption for each stage
- Total and stage-wise environmental impacts


## 📚 Code Structure

- `optimizeLCA()` - Main function for parameter input, optimization execution, and results display
- `objectiveFunction(x, params)` - Computes total impact and breakdown
- `constraintFunction(x, params)` - Ensures all process constraints are satisfied
- `calculateDependentVariables(x_opt, params)` - Computes physical and process parameters

## 🔐 License

This project is licensed under the **GNU GPL v3.0**.  
You may freely use, distribute, and modify it under the license terms.


For full license details, see [LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html).

## 🙋‍♀️ Authors

- **Christian Spreafico**  
- **Baris Ördek**

## 📩 Contact

For academic collaborations or questions, feel free to open an issue or reach out to the authors.

---

🔬 *Optimizing sustainability for a cleaner titanium future.*
