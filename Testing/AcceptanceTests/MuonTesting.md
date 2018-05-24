# Muon Testing



## Ionic Diffusion test

*Preparation*

* Files `51341.nxs`, `51342.nxs` and `51343.nxs`
* Make sure the location of these files is included in your search directories. Full instructions [are available online](http://www.mantidproject.org/MBC_Getting_set_up#MantidPlot_First-Time_Setup).

**Time required 5 - 10 minutes**

---

* Make sure that the `.nxs` files are in one of the directories on Mantid's search path.
* Open Muon Analysis
* Choose `Machine`: `EMU`
* In `Load Run` enter `51341` 
* Under `Settings` tab, ensure the `Enable Multiple Fitting` is checked
* In `Data Analysis Tab`
	* Beside `Runs` enter `51341-3` and check `Simultaneous` below
	* Right click the `Fit Function` table area; select `Add Function`
	* Add `Background` > `Flat Background`
	* Add `Muon` > `DynamicKuboToyabe`
	* Check the `Global` box for: `A0`, `Asym`, `Delta`, `Nu`
	* Click on the value `0.0000` for `Field` - a `...` box should appear, click it
		* On the `Set` dropdown, select `Fix All`
		* Enter values of `0`, `5`, `10` for the three runs
	* Set `A0 = 0.05`, `Asym = 0.15`, `Delta = 0.2`, `Nu = 0.1`
	* Click `Fit` > `Fit`
* In the `Results Table` tab
	* Under `Fitting Results` check the `Simultaneous` button
	* Click `Create Table`
	* A table should appear in the main Mantid work area
* In the `Results Table` workspace
	* Expected values (similar to) ..
	* f0.A0: 0.18235938124844153
	* f1.Asym: 0.0050597127151507
	* f1:Delta: 0.3301445010124604
	* f1:Nu: 0.33908659864509999

	


#### Possible problems

* Loading local data for 51341-3; these may not be in recognised directories - make sure that this is the case.


## Superconducting Copper Test

*Preparation*

* Files `EMU00020882.nxs` -- `EMU00020917.nxs`
* Make sure the location of these files is included in your search directories

**Time required 5 - 10 minutes**

---

* Make sure that the `.nxs` files are in one of the directories on Mantid's search path.
* Open Muon Analysis
* Choose `Machine`: `EMU`
* In `Load Run` enter `20882`
* Click the `>` tab next to `Load Run` a few times, this should produce new plots each time
* After loading `20884` switch to the `Data Analysis` tab
* Right click under `Property` and select `Add Function`
* Choose `Muon` > `ExpDecayMuon`
* Click `Fit` > `Fit`
* Click `Fit` > `Sequential Fit`
* Beside `Runs` enter `20889-20900`
* Click `start`
* Close the `Sequential Fits` window
* Switch to the `Results Table` tab
* In the `Values` table check `Run Number` and `Field_Danfysik`
* In the `Fitting Results` table check the `Sequential` radio button
* Ensure that all runs are checked in the `Fitting Results` table
* Click `Create Table`
* In the main Mantid GUI select the `ResultsTable` in the `Workspaces` tab, this is the one you just created
* Right click and select `Convert to MatrixWorkspace`
* In the resulting dialog select `Column X` as `Field_Danfysik` and `Column Y` as `lambda`
* Click `Run`
* Right click the new workspace created in the `Workspaces` tab
* Select `Plot Spectrum`
* You should get something like the following plot:
![](Cu-fitting.png)
