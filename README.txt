This code implements explicit modelling of deformation twinning within a crystal plasticity framework. 
The framework is implemented in a user element subroutine VUEL within Abaqus Explicit.

Run:
To run this code run following command:
	abaqus job=twnbndBfull input=1000_Elem_Vuel_New.inp user=VUEL_Twin.f double=both mp_mode=threads cpus=10 inter

Post Process:
There is a post processing script written in python which writes vtk files based on outputs of simulation.

python post_process.py


Hooking up another UMAT with the twinning VUEL
- In the umat file, change the name of subroutine umat to umat_crys.
- Add rW, iterErr, Fe to arguments of umat_crys(......., kinc, rW, iterErr, Fe).
- Define additional variables that are missing. iterErr and Fe are already defined so define real*8 rW(3, 3, 3, 3). In addition also define real*8 dGamdTau(nslip) required by the new jacobian for vuel.
- Add dGamdTau as argument to solvcrys(......,dGamdTau). Make sure you do this at both implementation and call.
- Comment the call calc_jac(...) and add the following on the next line to call new jacobian:
        call calc_new_jac(dfgrd1, FpInv_n, props, nprops, C, nslip, Smat, dGamdTau, Cmat, dgamma, Fe, spk2, rW)
	(This is dsigma/dF required by vuel. It is based on Fe Fp type split)
- In VUEL_Twin.f, specify indices of euler angles in your material properties. Note this index corresponds solely to the material props and not the vuel props.
- In VUEL_Twin.f, change gb_elem_x, input_elmt_x, input_node_x to right files for the problem. Also comment the include to any existing umat (!include fcc_umat_lbl34_Ti.f).
- In the input file, remove all properties in the *UEL PROPERTY that don't have angled brackets. Add your material props in between. Make sure you maintain 8 props per line requried by Abaqus.
- Adjust the number of Properties and Variables in this line:
*USER ELEMENT, Nodes=8, Type=VU1, Iproperties=4, Properties=XX, Coordinates=3, Variables=XX
	Properties value should be 13 + number of your materials props
	Variables value should be 217 + 16 * number of state variables for material
- Lastly open the Makefile and add your umat file name
- Enter following command on the terminal to build your combined object file:
	make vuel_comb.o
- Finally run abaqus with the user combined object file
	abaqus job=twnbndBfull input=1000_Elem_Vuel_New.inp user=vuel_comb.o double=both mp_mode=threads cpus=1 gpus=0 inter

