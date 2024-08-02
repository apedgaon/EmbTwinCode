umat_file=hcp_umat_TDT_pq_Ti.f
abq_ver=2023

VUEL_Twin.o: VUEL_Twin.f GSPT_3D_8N.f SHP_FUNC_3D_8N.f prep_umat.f gdefgrad_3D.f Gstrain_3D.f calc_new_jac.f
	ifort -c -fPIC -auto -extend-source -WB -I%I -I/usr/SIMULIA/EstProducts/${abq_ver}/SMAUsubs/PublicInterfaces -fpp -O2 -qopenmp VUEL_Twin.f
${umat_file}.o: ${umat_file}
	ifort -c -fPIC -auto -extend-source -WB -I%I -I/usr/SIMULIA/EstProducts/${abq_ver}/SMAUsubs/PublicInterfaces -fpp -O2 -qopenmp ${umat_file} -o ${umat_file}.o
vuel_comb.o: VUEL_Twin.o ${umat_file}.o
	xiar rc vuel_comb.o VUEL_Twin.o ${umat_file}.o
clean:
	rm *.o