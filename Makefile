umat_file=fcc_umat_lbl34_Ti.f

VUEL_Twin.o: VUEL_Twin.f
	ifort -c -fPIC -auto -extend-source -WB -I%I -I/usr/SIMULIA/EstProducts/2022/SMAUsubs/PublicInterfaces -fpp -O2 -qopenmp VUEL_Twin.f
${umat_file}.o: ${umat_file}
	ifort -c -fPIC -auto -extend-source -WB -I%I -I/usr/SIMULIA/EstProducts/2022/SMAUsubs/PublicInterfaces -fpp -O2 -qopenmp ${umat_file} -o ${umat_file}.o
vuel_comb.o: VUEL_Twin.o ${umat_file}.o
	xiar rc vuel_comb.o VUEL_Twin.o ${umat_file}.o
clean:
	rm *.o
