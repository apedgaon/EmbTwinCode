VUEL_Twin.o: VUEL_Twin.f
	ifort -c -fPIC -auto -extend-source -WB -I%I -I/usr/SIMULIA/EstProducts/2022/SMAUsubs/PublicInterfaces -fpp -O2 -qopenmp VUEL_Twin.f
fcc_umat_lbl34_Ti.o: fcc_umat_lbl34_Ti.f
	ifort -c -fPIC -auto -extend-source -WB -I%I -I/usr/SIMULIA/EstProducts/2022/SMAUsubs/PublicInterfaces -fpp -O2 -qopenmp fcc_umat_lbl34_Ti.f
vuel_comb.o: VUEL_Twin.o fcc_umat_lbl34_Ti.o
	xiar rc vuel_comb.o VUEL_Twin.o fcc_umat_lbl34_Ti.o
clean:
	rm *.o
