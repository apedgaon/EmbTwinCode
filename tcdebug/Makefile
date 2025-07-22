FC = ifort
FFLAGS = -g

tcnr.o: tcnr.f
	$(FC) $(FFLAGS) -c tcnr.f
trac_cont.o: ../trac_cont.f
	$(FC) $(FFLAGS) -c ../trac_cont.f
hcp_umat_TDT_pq_Ti.o: ../hcp_umat_TDT_pq_Ti_5.f
	$(FC) $(FFLAGS) -c ../hcp_umat_TDT_pq_Ti_5.f -o hcp_umat_TDT_pq_Ti.o
calc_new_jac.o: ../calc_new_jac.f
	$(FC) $(FFLAGS) -c ../calc_new_jac.f
a.out: tcnr.o trac_cont.o hcp_umat_TDT_pq_Ti.o calc_new_jac.o
	$(FC) $(FFLAGS) -o a.out tcnr.o trac_cont.o hcp_umat_TDT_pq_Ti.o calc_new_jac.o
clean:
	rm -f *.o a.out