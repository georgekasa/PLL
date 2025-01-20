//"verilogams"

`include "constants.vams"
`include "disciplines.vams"

`timescale 1s/1fs


module reference_veriloga #( 
	parameter real F_REF = 192e6,//ref
	parameter real PN_100KHz = -150,
	parameter real DLY_INIT = 100e-9,
	parameter integer MY_SEED_ref = 29
	)
	(
	output reg  out = 1'b0
	);
	
	
	
	//local parameter
	localparam real T_REF = 1.0/F_REF;
	localparam real Jitter_rms = T_REF/2.0/`M_PI * sqrt(F_REF * 10**(PN_100KHz/10.0));
	
	integer  seed_ref = MY_SEED_ref;
	real randn = 0;
	real jitter = 0;
	real jitter_pre = 0;
	real jitter_non_accu = 0;
	real ref_period;
	
	
	
	initial begin
		#DLY_INIT
		forever begin
		// non - accumulative jitter
		randn = $dist_normal(seed_ref, 0, 1000) * 1e-3;
		jitter = randn*Jitter_rms;
		jitter_non_accu = jitter -jitter_pre;
		jitter_pre = jitter;
		
		
		// clk with accumulative jitter
		ref_period = 1.0/F_REF + jitter_non_accu;
		
		#(0.5*ref_period)
		out = 1'b1;
		
		#(0.5*ref_period)
		out = 1'b0;
		end
	end
		
		
endmodule
