`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: NITT
// Engineer: 
// 
// Create Date: 12.03.2025 12:50:48
// Design Name: 
// Module Name: FrequencyAssignment
// Project Name: Task Schedhuler
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision: 2
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////
module FrequencyAssignment #(
    parameter N = 10,           // Number of tasks
    parameter L = 8,            // Number of frequency levels
    parameter Ebudget = 1000,   // Energy budget
    parameter Lg = 0,           // Leakage factor
    parameter I_sub = 0,        // Substrate current
    parameter V_bs = 0,         // Bias voltage
    parameter I_j = 0,          // Junction current
    parameter C_eff = 1,        // Effective capacitance
    // Fixed-point parameters
    parameter FRAC_BITS = 16,    // Number of fractional bits for fixed-point representation
    parameter NUM_PROCESSORS = 3
)(
    input clk,
    input reset,
    input [0:32*N-1] M,             // Execution cycles for each task
    input [0:32*L-1] f,             // Available frequency levels
    input [0:32*L-1] v,             // Corresponding voltage level
    input [0:32*N*L-1] gamma,       // Efficiency factor for each task-frequency pair
    input [0:32*N-1] D,             // Deadlines
    input [0:32*N-1] CPT,           // Critical path tasks
    input [0: NUM_PROCESSORS*N*32-1] exec_time_in,
    output reg [0:32*N-1] assigned_f, // Assigned frequency index
    output reg [31:0] Etotal,       // Total Energy consumption for this freq_assignment
    output reg valid                 // Indicates successful assignment
);

    // Fixed-point registers (Q16.16 format)
    reg [31:0] Emin, tcp, t_star;
    reg [31:0] E_mid, T_mid;
    reg [31:0] P_sta [0:L-1], P_dyn [0:L-1];
    reg [31:0] M_i;              // Execution cycles for each task
    reg [31:0] f_mid, v_mid, f_star;
    reg [31:0] gamma_il;         // Efficiency factor for each task-frequency pair
    reg [31:0] D_i;              // Deadlines
    reg [31:0] CPT_j;            // Critical path tasks
    reg [31:0] avg_exec_time [0: N-1];
    reg [31:0] exec_time [0:NUM_PROCESSORS-1][0:N-1];
    
    // Loop counters and control variables
    integer i, j;
    reg [31:0] low, high, mid, l_star;
    
    // Binary search state machine states
    localparam IDLE = 0,
               LOAD_TASK = 1,
               BINARY_SEARCH = 2,
               CALC_PARAMS = 3,
			   CHECK_AND_UPDATE = 4,
               FINALIZE_TASK = 5,
               COMPLETE = 6;
    
    reg [2:0] state, next_state;
    
    // Fixed-point multiplication function (Q16.16 format)
    function [31:0] fp_mult;
        input [31:0] a;
        input [31:0] b;
        reg [63:0] temp;
        begin
            temp = a*b;
            fp_mult = temp[47:16];
        end
      endfunction
    
    // Fixed-point division function (Q16.16 format)
    function [31:0] fp_div;
        input [31:0] a;
        input [31:0] b;
        begin
            // Scale up numerator before division to maintain precision
            fp_div =(a  / b) << FRAC_BITS;
        end
    endfunction
    
    // Convert integer to fixed-point
    function [31:0] int_to_fp;
        input [31:0] int_val;
        begin
            int_to_fp = int_val << FRAC_BITS;
        end
    endfunction
    
    always @(posedge clk or posedge reset) begin
        if (reset) begin
            
            for (i = 0; i < NUM_PROCESSORS; i = i + 1) begin
                for (j = 0; j < N; j = j + 1) begin
                    exec_time[i][j] = exec_time_in[32*(N*i+j)+:32];
                end
            end
            state <= IDLE;
            valid <= 0;
            for (i = 0; i < N; i = i + 1) begin
                assigned_f[32*i+:32] <= 0;
            end
            Emin <= 32'hFFFFFFFF; // Infinite
            Etotal <= 0;
            tcp <= 0;
            i <= 0;
        end else begin
            state <= next_state;
            end
   end
   
   always @(state) begin     
            case (state)
                IDLE: begin
                    i <= 0;
//                    valid <= 0;
                    if(valid == 0) begin
                        Etotal <= 0;
                        tcp <= 0;
                        next_state<= LOAD_TASK;
                    end
                end
                
                LOAD_TASK: begin
                    low <= 0;
                    high <= L-1;
                    Emin <= 32'hFFFFFFFF;
                    l_star <= 0;
                    M_i <= M[32*i+:32] << FRAC_BITS ; // Convert to fixed-point
                    D_i <= D[32*i+:32] ; // Convert to fixed-point
                    next_state<=BINARY_SEARCH;
                end
                
                BINARY_SEARCH: begin
                    if (low < high ) begin
                        mid <= (low + high) >> 1;
                        next_state <= CALC_PARAMS;
                    end else begin
                        next_state <= FINALIZE_TASK;
                    end
                end
                
                CALC_PARAMS: begin
                    f_mid = f[32*mid+:32] << FRAC_BITS ;
                    v_mid = v[32*mid+:32] ;
                    gamma_il = gamma[32*(L*i+mid)+:32] ;
//                    avg_exec_time[i] = 0;
//                    for (j = 0; j < NUM_PROCESSORS; j = j + 1) begin
//                        avg_exec_time[i] = avg_exec_time[i] + exec_time[j][i];
//                    end
//                    avg_exec_time[i] = avg_exec_time[i] / NUM_PROCESSORS;
//                        // Calculate T_mid = M_i / (gamma_il * f_mid)
                    T_mid = fp_div(M_i, fp_mult(gamma_il, f_mid));
                    
                    // Calculate P_sta and P_dyn
                    P_sta[mid] = fp_mult(int_to_fp(Lg), 
                                 fp_mult(v_mid, int_to_fp(I_sub)) + fp_mult(int_to_fp(V_bs), int_to_fp(I_j)));
                    //Ceff in order of 10^-8 and freq in MHZ -> /100 (*0.01)
                    P_dyn[mid] = fp_mult(int_to_fp(C_eff), 
                                 fp_mult(v_mid, fp_mult(v_mid, f_mid)))*0.01;
                    
                    // Calculate E_mid = T_mid * (P_sta[mid] + P_dyn[mid])
                    E_mid = fp_mult(T_mid, (P_sta[mid] + P_dyn[mid]))*0.01; //T in order of 10^-2
					next_state <= CHECK_AND_UPDATE;
				end
				
				CHECK_AND_UPDATE: begin
                  
                    if (Etotal + E_mid <= int_to_fp(Ebudget))begin
                        if(T_mid + tcp <= D_i) begin
                            if (E_mid < Emin) begin
                                l_star <= mid;
                                f_star <= f[32*mid+:32] << FRAC_BITS;
                                Emin <= E_mid;
                                t_star <= T_mid;
                                low <= mid + 1;
                                
                            end else begin
                                high <= mid;
                            end
                        end else begin
                            high <= mid;
                        end
                    end else begin
                        low <= mid + 1;
                    end
                    next_state <= BINARY_SEARCH;
                end
                
                
                FINALIZE_TASK: begin
                    if (Etotal+Emin <= int_to_fp(Ebudget)) begin
                        Etotal <= Etotal + Emin;
                        assigned_f[32*i+:32] <= l_star; // Convert back to integer
                        
                        // Update critical path time
                        for (j = 0; j < N; j = j + 1) begin
                            CPT_j = CPT[32*j+:32];
                            if (i == CPT_j) begin
                                tcp <= tcp + t_star;
                            end
                        end
                        
                        i <= i + 1;
                        if (i == N-1) begin
                            valid <= 1;
                            next_state <= COMPLETE;
                            end
                        else begin
                            next_state <= LOAD_TASK;
                        end
                     end
                     else begin
                        valid <= 0;
                        next_state <= COMPLETE;
                    end
                end
                
                COMPLETE: begin
                    next_state <= IDLE;
                end
                
                default: next_state <= IDLE;
            endcase
            
        end
endmodule




/*

div_gen_0 your_instance_name (
  .aclk(aclk),                                      // input wire aclk
  .s_axis_divisor_tvalid(s_axis_divisor_tvalid),    // input wire s_axis_divisor_tvalid
  .s_axis_divisor_tready(s_axis_divisor_tready),    // output wire s_axis_divisor_tready
  .s_axis_divisor_tdata(s_axis_divisor_tdata),      // input wire [31 : 0] s_axis_divisor_tdata
  .s_axis_dividend_tvalid(s_axis_dividend_tvalid),  // input wire s_axis_dividend_tvalid
  .s_axis_dividend_tready(s_axis_dividend_tready),  // output wire s_axis_dividend_tready
  .s_axis_dividend_tdata(s_axis_dividend_tdata),    // input wire [31 : 0] s_axis_dividend_tdata
  .m_axis_dout_tvalid(m_axis_dout_tvalid),          // output wire m_axis_dout_tvalid
  .m_axis_dout_tdata(m_axis_dout_tdata)            // output wire [47 : 0] m_axis_dout_tdata
);

*/