`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 17.03.2025 14:26:12
// Design Name: 
// Module Name: TaskScheduler
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module TaskScheduler #(
    // Task Assignment parameters
    parameter NUM_TASKS = 10,
    parameter NUM_PROCESSORS = 3,
    parameter MAX_FREQ = 1500,
    parameter DATA_WIDTH = 32,
    parameter RANK_WIDTH = 16,
    
    // Frequency Assignment parameters
    parameter L = 8,              // Number of frequency levels
    parameter Ebudget = 10,     // Energy budget
    parameter Lg = 0,             // Leakage factor
    parameter I_sub = 0,          // Substrate current
    parameter V_bs = 0,           // Bias voltage
    parameter I_j = 0,            // Junction current
    parameter C_eff = 1,          // Effective capacitance
    parameter FRAC_BITS = 16      // Number of fractional bits for fixed-point representation
)(
    input clk,
    input reset,
    // Frequency Assignment inputs
    input [0:DATA_WIDTH*NUM_TASKS-1] M,            // Execution cycles for each task
    input [0:DATA_WIDTH*L-1] f,                    // Available frequency levels
    input [0:DATA_WIDTH*L-1] v,                    // Corresponding voltage level
    input [0:DATA_WIDTH*NUM_TASKS*L-1] gamma,      // Efficiency factor for each task-frequency pair
    input [0:DATA_WIDTH*NUM_TASKS-1] D,            // Deadlines
    input [0:DATA_WIDTH*NUM_TASKS-1] CPT,          // Critical path tasks
    
    // Task Assignment inputs
    input [0:NUM_TASKS*NUM_TASKS*DATA_WIDTH-1] comm_cost_in,
    input [0:NUM_PROCESSORS*NUM_TASKS*DATA_WIDTH-1] exec_time_in,
    
    // Outputs
    output [31:0] Etotal,
    output done,
    output [0:NUM_PROCESSORS*NUM_TASKS-1] processor_assignment,
    output [DATA_WIDTH-1:0] Energy_total 
);

    // Internal signals
    wire freq_assignment_valid;
    wire [0:DATA_WIDTH*NUM_TASKS-1] assigned_frequencies;
    
    // Instantiate Frequency Assignment module
    FrequencyAssignment #(
        .N(NUM_TASKS),
        .L(L),
        .Ebudget(Ebudget),
        .Lg(Lg),
        .I_sub(I_sub),
        .V_bs(V_bs),
        .I_j(I_j),
        .C_eff(C_eff),
        .FRAC_BITS(FRAC_BITS)
        
    ) freq_assign_inst (
        .clk(clk),
        .reset(reset),
        .M(M),
        .f(f),
        .v(v),
        .gamma(gamma),
        .D(D),
        .CPT(CPT),
        .assigned_f(assigned_frequencies),
        .Etotal(Etotal),
        .valid(freq_assignment_valid),
        .exec_time_in(exec_time_in)
    );
    
    // Instantiate Task Assignment module
    task_assignment #(
        .NUM_TASKS(NUM_TASKS),
        .NUM_PROCESSORS(NUM_PROCESSORS),
        .MAX_FREQ(MAX_FREQ),
        .DATA_WIDTH(DATA_WIDTH),
        .RANK_WIDTH(RANK_WIDTH),
         .Lg(Lg),
        .I_sub(I_sub),
        .V_bs(V_bs),
        .I_j(I_j),
        .C_eff(C_eff),
        .FRAC_BITS(FRAC_BITS),
        .L(L)
    ) task_assign_inst (
        .clk(clk),
        .reset(reset),
        .v(v),
        .f(f),
        .start(freq_assignment_valid),  // Start task assignment when frequency assignment is valid
        .comm_cost_in(comm_cost_in),
        .exec_time_in(exec_time_in),
        .freq_in(assigned_frequencies),  // Use assigned frequencies from the frequency assignment module
        .done(done),
        .processor_assignment(processor_assignment),
        .Energy_total(Energy_total)
    );

endmodule
