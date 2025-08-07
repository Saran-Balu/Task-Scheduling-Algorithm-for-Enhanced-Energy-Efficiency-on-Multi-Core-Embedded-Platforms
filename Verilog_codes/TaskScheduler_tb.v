






`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 17.03.2025 15:30:00
// Design Name: 
// Module Name: TaskScheduler_tb
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: Testbench for the integrated TaskScheduler module
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////

module TaskScheduler_tb();
    // Parameters
    parameter NUM_TASKS = 10;
    parameter NUM_PROCESSORS = 3;
    parameter MAX_FREQ = 1500;
    parameter DATA_WIDTH = 32;
    parameter RANK_WIDTH = 16;
    parameter L = 8;
    parameter Ebudget = 25;  // Using the value from FrequencyAssignment_tb
    parameter Lg = 1;
    parameter I_sub = 1;
    parameter V_bs = 1;
    parameter I_j = 1;
    parameter C_eff = 1;
    parameter FRAC_BITS = 16;
    
    // Testbench signals
    reg clk;
    reg reset;
    reg start;
    
    // Frequency Assignment inputs
    reg [0:DATA_WIDTH*NUM_TASKS-1] M;
    reg [0:DATA_WIDTH*L-1] f;
    reg [0:DATA_WIDTH*L-1] v;
    reg [0:DATA_WIDTH*NUM_TASKS*L-1] gamma;
    reg [0:DATA_WIDTH*NUM_TASKS-1] D;
    reg [0:DATA_WIDTH*NUM_TASKS-1] CPT;
    
    // Task Assignment inputs
    reg [0:NUM_TASKS*NUM_TASKS*DATA_WIDTH-1] comm_cost_in;
    reg [0:NUM_PROCESSORS*NUM_TASKS*DATA_WIDTH-1] exec_time_in;
    
    //Counter for time when processes are done
    reg FA_done, TA_done;
    // Outputs
    wire [31:0] Etotal;
    wire done;
    wire [0:NUM_PROCESSORS*NUM_TASKS-1] processor_assignment;
    wire [DATA_WIDTH-1:0] Energy_total;
    
    // Function to convert floating point to fixed point
    function [31:0] real_to_fixed;
        input real real_val;
        begin
            real_to_fixed = $rtoi(real_val * (2.0 ** FRAC_BITS));
        end
    endfunction
    
    // Function to convert fixed point to real
    function real fixed_to_real;
        input [31:0] fixed_val;
        begin
            fixed_to_real = $itor(fixed_val) / (2.0 ** FRAC_BITS);
        end
    endfunction
    
    // Instantiate the Unit Under Test (UUT)
    TaskScheduler #(
        .NUM_TASKS(NUM_TASKS),
        .NUM_PROCESSORS(NUM_PROCESSORS),
        .MAX_FREQ(MAX_FREQ),
        .DATA_WIDTH(DATA_WIDTH),
        .RANK_WIDTH(RANK_WIDTH),
        .L(L),
        .Ebudget(Ebudget),
        .Lg(Lg),
        .I_sub(I_sub),
        .V_bs(V_bs),
        .I_j(I_j),
        .C_eff(C_eff),
        .FRAC_BITS(FRAC_BITS)
    ) uut (
        .clk(clk),
        .reset(reset),
        .M(M),
        .f(f),
        .v(v),
        .gamma(gamma),
        .D(D),
        .CPT(CPT),
        .comm_cost_in(comm_cost_in),
        .exec_time_in(exec_time_in),
        .Etotal(Etotal),
        .Energy_total(Energy_total),
        .done(done),
        .processor_assignment(processor_assignment)
    );
    
    // Clock generation
    initial begin
        clk = 0;
        forever #5 clk = ~clk; // 100MHz clock
    end
    
    // Task to display results
    task display_results;
        begin
            $display("TaskScheduler Results:");
            $display("------------------------");
            $display("Total Energy after TA: %f", fixed_to_real(Energy_total));
            $display("Total Energy: %f", fixed_to_real(Etotal));
            $display("Task | Frequency | Processor");
            $display("-----|-----------|----------");
            for (integer i = 0; i < NUM_TASKS; i = i + 1) begin
                $display("  %d  |    %d    |    %d", 
                         i+1, 
                         f[DATA_WIDTH*(assigned_frequencies[DATA_WIDTH*i+: DATA_WIDTH])+:DATA_WIDTH], 
                         processor_assignment[NUM_PROCESSORS*i+:NUM_PROCESSORS]);
            end
        end
    endtask
    
    // For debugging - monitor assigned frequencies from the Frequency Assignment module
    wire [0:DATA_WIDTH*NUM_TASKS-1] assigned_frequencies;
    assign assigned_frequencies = uut.assigned_frequencies;
    
    integer task_idx, freq_idx;
    real gamma_value;
    
    // Test stimulus
    initial begin
        // Initialize inputs
        reset = 1;
        start = 0;
        FA_done = 0;
        TA_done = 0;
        // Initialize FrequencyAssignment inputs
        // Execution cycles for each task from FrequencyAssignment_tb
        // Execution cycles for each task in order of 10^4
        M = {
            32'd19500,  // Task n1 
            32'd25000,  // Task n2
            32'd21500,  // Task n3
            32'd19000,  // Task n4
            32'd17500,  // Task n5
            32'd19000,  // Task n6
            32'd16500,  // Task n7
            32'd15000,  // Task n8
            32'd20000,  // Task n9
            32'd22000   // Task n10 
        };
        
        // Available frequency levels in MHz
        f = {
            32'd1500, 32'd1300, 32'd1100, 32'd900,
            32'd667, 32'd600, 32'd500, 32'd400             
        };
        
        // Corresponding voltage levels
        v = {
            real_to_fixed(1.35), // 1500 MHz
            real_to_fixed(1.30), // 1300 MHz
            real_to_fixed(1.25), // 1100 MHz
            real_to_fixed(1.15), // 900 MHz
            real_to_fixed(1.05), // 667 MHz
            real_to_fixed(1.00), // 600 MHz
            real_to_fixed(0.95), // 500 MHz
            real_to_fixed(0.90)  // 400 MHz
        };
        
        // Initialize gamma matrix with values between 0.4 and 1.0
        for (task_idx = 0; task_idx < NUM_TASKS; task_idx = task_idx + 1) begin
            for (freq_idx = 0; freq_idx < L; freq_idx = freq_idx + 1) begin
                // Gamma value between 0.4 and 1.0
                gamma_value = 0.8 + (0.2 * (L - freq_idx) / L) + (0.1 * (task_idx % 3) / 3);
                if (gamma_value > 1.0) gamma_value = 1.0;
                
                // Convert to fixed point and assign
                gamma[DATA_WIDTH*(task_idx*L + freq_idx)+:DATA_WIDTH] = real_to_fixed(gamma_value);
                
//                // Display for debugging
//                $display("Task %0d, Freq %0d: gamma = %f (fixed: %h)", 
//                         task_idx+1, freq_idx+1, gamma_value, 
//                         real_to_fixed(gamma_value));
            end
        end
        
        // Deadlines from FrequencyAssignment_tb
        D = {
//            real_to_fixed(27),  // Task 1
//            real_to_fixed(47),  // Task 2
//            real_to_fixed(79),  // Task 3
//            real_to_fixed(126),  // Task 4
//            real_to_fixed(123),  // Task 5
//            real_to_fixed(126),  // Task 6
//            real_to_fixed(91),  // Task 7
//            real_to_fixed(117),  // Task 8
//            real_to_fixed(121),  // Task 9
//            real_to_fixed(130)   // Task 10

            real_to_fixed(26),  // Task 1
            real_to_fixed(55),  // Task 2
            real_to_fixed(53),  // Task 3
            real_to_fixed(80),  // Task 4
            real_to_fixed(79),  // Task 5
            real_to_fixed(80),  // Task 6
            real_to_fixed(79),  // Task 7
            real_to_fixed(105),  // Task 8
            real_to_fixed(106),  // Task 9
            real_to_fixed(106)   // Task 10
        };
        
        // Critical path tasks
        CPT = {
            32'd0,  // Task 1 is on critical path
            32'd2,  // Task 3 is on critical path
            32'd6,  // Task 7 is on critical path
            32'd9,  // Task 10 is on critical path
            32'd0, 32'd0, 32'd0, 32'd0, 32'd0, 32'd0  // Padding
        };
        
        // Communication costs matrix from task_assignment_test
        comm_cost_in = {
            // n1 -> n2, n3, n4, n5, n6
            32'd0,  32'd18, 32'd12, 32'd9,  32'd11, 32'd14, 32'd0,  32'd0,  32'd0,  32'd0,
            // n2 -> n8, n9
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd19, 32'd16, 32'd0,  32'd0,
            // n3 -> n7
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd23,  32'd0, 32'd0, 32'd0,
            // n4 -> n8,n9
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd27, 32'd23,  32'd0,
            // n5 -> n9
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0, 32'd13, 32'd0,
            // n6 -> n8
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd15,  32'd0, 32'd0,
            // n7 -> n10
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd17,
            // n8 -> n10
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd11,
            // n9 -> n10
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd13,
            // n10 has no successors
            32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0,  32'd0
        };
        
        // Execution times from task_assignment_test
        exec_time_in = {
            // Processor u1 (execution times for tasks n1-n10)
            32'd14, 32'd13, 32'd11, 32'd13, 32'd12, 32'd13, 32'd7,  32'd5,  32'd18, 32'd21,
            // Processor u2 (execution times for tasks n1-n10)
            32'd16, 32'd19, 32'd13, 32'd8,  32'd13, 32'd16, 32'd15, 32'd11, 32'd12, 32'd7,
            // Processor u3 (execution times for tasks n1-n10)
            32'd9,  32'd18, 32'd19, 32'd17, 32'd10, 32'd9,  32'd11, 32'd14, 32'd20, 32'd16
        };
        
        // Apply reset
        #20;
        reset = 0;
        
        // Wait for the process to complete
        wait(done);
        #50;
        
        // Display results
        display_results();
        
        // End simulation
        #10;
        $finish;
    end
    
    // Monitor state changes
    always @(posedge clk) begin
        
        if (uut.freq_assign_inst.valid)begin
            if(FA_done == 0)begin
                FA_done = 1; 
                $display("Frequency Assignment completed at time %t", $time);
            end
        end
        
        if (done) begin
            if (TA_done == 0) begin
                TA_done = 1;
                $display("Task Assignment completed at time %t", $time);
            end
        end
    end
    
    // Main monitor
    initial begin
        $monitor("Time: %t, Reset: %b, Start: %b, FreqValid: %b, Done: %b", 
                 $time, reset, start, uut.freq_assignment_valid, done);
    end
    
endmodule