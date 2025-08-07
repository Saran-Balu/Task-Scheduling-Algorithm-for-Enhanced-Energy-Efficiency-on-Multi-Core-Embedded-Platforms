`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 15.04.2025 12:19:17
// Design Name: 
// Module Name: isaecc_algo_tb
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


module isaecc_algo_tb();
    // Parameters from the original module
    parameter NUM_TASKS = 10;
    parameter NUM_PROCESSORS = 3;
    parameter MAX_FREQ = 1500;
    parameter DATA_WIDTH = 64;
    parameter RANK_WIDTH = 64;
    
    parameter L = 8;
    parameter E_given_G =50;
    parameter Lg = 1;
    parameter I_sub = 1;
    parameter V_bs = 1;
    parameter I_j = 1;
    parameter C_eff = 558345748e1;  // Adjusted from 0.000000001 to avoid precision issues in fixed point
    parameter FRAC_BITS = 32;  // Same as in the design
    
    
    //reg [0:32*NUM_TASKS-1] M;
    reg [0:DATA_WIDTH*L-1] f;
    reg [0:DATA_WIDTH*L-1] v;
    //reg [0:32*NUM_TASKS-1] D;
    

    
    // Testbench signals
    reg clk;
    reg reset;
    reg start;
    reg [0: NUM_TASKS*NUM_TASKS*DATA_WIDTH-1] comm_cost_in;
    reg [0: NUM_PROCESSORS*NUM_TASKS*DATA_WIDTH-1] exec_time_in;
    //reg [0: NUM_TASKS*DATA_WIDTH-1] freq_in;
    wire done;
    wire [0:3*NUM_TASKS-1] processor_assignment;
    wire [0:DATA_WIDTH*NUM_TASKS-1] task_assignment;
    wire [0:DATA_WIDTH*NUM_TASKS-1] frequency_assignment;
    wire [DATA_WIDTH-1:0] E_total;
    //wire valid;
    
    // Instantiate the Unit Under Test (UUT)
    isaecc_algo #(
        .NUM_TASKS(NUM_TASKS),
        .NUM_PROCESSORS(NUM_PROCESSORS),
        .MAX_FREQ(MAX_FREQ),
        .DATA_WIDTH(DATA_WIDTH),
        .RANK_WIDTH(RANK_WIDTH),
        .L(L),
        .E_given_G(E_given_G),
        .C_eff(C_eff)
    ) uut (
        .clk(clk),
        .reset(reset),
        .start(start),
        .comm_cost_in(comm_cost_in),
        .exec_time_in(exec_time_in),
        .f(f),
        .v(v),
        //.freq_in(freq_in),
        .done(done),
        .processor_assignment(processor_assignment),
        .frequency_assignment(frequency_assignment),
        .energy_assignment(energy_assignment),
        .E_total(E_total)
    );
    
    // Clock generation
    initial begin
        clk = 0;
        forever #5 clk = ~clk; // 100MHz clock
    end
    
     function [DATA_WIDTH-1:0] real_to_fixed;
        input real real_val;
        begin
            real_to_fixed = real_val * (2.0 ** FRAC_BITS);
        end
    endfunction
     function real fixed_to_real;
        input [DATA_WIDTH-1:0] fixed_val;
        begin
            fixed_to_real = $itor(fixed_val) / (2.0 ** FRAC_BITS);
        end
    endfunction
    
    // Function to convert floating point to fixed point integer
   
    
    // Task to display results
    task display_results;
        begin
            $display("TaskScheduler Results:");
            $display("------------------------");
            $display("Total Energy: %f", fixed_to_real(E_total>>7));
            $display("Task | Frequency | Processor");
            $display("-----|-----------|----------");
            for (integer i = 0; i < NUM_TASKS; i = i + 1) begin
                $display("  %d  |    %d    |    %d", 
                         i+1, 
                         fixed_to_real(frequency_assignment[DATA_WIDTH*i+:DATA_WIDTH]), 
                         processor_assignment[NUM_PROCESSORS*i+:NUM_PROCESSORS]);
            end
        end
    endtask
    
    // Test stimulus
    initial begin
        // Initialize inputs
        f=0;
        v=0;
       // D=0;
       
       
        reset = 1;
        start = 0;
        
        // Available frequency levels in MHz
        f = {
            64'd1500, 64'd1300, 64'd1100, 64'd900,
            64'd667, 64'd600, 64'd500, 64'd400             
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
        
        // Communication costs matrix from the task graph in Image 1
        // Format: comm_cost[predecessor][successor]
        // Each row represents a task, each column represents a successor task
        comm_cost_in = {
            // n1 -> n2, n3, n4, n5, n6
            64'd0,  64'd18, 64'd12, 64'd9,  64'd11, 64'd14, 64'd0,  64'd0,  64'd0,  64'd0,
            // n2 -> n8, n9
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0, 64'd19, 64'd16,  64'd0,
            // n3 -> n7
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd23,  64'd0, 64'd0, 64'd0,
            // n4 -> n8,n9
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd27, 64'd23,  64'd0,
            // n5 -> n9
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0, 64'd13, 64'd0,
            // n6 -> n8
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd15,  64'd0, 64'd0,
            // n7 -> n10
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd17,
            // n8 -> n10
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd11,
            // n9 -> n10
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd13,
            // n10 has no successors
            64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0,  64'd0
        };
        
        // Execution times from Image 2 (Table 2)
        // Format: exec_time[processor][task]
        exec_time_in = {
            // Processor u1 (execution times for tasks n1-n10)
             64'd15, 64'd14, 64'd12, 64'd14, 64'd13, 64'd14, 64'd8,  64'd6,  64'd19, 64'd22,
            // Processor u2 (execution times for tasks n1-n10)
            64'd17, 64'd20, 64'd14, 64'd9,  64'd14, 64'd17, 64'd16, 64'd12, 64'd13, 64'd8,
            // Processor u3 (execution times for tasks n1-n10)
            64'd10,  64'd19, 64'd20, 64'd18, 64'd11, 64'd10,  64'd12, 64'd15, 64'd21, 64'd17

        };
        
        // Frequency for each task (assuming 100% maximum frequency for all tasks initially)
        // In a real application, this would vary by task requirements
        
        
        // Apply reset
        #20;
        reset = 0;
        
        // Start the task assignment process
        #5;
        start = 1;
        #10;
        start = 0;
        
        // Wait for the process to complete
        wait(done);
        #10;
        
        // Display results
        display_results();
        
        
        
        $finish;
    end
    
    // Monitor changes in the state for debugging
    initial begin
        $monitor("Time: %0t, Reset: %b, Start: %b, Done: %b", 
                 $time, reset, start, done);
    end
    
endmodule
