`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 15.04.2025 14:42:24
// Design Name: 
// Module Name: GDES_tb
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


module GDES_tb( );
    // Parameters from the original module
    parameter NUM_TASKS = 10;
    parameter NUM_PROCESSORS = 3;
    parameter MAX_FREQ = 1500;
    parameter DATA_WIDTH = 32;
    parameter RANK_WIDTH = 16;
    
    parameter L = 8;
    parameter E_given_G = 100;
    parameter Lg = 1;
    parameter I_sub = 1;
    parameter V_bs = 1;
    parameter I_j = 1;
    parameter C_eff = 1;  // Adjusted from 0.000000001 to avoid precision issues in fixed point
    parameter FRAC_BITS = 16;  // Same as in the design
    
    
    //reg [0:32*NUM_TASKS-1] M;
    reg [0:32*L-1] f;
    reg [0:32*L-1] v;
    //reg [0:32*NUM_TASKS-1] D;
    

    
    // Testbench signals
    reg clk;
    reg reset;
    reg start;
    reg [DATA_WIDTH-1:0] D_G;
    reg [0: NUM_TASKS*NUM_TASKS*DATA_WIDTH-1] comm_cost_in;
    reg [0: NUM_PROCESSORS*NUM_TASKS*DATA_WIDTH-1] exec_time_in;
    //reg [0: NUM_TASKS*DATA_WIDTH-1] freq_in;
    wire done;
    wire [0:3*NUM_TASKS-1] processor_assignment;
    wire [0:DATA_WIDTH*NUM_TASKS-1] task_assignment;
    wire [0:DATA_WIDTH*NUM_TASKS-1] frequency_assignment;
    wire [31:0] E_total;
    //wire valid;
    
    // Instantiate the Unit Under Test (UUT)
    gdes #(
        .NUM_TASKS(NUM_TASKS),
        .NUM_PROCESSORS(NUM_PROCESSORS),
        .MAX_FREQ(MAX_FREQ),
        .DATA_WIDTH(DATA_WIDTH),
        .RANK_WIDTH(RANK_WIDTH),
        .L(L),
        .E_given_G(E_given_G),
        .C_eff(C_eff),
        .Lg(Lg),
        .I_sub(I_sub),
        .V_bs(V_bs),
        .I_j(I_j),
        .FRAC_BITS(FRAC_BITS)
    ) uut (
        .clk(clk),
        .reset(reset),
        .start(start),
        .f(f),
        .v(v),
        .comm_cost_in(comm_cost_in),
        .exec_time_in(exec_time_in),
        .D_G(D_G),
        //.freq_in(freq_in),
        .done(done),
        .frequency_assignment(frequency_assignment),
        .energy_assignment(energy_assignment),
        .processor_assignment(processor_assignment),
        .E_total(E_total)
    );
    
    // Clock generation
    initial begin
        clk = 0;
        forever #5 clk = ~clk; // 100MHz clock
    end
    
     function [31:0] real_to_fixed;
        input real real_val;
        begin
            real_to_fixed = $rtoi(real_val * (2.0 ** FRAC_BITS));
        end
    endfunction
     function real fixed_to_real;
        input [31:0] fixed_val;
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
            $display("Total Energy: %f", fixed_to_real(E_total));
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
        
        D_G = 106;
        
        // Communication costs matrix from the task graph in Image 1
        // Format: comm_cost[predecessor][successor]
        // Each row represents a task, each column represents a successor task
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
        
        // Execution times from Image 2 (Table 2)
        // Format: exec_time[processor][task]
        exec_time_in = {
            // Processor u1 (execution times for tasks n1-n10)
            32'd15, 32'd14, 32'd12, 32'd14, 32'd13, 32'd14, 32'd8,  32'd6,  32'd19, 32'd22,
            // Processor u2 (execution times for tasks n1-n10)
            32'd17, 32'd20, 32'd14, 32'd9,  32'd14, 32'd17, 32'd16, 32'd12, 32'd13, 32'd8,
            // Processor u3 (execution times for tasks n1-n10)
            32'd10,  32'd19, 32'd20, 32'd18, 32'd11, 32'd10,  32'd12, 32'd15, 32'd21, 32'd17
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
