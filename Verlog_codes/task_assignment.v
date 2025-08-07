`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 25.02.2025 19:43:00
// Design Name: 
// Module Name: task_assignment
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


module task_assignment #(
    parameter NUM_TASKS = 10,
    parameter NUM_PROCESSORS = 3,
    parameter L = 8,            // Number of frequency levels
    parameter MAX_FREQ = 1500, // Representing 100% of maximum frequency
    parameter DATA_WIDTH = 32,  // Width for execution time and comm cost
    parameter RANK_WIDTH = 16,  // Width for rank and time calculations
    parameter Lg = 1,             // Leakage factor
    parameter I_sub = 1,          // Substrate current
    parameter V_bs = 1,           // Bias voltage
    parameter I_j = 1,            // Junction current
    parameter C_eff = 1,          // Effective capacitance
    parameter FRAC_BITS = 16  
)(
    input clk,
    input reset,
    input start,
    input [0:DATA_WIDTH*L-1] f,             // Available frequency levels
    input [0:L*DATA_WIDTH-1] v,             // Corresponding voltage level
    input [0: NUM_TASKS*NUM_TASKS*DATA_WIDTH-1] comm_cost_in,
    input [0: NUM_PROCESSORS*NUM_TASKS*DATA_WIDTH-1] exec_time_in,
    input [0: NUM_TASKS*DATA_WIDTH-1] freq_in,
    output reg done,
    output reg [0:3*NUM_TASKS-1] processor_assignment, // Stores which processor each task is assigned to
    output reg [DATA_WIDTH-1:0] Energy_total
);

    // Graph structure - Communication costs between tasks
    // Format: comm_cost[predecessor][successor]
    reg [DATA_WIDTH-1:0] comm_cost [0:NUM_TASKS-1][0:NUM_TASKS-1];
    
    // Execution times for each task on each processor
    reg [DATA_WIDTH-1:0] exec_time [0:NUM_PROCESSORS-1][0:NUM_TASKS-1];
    
    // Stores the rank of each task
    reg [RANK_WIDTH-1:0] task_rank [0:NUM_TASKS-1];
    
    // Sorted tasks based on rank (indices of tasks)
    reg [3:0] sorted_tasks [0:NUM_TASKS-1];
    
    // Processor availability time
    reg [RANK_WIDTH-1:0] processor_avail [0:NUM_PROCESSORS-1];
    
    // Actual Finish Time for each task
    reg [DATA_WIDTH-1:0] AFT [0:NUM_TASKS-1];
    
    reg [DATA_WIDTH-1:0] freq [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] volt [0:NUM_TASKS-1];
    
    // State machine states
    reg [3:0] state;
    parameter IDLE = 4'd0;
    parameter LOAD_DATA = 4'd1;
    parameter CALC_RANKS = 4'd2;
    parameter SORT_TASKS = 4'd3;
    parameter ASSIGN_TASKS = 4'd4;
    parameter DONE = 4'd5;
    
    // Temporary variables for computation
    reg [3:0] current_task;
    reg [3:0] current_processor;
    reg [RANK_WIDTH-1:0] est, eft, min_eft;
    reg [3:0] best_processor;
    reg [3:0] task_idx;
    reg signed [4:0] i, j, k;
    reg [RANK_WIDTH-1:0] avg_exec_time [0: NUM_TASKS-1];
    reg [RANK_WIDTH-1:0] max_successor_rank;
    reg [RANK_WIDTH-1:0] temp_rank;
    reg [3:0] temp_task;
    reg [RANK_WIDTH-1:0] max_predecessor_time;
    reg [RANK_WIDTH-1:0] latest_task;
    reg [RANK_WIDTH-1:0] latest_finish_time;
    
    reg [DATA_WIDTH-1:0] P_static [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] P_dynamic [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] Energy [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] freq_mult_ratio;
    reg [DATA_WIDTH-1:0] freq_mult;
    //reg [DATA_WIDTH-1:0] Energy_total;
    
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
    
    // Load data from inputs
    task load_data;
        begin
            // Load communication costs
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                for (j = 0; j < NUM_TASKS; j = j + 1) begin
                    comm_cost[i][j] = comm_cost_in[DATA_WIDTH*(NUM_TASKS*i+j)+:DATA_WIDTH];
                end
            end
            
            // Load execution times
            for (i = 0; i < NUM_PROCESSORS; i = i + 1) begin
                for (j = 0; j < NUM_TASKS; j = j + 1) begin
                    exec_time[i][j] = exec_time_in[DATA_WIDTH*(NUM_TASKS*i+j)+:DATA_WIDTH];
                end
            end
            
            //Load the frequency of each task
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                freq[i] = f[DATA_WIDTH*(freq_in[DATA_WIDTH*i+: DATA_WIDTH])+:DATA_WIDTH];
            end
            
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                volt[i] = v[DATA_WIDTH*(freq_in[DATA_WIDTH*i+: DATA_WIDTH])+:DATA_WIDTH];
            end
            
            
            // Initialize processor availability
            for (i = 0; i < NUM_PROCESSORS; i = i + 1) begin
                processor_avail[i] = 0;
            end
            
            // Initialize AFT
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                AFT[i] = 0;
            end
            
            
        end
    endtask
    
    // Task to calculate ranks of all tasks
    task calculate_ranks;
        begin
            // First calculate average execution time for each task
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                avg_exec_time[i] = 0;
                for (j = 0; j < NUM_PROCESSORS; j = j + 1) begin
                    avg_exec_time[i] = avg_exec_time[i] + exec_time[j][i];
                end
                avg_exec_time[i] = avg_exec_time[i] / NUM_PROCESSORS;
                
                // Initial rank is just the average execution time (will be updated)
                task_rank[i] = avg_exec_time[i];
            end
            
            // Start with exit node and work backwards
            // Multiple passes for rank calculation (bottom-up approach)
//            for (k = 0; k < 3; k = k + 1) begin // 3 passes should be enough for most graphs
                for (i = NUM_TASKS-1; i >= 0; i = i - 1) begin
                    avg_exec_time[i] = 0;
                    for (j = 0; j < NUM_PROCESSORS; j = j + 1) begin
                        avg_exec_time[i] = avg_exec_time[i] + exec_time[j][i];
                    end
                    avg_exec_time[i] = avg_exec_time[i] / NUM_PROCESSORS;
                    
                    max_successor_rank = 0;
                    
                    // Find maximum (comm_cost + rank) among successors
                    for (j = 0; j < NUM_TASKS; j = j + 1) begin
                        if (comm_cost[i][j] > 0) begin // If j is a successor of i
                            if (comm_cost[i][j] + task_rank[j] > max_successor_rank) begin
                                max_successor_rank = comm_cost[i][j] + task_rank[j];
                            end
                        end
                    end
                    
                    // Update rank
                    task_rank[i] = avg_exec_time[i] + max_successor_rank;
                end
//            end
        end
    endtask
    
    // Task to sort tasks based on their ranks (descending order)
    task sort_tasks_by_rank;
        begin
            // Initialize sorted_tasks
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                sorted_tasks[i] = i;
            end
            
            // Simple bubble sort
            for (i = 0; i < NUM_TASKS-1; i = i + 1) begin
                for (j = 0; j < NUM_TASKS-1-i; j = j + 1) begin
                    if (task_rank[sorted_tasks[j]] < task_rank[sorted_tasks[j+1]]) begin
                        temp_task = sorted_tasks[j];
                        sorted_tasks[j] = sorted_tasks[j+1];
                        sorted_tasks[j+1] = temp_task;
                    end
                end
            end
        end
    endtask
    
    // Assign tasks to processors
    task assign_tasks;
        begin
            for (task_idx = 0; task_idx < NUM_TASKS; task_idx = task_idx + 1) begin
                current_task = sorted_tasks[task_idx];
                min_eft = {RANK_WIDTH{1'b1}}; // Maximum value for the bit width
                best_processor = 0;
                
                // Try each processor
                for (current_processor = 0; current_processor < NUM_PROCESSORS; current_processor = current_processor + 1) begin
                    
                    latest_task = 0;
                    latest_finish_time = 0;
                    
                    for (i = 0; i < NUM_TASKS; i = i + 1) begin
                        if (processor_assignment[NUM_PROCESSORS*i+: NUM_PROCESSORS] == current_processor) begin
                            if (AFT[i] > latest_finish_time) begin
                                latest_finish_time = AFT[i];
                                latest_task = i;
                            end
                        end
                    end
                    
                    // Update availability time for processor y
                    processor_avail[current_processor] = AFT[latest_task];
                    
                    // Calculate EST
                    if (current_task == 0) begin
                        // Entry task
                        est = 0;
                    end else begin
                    
                    max_predecessor_time = 0;

                    // Find maximum of (AFT[predecessor] + comm_cost[predecessor][current_task]) among all predecessors
                    for (i = 0; i < NUM_TASKS; i = i + 1) begin
                        if (comm_cost[i][current_task] > 0) begin // If i is a predecessor of current_task
                        // Add communication cost only if tasks run on different processors
                            if (processor_assignment[NUM_PROCESSORS*i+: NUM_PROCESSORS] != current_processor) begin
                                if (AFT[i] + comm_cost[i][current_task] > max_predecessor_time) begin
                                    max_predecessor_time = AFT[i] + comm_cost[i][current_task];
                                end
                            end else begin
                            // No communication cost if on same processor
                            if (AFT[i] > max_predecessor_time) begin
                                max_predecessor_time = AFT[i];
                            end
                            end
                        end
                    end
                    
                    
                    
                    // Final EST calculation - max of processor availability and max predecessor time
                    est = (processor_avail[current_processor] > max_predecessor_time) ?  processor_avail[current_processor] : max_predecessor_time;
                  
                    end
                    
                    // Calculate EFT
                    eft = est + ((exec_time[current_processor][current_task] * MAX_FREQ) / freq[current_task]);//check
                    
                    // Check if this is the best processor so far
                    if (eft < min_eft) begin
                        min_eft = eft;
                        best_processor = current_processor;
                    end
                end
                
                // Assign the task to the best processor
                processor_assignment[NUM_PROCESSORS*current_task+: NUM_PROCESSORS] = best_processor;
                
                P_static[current_task] = fp_mult(int_to_fp(Lg), fp_mult(volt[current_task], int_to_fp(I_sub)) + fp_mult(int_to_fp(V_bs), int_to_fp(I_j)));
                    //Ceff in order of 10^-8 and freq in MHZ -> /100 (*0.01)
                P_dynamic[current_task] = fp_mult(int_to_fp(C_eff), fp_mult(volt[current_task], fp_mult(volt[current_task], freq[current_task]<<FRAC_BITS)))*0.01;
                
                freq_mult = (MAX_FREQ<<FRAC_BITS)*10;
                freq_mult_ratio = fp_div(freq_mult,freq[current_task]<<FRAC_BITS);
                
                Energy[current_task] = fp_mult((P_static[current_task]+P_dynamic[current_task]),fp_mult(exec_time [best_processor][current_task]<<FRAC_BITS, freq_mult_ratio))*0.01 *0.1;    
                Energy_total = Energy_total + Energy[current_task]; 
                // Update processor availability
                // processor_avail[best_processor] = min_eft;
                
                // Update AFT for this task
                
                AFT[current_task] = min_eft; // change - find what
            end
        end
    endtask
    
    initial begin
        // Initialize processor assignments
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                processor_assignment[NUM_PROCESSORS*i+: NUM_PROCESSORS] = 0;
            end
    end
    
    
    // Main state machine
    always @(posedge clk or posedge reset) begin
        if (reset) begin
            Energy_total = 0;
            state <= IDLE;
            done <= 0;
        end else begin
            case (state)
                IDLE: begin
                    if (start) begin
                        state <= LOAD_DATA;
                        done <= 0;
                    end
                end
                
                LOAD_DATA: begin
                    load_data();
                    state <= CALC_RANKS;
                end
                
                CALC_RANKS: begin
                    calculate_ranks();
                    state <= SORT_TASKS;
                end
                
                SORT_TASKS: begin
                    sort_tasks_by_rank();
                    state <= ASSIGN_TASKS;
                end
                
                ASSIGN_TASKS: begin
                    assign_tasks();
                    state <= DONE;
                end
                
                DONE: begin
                    done <= 1;
                    // Stay in DONE state until reset
                end
                
                default: state <= IDLE;
            endcase
        end
    end

endmodule