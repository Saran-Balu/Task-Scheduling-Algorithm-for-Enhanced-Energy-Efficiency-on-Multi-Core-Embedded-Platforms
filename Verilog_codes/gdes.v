`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 15.04.2025 12:17:24
// Design Name: 
// Module Name: gdes
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


module gdes #(
    parameter NUM_TASKS = 10,
    parameter NUM_PROCESSORS = 3,
    parameter L = 8,            // Number of frequency levels
    parameter E_given_G = 1000,   // Energy budget
    parameter C_eff = 1,        // Effective capacitance
    parameter MAX_FREQ = 1500, // Representing 100% of maximum frequency
    parameter DATA_WIDTH = 32,  // Width for execution time and comm cost
    parameter RANK_WIDTH = 16,  // Width for rank and time calculations
     parameter FRAC_BITS = 16,
     parameter Lg = 1,
    parameter I_sub = 1,
    parameter V_bs = 1,
    parameter I_j = 1 
)(
    input clk,
    input reset,
    input start,
    input [0:DATA_WIDTH*L-1] f,             // Available frequency levels
    input [0:DATA_WIDTH*L-1] v,             // Corresponding voltage level
    //input [0:32*NUM_TASKS-1] D,             // Deadlines
    input [0: NUM_TASKS*NUM_TASKS*DATA_WIDTH-1] comm_cost_in,
    input [0: NUM_PROCESSORS*NUM_TASKS*DATA_WIDTH-1] exec_time_in,
//    input [0:32*N-1] D,             // Deadlines
    input [DATA_WIDTH-1:0] D_G,               //Deadline of application
    output reg done,
    output reg [0:DATA_WIDTH*NUM_TASKS-1] frequency_assignment,
    output reg [0:DATA_WIDTH*NUM_TASKS-1] energy_assignment,
    output reg [0:3*NUM_TASKS-1] processor_assignment, // Stores which processor each task is assigned to
    output reg [31:0] E_total           // Energy consumption
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
    reg [DATA_WIDTH-1:0] processor_avail [0:NUM_PROCESSORS-1];
    
    // Actual Finish Time for each task
    reg [DATA_WIDTH-1:0] AFT [0:NUM_TASKS-1], AST [0:NUM_TASKS-1];
    
    reg [DATA_WIDTH-1:0] freq [0:L-1];
    reg [DATA_WIDTH-1:0] volt [0:L-1];
    reg [DATA_WIDTH-1:0] E [0:NUM_TASKS-1];
    
    // State machine states
    reg [3:0] state;
    parameter IDLE = 4'd0;
    parameter LOAD_DATA = 4'd1;
    parameter CALC_RANKS = 4'd2;
    parameter SORT_TASKS = 4'd3;
    parameter CALC_PARAM = 4'd4;
    parameter ASSIGN_TASKS = 4'd5;
    parameter DONE = 4'd6;
    
    // Temporary variables for computation
    reg [3:0] current_task;
    //reg [3:0] current_processor;
    reg [DATA_WIDTH-1:0] est [0:NUM_TASKS -1][0:NUM_PROCESSORS-1], lft [0:NUM_TASKS-1][0:NUM_PROCESSORS-1], eft [0:NUM_TASKS-1][0:NUM_PROCESSORS-1], min_eft;
    reg [DATA_WIDTH-1:0] MET [0:NUM_TASKS -1][0:NUM_PROCESSORS-1],UBET [0:NUM_TASKS -1][0:NUM_PROCESSORS-1], MET_low [0:NUM_TASKS -1][0:NUM_PROCESSORS-1];
    reg [DATA_WIDTH-1:0] f_low [0:NUM_TASKS -1][0:NUM_PROCESSORS-1], E_low [0:NUM_TASKS -1][0:NUM_PROCESSORS-1];
   
   // reg [RANK_WIDTH-1:0] E_given_G;
    
    reg [3:0] best_processor = -1;
    reg [3:0] proc = -1;
    reg [3:0] task_idx;
    reg signed [4:0] i, j, k, u;
    reg [DATA_WIDTH-1:0] avg_exec_time [0: NUM_TASKS-1];
    reg [RANK_WIDTH-1:0] max_successor_rank;
    reg [RANK_WIDTH-1:0] temp_rank;
    reg [3:0] temp_task;
    reg [DATA_WIDTH-1:0] max_predecessor_time = 0, min_successor_time = 32'hFFFFFFFF;
    reg [RANK_WIDTH-1:0] latest_task;
    reg [DATA_WIDTH-1:0] latest_finish_time;
    reg [DATA_WIDTH-1:0] E_min;
    
    reg [DATA_WIDTH-1:0] P_sta [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] P_dyn [0:NUM_TASKS-1];
     
    reg [DATA_WIDTH-1:0] E_given [0:NUM_TASKS-1]; 
    
    
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
    
    
    function [DATA_WIDTH-1:0] nearest_freq;
        input [DATA_WIDTH-1:0] input_f;
        begin
            for(k=0;k<L;k=k+1)begin
                if (freq[k]>=input_f)
                    nearest_freq = freq[k];
            end
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
            for (i = 0; i < L; i = i + 1) begin
                freq[i] = f[DATA_WIDTH*i+: DATA_WIDTH]  << FRAC_BITS;
            end
            
            //Load the Voltage of each task
            for (i = 0; i < L; i = i + 1) begin
                volt[i] = v[DATA_WIDTH*i+: DATA_WIDTH];
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
                    if (task_rank[sorted_tasks[j]] > task_rank[sorted_tasks[j+1]]) begin
                        temp_task = sorted_tasks[j];
                        sorted_tasks[j] = sorted_tasks[j+1];
                        sorted_tasks[j+1] = temp_task;
                    end
                end
            end
        end
    endtask
    
    //Calculate energies
    task calculate_est;
        begin
            for (task_idx = 0; task_idx < NUM_TASKS; task_idx = task_idx + 1) begin
                min_eft = {DATA_WIDTH{1'b1}};
                for (u=0; u< NUM_PROCESSORS ; u=u+1)begin
                    latest_task = 0;
                    latest_finish_time = 0;
                    
//                    for (i = 0; i < NUM_TASKS; i = i + 1) begin
//                        if (processor_assignment[NUM_PROCESSORS*i+: NUM_PROCESSORS] == u) begin
//                            if (AFT[i] > latest_finish_time) begin
//                                latest_finish_time = AFT[i];
//                                latest_task = i;
//                            end
//                        end
//                    end
                    
                    // Update availability time for processor y
//                    processor_avail[u] = AFT[latest_task];
                    
                    // Calculate EST
                    if (task_idx == 0) begin
                        // Entry task
                        est[task_idx][u] = 0;
                    end 
                    else begin
                    
                        max_predecessor_time = 0;
    
                        // Find maximum of (AFT[predecessor] + comm_cost[predecessor][task_idx]) among all predecessors
                        for (i = 0; i < NUM_TASKS; i = i + 1) begin
                            if (comm_cost[i][task_idx] > 0) begin // If i is a predecessor of current_task
                            // Add communication cost only if tasks run on different processors
                                if (processor_assignment[NUM_PROCESSORS*i+: NUM_PROCESSORS] != u) begin
                                    if (AFT[i] + comm_cost[i][task_idx] > max_predecessor_time) begin
                                        max_predecessor_time = AFT[i] + comm_cost[i][task_idx];
                                    end
                                end
                                else begin
                                // No communication cost if on same processor
                                    if (AFT[i] > max_predecessor_time) begin
                                        max_predecessor_time = AFT[i];
                                    end
                                end
                            end
                        end
                    end
                    
                    
                    // Final EST calculation - max of processor availability and max predecessor time
                    est[task_idx][u] =  max_predecessor_time;
                  
                    
                    // Calculate EFT
                    eft[task_idx][u] = est[task_idx][u] + exec_time[u][task_idx];//check
                    
                    // Check if this is the best processor so far
                    if (eft[task_idx][u] < min_eft) begin
                        min_eft = eft[task_idx][u];
                        best_processor = u;
                    end
                end
                processor_assignment[3*task_idx+:3] = best_processor;
                AFT[task_idx] = min_eft;
                 
            end
                
                    
                
            for (task_idx = NUM_TASKS; task_idx > 0; task_idx = task_idx - 1) begin
//                for (u=0; u< NUM_PROCESSORS ; u=u+1)begin
                    proc = processor_assignment[NUM_PROCESSORS*(task_idx-1)+:NUM_PROCESSORS];
                    AST[task_idx-1] = AFT[task_idx-1] - exec_time[proc][task_idx-1];
//                end
            end    
                
                
            for (task_idx = NUM_TASKS; task_idx > 0; task_idx = task_idx - 1) begin
                for (u=0; u< NUM_PROCESSORS ; u=u+1)begin
                    if (task_idx-1 == NUM_TASKS-1) begin
                            // Entry task
                            lft[task_idx-1][u] = D_G;
                    end 
                    else begin
                        min_successor_time = 32'hFFFFFFFF;

                        // Find maximum of (AFT[predecessor] - comm_cost[task_idx][successor]) among all predecessors
                            for (j = 0; j < NUM_TASKS; j = j + 1) begin
                                if (comm_cost[task_idx-1][j] > 0) begin // If i is a successor of task_idx
                                // Add communication cost only if tasks run on different processors
                                    if (processor_assignment[NUM_PROCESSORS*j+: NUM_PROCESSORS] != u) begin
                                        if (AST[j] - comm_cost[task_idx-1][j] < min_successor_time) begin
                                            min_successor_time = AST[j] - comm_cost[task_idx-1][j];
                                        end 
                                    end
                                    else begin
                                        if (AST[j] < min_successor_time) begin
                                            min_successor_time = AST[j];
                                        end 
                                    end
                                end
                            end
                            lft[task_idx-1][u] =  min_successor_time;
                    end
                    
                
                end
                        
                
                
                
            end//end of task_idx for
            
            for (task_idx = 0; task_idx < NUM_TASKS; task_idx = task_idx + 1) begin
                for (u=0; u< NUM_PROCESSORS ; u=u+1)begin
                    if (lft[task_idx][u] > est[task_idx][u]) begin
                        MET[task_idx][u] = lft[task_idx][u] - est[task_idx][u];
                    end
                    else begin
                        MET[task_idx][u] = 0;
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
                best_processor = -1;//changed
                E_min = 32'hFFFFFFFF;
                for (u=0; u< NUM_PROCESSORS ; u=u+1)begin
                    if(MET[current_task][u] >= exec_time[u][current_task])begin
                        UBET[current_task][u] = fp_div(fp_mult(freq[0], int_to_fp(exec_time[u][current_task])),freq[L-1]);
                        MET_low[current_task][u] = (int_to_fp(MET[current_task][u])<UBET[current_task][u])?int_to_fp(MET[current_task][u]):UBET[current_task][u];
                        f_low[current_task][u]= fp_div(fp_mult(freq[0], int_to_fp(exec_time[u][current_task])),MET_low[current_task][u]);                        
                        P_sta[current_task] = fp_mult(int_to_fp(Lg), fp_mult(volt[0], int_to_fp(I_sub)) + fp_mult(int_to_fp(V_bs), int_to_fp(I_j)));
                        //Ceff in order of 10^-8 and freq in MHZ -> /100 (*0.01)
                        P_dyn[current_task] = fp_mult(int_to_fp(C_eff), fp_mult(volt[0], fp_mult(volt[0], nearest_freq(f_low[current_task][u]))))*0.01;
                        
                        E_low[current_task][u] = fp_mult((P_sta[current_task] + P_dyn[current_task]),MET_low[current_task][u])*0.01;//for time units being in order of 10^-2
//                        best_processor = u;
                    end
                    if (E_low[current_task][u] < E_min) begin
                            E_min = E_low[current_task][u];
                            best_processor = u;
                        end
                end
                E_total = E_total + E_low[current_task][best_processor];
                
                frequency_assignment[DATA_WIDTH*current_task+:DATA_WIDTH] = nearest_freq(f_low[current_task][best_processor]);
                processor_assignment[3*current_task+:3] = best_processor;
                AFT[current_task] = lft[current_task][best_processor];
                AST[current_task] = lft[current_task][best_processor]-MET_low[current_task][best_processor];
           end
        end
    endtask
    
    initial begin
        // Initialize processor assignments
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                processor_assignment[NUM_PROCESSORS*i+: NUM_PROCESSORS] = -1;
                E[i] = 0;
            end
    end
    
    initial begin
        // Initialize processor assignments
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                frequency_assignment[DATA_WIDTH*i+: DATA_WIDTH] = 0;
            end
    end
    
    initial begin
        // Initialize processor assignments
            for (i = 0; i < NUM_TASKS; i = i + 1) begin
                energy_assignment[DATA_WIDTH*i+: DATA_WIDTH] = 0;
            end
    end    
    
    // Main state machine
    always @(posedge clk or posedge reset) begin
        if (reset) begin
            E_total = 0;
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
                    state <= CALC_PARAM;
                end
                
                CALC_PARAM: begin
                    calculate_est();
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