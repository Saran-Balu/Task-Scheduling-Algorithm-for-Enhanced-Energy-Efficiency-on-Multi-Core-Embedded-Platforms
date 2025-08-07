`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 10.04.2025 12:53:54
// Design Name: 
// Module Name: isaecc_algo
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


module isaecc_algo #(
    parameter NUM_TASKS = 10,
    parameter NUM_PROCESSORS = 3,
    parameter L = 8,            // Number of frequency levels
    parameter E_given_G = 37,   // Energy budget
    parameter C_eff = 1,        // Effective capacitance
    parameter MAX_FREQ = 1500, // Representing 100% of maximum frequency
    parameter DATA_WIDTH = 64,  // Width for execution time and comm cost
    parameter RANK_WIDTH = 64,  // Width for rank and time calculations
     parameter FRAC_BITS = 32,
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
    
    output reg done,
    output reg [0:DATA_WIDTH*NUM_TASKS-1] frequency_assignment,
    output reg [0:DATA_WIDTH*NUM_TASKS-1] energy_assignment,
    output reg [0:3*NUM_TASKS-1] processor_assignment, // Stores which processor each task is assigned to
    output reg [DATA_WIDTH-1:0] E_total           // Energy consumption
);

    // Graph structure - Communication costs between tasks
    // Format: comm_cost[predecessor][successor]
    reg [DATA_WIDTH-1:0] comm_cost [0:NUM_TASKS-1][0:NUM_TASKS-1];
    
    // Execution times for each task on each processor
    reg [DATA_WIDTH-1:0] exec_time [0:NUM_PROCESSORS-1][0:NUM_TASKS-1];
    
    // Stores the rank of each task
    reg [DATA_WIDTH-1:0] task_rank [0:NUM_TASKS-1];
    
    // Sorted tasks based on rank (indices of tasks)
    reg [3:0] sorted_tasks [0:NUM_TASKS-1];
    
    // Processor availability time
    reg [DATA_WIDTH-1:0] processor_avail [0:NUM_PROCESSORS-1];
    
    // Actual Finish Time for each task
    reg [DATA_WIDTH-1:0] AFT [0:NUM_TASKS-1];
    
    reg [DATA_WIDTH-1:0] freq [0:L-1];
    reg [DATA_WIDTH-1:0] volt [0:L-1];
    reg [DATA_WIDTH-1:0] E [0:NUM_TASKS-1];
    
    // State machine states
    reg [3:0] state;
    parameter IDLE = 4'd0;
    parameter LOAD_DATA = 4'd1;
    parameter CALC_RANKS = 4'd2;
    parameter SORT_TASKS = 4'd3;
    parameter ENERGY = 4'd4;
    parameter ASSIGN_TASKS = 4'd5;
    parameter DONE = 4'd6;
    
    // Temporary variables for computation
    reg [3:0] current_task;
    //reg [3:0] current_processor;
    reg [DATA_WIDTH-1:0] est, eft, min_eft;
    reg [DATA_WIDTH-1:0] E_min_G = 0;
    reg [DATA_WIDTH-1:0] E_max_G = 0;
    reg [DATA_WIDTH-1:0] E_so_far = 0;
    reg [DATA_WIDTH-1:0] E_after =0;
    reg [DATA_WIDTH-1:0] E_avg_G;
    reg [DATA_WIDTH-1:0] E_ie_G;
    reg [DATA_WIDTH-1:0] minimum;
    reg [DATA_WIDTH-1:0] maximum;
   // reg [RANK_WIDTH-1:0] E_given_G;
    
    reg [3:0] best_processor;
    reg [3:0] task_idx;
    reg signed [4:0] i, j, k, u;
    reg [DATA_WIDTH-1:0] avg_exec_time [0: NUM_TASKS-1];
    reg [RANK_WIDTH-1:0] max_successor_rank;
    reg [DATA_WIDTH-1:0] temp_rank;
    reg [3:0] temp_task;
    reg [DATA_WIDTH-1:0] max_predecessor_time = 0;
    reg [RANK_WIDTH-1:0] latest_task;
    reg [DATA_WIDTH-1:0] latest_finish_time;
    reg [DATA_WIDTH-1:0] freq_mult_ratio;
    reg [DATA_WIDTH-1:0] freq_mult;
    
    
    
    reg [DATA_WIDTH-1:0] E_max [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] E_min [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] E_pro_min [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] E_pro_max [0:NUM_TASKS-1];
    
    reg [DATA_WIDTH-1:0] P_sta [0:L-1];
    reg [DATA_WIDTH-1:0] P_dyn [0:L-1];
     
    reg [DATA_WIDTH-1:0] E_avg [0:NUM_TASKS-1]; 
    reg [DATA_WIDTH-1:0] recip_el [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] el [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] E_wa [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] E_pre [0:NUM_TASKS-1];
    reg [DATA_WIDTH-1:0] E_given [0:NUM_TASKS-1]; 
    
    
    // Fixed-point multiplication function (Q16.16 format)
    function [DATA_WIDTH-1:0] fp_mult;
        input [DATA_WIDTH-1:0] a;
        input [DATA_WIDTH-1:0] b;
        reg [2*DATA_WIDTH-1:0] temp;
        begin
            temp = a*b;
            fp_mult = temp[95:32];
        end
      endfunction
    
    // Fixed-point division function (Q16.16 format)
    function [DATA_WIDTH-1:0] fp_div;
        input [DATA_WIDTH-1:0] a;
        input [DATA_WIDTH-1:0] b;
        begin
            // Scale up numerator before division to maintain precision
            fp_div =(a  / b) << FRAC_BITS;
        end
    endfunction
    
    // Convert integer to fixed-point
    function [DATA_WIDTH-1:0] int_to_fp;
        input [DATA_WIDTH-1:0] int_val;
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
                    exec_time[i][j] = (exec_time_in[DATA_WIDTH*(NUM_TASKS*i+j)+:DATA_WIDTH] <<FRAC_BITS) *0.01;
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
                AFT[i] = 64'hFFFFFFFFFFFFFFFF;
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
                    avg_exec_time[i] = avg_exec_time[i] + ((exec_time[j][i]>>FRAC_BITS)*100);
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
                        avg_exec_time[i] = avg_exec_time[i] + ((exec_time[j][i]>>FRAC_BITS)*100);
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
    
    //Calculate energies
    task calculate_energy;
        begin
            for (i=0; i< NUM_TASKS ; i=i+1)begin
                E_min[i] = 64'hFFFFFFFFFFFFFFFF;
                E_max[i] = 64'h00000000;
                for(u=0; u<NUM_PROCESSORS; u=u+1) begin
                    E_pro_min[u] = fp_mult(int_to_fp(Lg), fp_mult(volt[L-1], int_to_fp(I_sub)) +
                         fp_mult(int_to_fp(V_bs), int_to_fp(I_j))) + fp_mult(fp_mult(C_eff, (fp_mult (volt[L-1], (fp_mult (volt[L-1],freq[L-1]))))),exec_time[u][i])*0.01; 
                    E_pro_max[u] = fp_mult(int_to_fp(Lg), fp_mult(volt[0], int_to_fp(I_sub)) + fp_mult(int_to_fp(V_bs), int_to_fp(I_j)))+
                            fp_mult(fp_mult(C_eff, (fp_mult (volt[0], (fp_mult (volt[0],freq[0]))))),exec_time[u][i])*0.01; 
                    
                    if (E_pro_min[u]<E_min[i]) begin
                        E_min[i] = E_pro_min[u];
                    end
                    if (E_pro_max[u] > E_max[i]) begin
                        E_max[i] = E_pro_max[u];
                    end
                end    
                    
                 E_avg[i] = (E_min[i] + E_max[i])>> 1; 
            
            end
            
            for (i=0; i< NUM_TASKS ; i=i+1) begin
                E_min_G = E_min_G + E_min[i];
                E_max_G = E_max_G + E_max[i];
            end
            
            E_avg_G = (E_min_G + E_max_G)>>1;
            
             
            for (i=0; i< NUM_TASKS ; i=i+1) begin
                E_ie_G = (E_given_G << FRAC_BITS) - E_min_G;
                
                el[i] = fp_div(E_avg[i]<<7,E_avg_G);
                E_wa[i] = fp_mult(E_ie_G, el[i]) + (E_min[i]<<7);
                //recip_el[i] = fp_div(E_avg_G*2**7,E_avg[i]);
                //E_wa[i] = fp_div(E_ie_G*2**14, recip_el[i]) + E_min[i]<<7;
                
                if (E_wa[i] < E_max[i]<<7)
                    E_pre[i] = E_wa[i];
                else
                    E_pre[i] = E_max[i]<<7;
            end
            
        end
    
    endtask
    
    
    // Assign tasks to processors
    task assign_tasks;
        begin
            for (task_idx = 0; task_idx < NUM_TASKS; task_idx = task_idx + 1) begin
                current_task = sorted_tasks[task_idx];
                AFT[current_task] = 64'hFFFFFFFFFFFFFFFF;
                E_after=0;
                E_so_far=0;
                for (k=0; k<task_idx; k=k+1) begin
                    E_so_far = energy_assignment[DATA_WIDTH*sorted_tasks[k]+: DATA_WIDTH] + E_so_far;
                end
                
                for (k=NUM_TASKS-1; k > task_idx; k= k-1)begin
                    E_after = E_pre[sorted_tasks[k]] + E_after;
                end 
                
                if ((E_given_G<<(FRAC_BITS+7)) > (E_so_far + E_after)) begin
                    E_given[current_task] = (E_given_G<<(FRAC_BITS+7)) - E_so_far - E_after;
                end
                
                
                for (u=0; u<NUM_PROCESSORS; u=u+1)begin
                
                    latest_task = NUM_TASKS+1;
                    latest_finish_time = 0;
                    
                    for (k = 0; k < NUM_TASKS; k = k + 1) begin
                        if (processor_assignment[NUM_PROCESSORS*k+: NUM_PROCESSORS] == u) begin
                            if (AFT[k] > latest_finish_time) begin
                                latest_finish_time = AFT[k];
                                latest_task = k;
                            end
                        end
                    end
                    
                    // Update availability time for processor u
                    if (latest_task == NUM_TASKS+1) begin
                        processor_avail[u] = 0;
                    end
                    
                    else begin
                        processor_avail[u] = AFT[latest_task];
                    end
                                
                    for (i=L ; i>0 ; i=i-1) begin
                        
                        P_sta[i-1] = fp_mult(int_to_fp(Lg), fp_mult(volt[i-1], int_to_fp(I_sub)) + fp_mult(int_to_fp(V_bs), int_to_fp(I_j)));
                    //Ceff in order of 10^-8 and freq in MHZ -> /100 (*0.01)
                       P_dyn[i-1] = fp_mult(C_eff, fp_mult(volt[i-1], fp_mult(volt[i-1], freq[i-1])))*0.01;
                        
                      // E[current_task] =  fp_mult (int_to_fp(C_eff), fp_mult (volt[i-1], fp_mult (volt[i-1], fp_mult (freq[0],exec_time [u][current_task]))))*0.01;
                       //E[current_task] = fp_div (fp_mult ((P_sta[i-1]+P_dyn[i-1]), fp_mult(exec_time [u][current_task], freq[0])),freq[i-1])*0.01;
                        
                        freq_mult = freq[0]<<7;
                        freq_mult_ratio = fp_div(freq_mult,freq[i-1]);
                       
                        E[current_task] = fp_mult((P_sta[i-1]+P_dyn[i-1]),fp_mult(exec_time [u][current_task], freq_mult_ratio)) ;
                        if (E[current_task] > E_given[current_task]) begin
                        
                        end
                        else begin    
                            if (current_task == 0) begin
                                // Entry task
                                est = 0;
                            end 
                            else begin
                                max_predecessor_time = 0;

                                // Find maximum of (AFT[predecessor] + comm_cost[predecessor][current_task]) among all predecessors
                                for (j = 0; j < NUM_TASKS; j = j + 1) begin
                                    if (comm_cost[j][current_task] > 0) begin // If i is a predecessor of current_task
                                    // Add communication cost only if tasks run on different processors
                                        if (processor_assignment[3*j+: 3] != u) begin
                                            if (AFT[j] + (comm_cost[j][current_task]<<FRAC_BITS) > max_predecessor_time) begin
                                                max_predecessor_time = AFT[j] + (comm_cost[j][current_task]<< FRAC_BITS);
                                            end
                                        end 
                                        else begin
                                        // No communication cost if on same processor
                                        if (AFT[j] > max_predecessor_time) begin
                                            max_predecessor_time = AFT[j];
                                        end
                                        end
                                    end
                                end
                                est = (processor_avail[u] > max_predecessor_time) ?  processor_avail[u] : max_predecessor_time;
                            end
                    
                        
                        
                        // Final EST calculation - max of processor availability and max predecessor time
                            
                        
                        
                        // Calculate EFT
                            eft = est + ((fp_mult(exec_time[u][current_task], freq_mult_ratio)*100)>>7);//check
                        
                            if (eft < AFT[current_task]) begin
                                processor_assignment[3*current_task+: 3] = u;
                                frequency_assignment[DATA_WIDTH*current_task+: DATA_WIDTH] = freq[i-1];
                                energy_assignment[DATA_WIDTH*current_task+: DATA_WIDTH] = E[current_task];
                                
                                AFT[current_task]=eft;
                            end
                                 
                        end
                        
                    end
                end
                
              
                E_total = E_total+energy_assignment[DATA_WIDTH*current_task+: DATA_WIDTH]; 
               
               
               // for (i=0; i<NUM_TASKS; i=i+1) begin
               //   E_total = E_total + E[i];
               //end
            
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
                    state <= ENERGY;
                end
                
                ENERGY: begin
                    calculate_energy();
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
