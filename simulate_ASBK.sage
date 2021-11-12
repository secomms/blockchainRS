##Program to simulate the adversarial success probability


reset();

m = 50; #number of light nodes connected to each full node
n = 100; #2D-RS code length
s_min = 2; #minimum number of samples required by each light node
s_max = 20; #maximum number of samples required by each light node
s_step = 2; #difference between two adjacent tested values of s
gamma_val = 0.5; #ratio of symbols required to decode

num_test = 10^4; #maximum number of protocol executions, for each tested value of s
max_err = 100; #maximum number of adversary successes, for each tested value of s

gamma_times_n = round(gamma_val*n);

#Prepare for simulation
Pn = Permutations(range(n));

sim_values = []; #simulated success probabilities
th_values = []; #theoretical success probabilities
s_values = []; #simulated values of s

pr_fail = max_err;
s = s_min;

while (pr_fail>= max_err)&(s<=s_max):    
    
    pr_fail = 0;
    num_tx = 0;
    sim_pr_no_decode = 0;
    s_values.append(s);
    
    while (num_tx<num_test)&(pr_fail<max_err):
        
       
        num_tx+=1; #update number of protocol executions

        #Define a length-n vector with all ones; every asked symbol is then set to 0
        this_symbols = vector(ZZ,n);

        positions_unknown = range(n); #index of positions which are hidden

        list_of_queried_symbols = []; #list containing the set of symbols asked from all the light nodes

        #For each light node, select s distinct positions
        for j in range(m):
            
            #Select s positions at random
            rand_perm = Pn.random_element(); 
            queried_positions = rand_perm[0:s];
            
            list_of_queried_symbols.append(queried_positions);
            
            #Everytime a symbol is requested, change the value of this_symbols to 1
            for ell in queried_positions:
                this_symbols[ell] = 1;

        #Compute the number x of new asked symbols
        x = this_symbols.list().count(1);
        
        #If the number of received symbols is not enough to decode, then the protocol has failed
        if (x<gamma_times_n):
            pr_fail += 1; #the number of new symbols is not large enough to allow for decoding
        else:
            
            #Compute the minimum number d of symbols which the malicious node cannot reveal
            d = x-gamma_times_n+1; #compute the minimum number of symbols which must not be revealed

            #Select d symbols among the asked ones, at random, which will not be given to the light nodes
            new_pos_asked = this_symbols.support();
            pos_to_deny = Set(new_pos_asked[0:d]);

            #See if there is at least a light node which has received all symbols
            flag_all = 0; #flag_all = 1 when adversary succeeds
            light_node_index = 0;
            
            while(flag_all == 0)&(light_node_index<m):
                pos_deny = Set(list_of_queried_symbols[light_node_index]).intersection(pos_to_deny);
                if len(pos_deny)==0:
                    flag_all = 1;
                light_node_index += 1;
            
            
            if flag_all:
                pr_fail = pr_fail+1;
                
        print(str(num_tx)+"   --  "+str(pr_fail/num_tx*1.)+", num failures = "+str(pr_fail));

#############################################################################

    ### Theoretical estimate
    average_val = round(n*(1-(1-s/n)^m));
    if round(average_val)<gamma_times_n:
        th_pr_fail = 1;
    else:
        th_pr_fail =  1 - (1 - binomial(gamma_times_n-1,s)/binomial(round(average_val),s))^m ;

    print("s = "+str(s)+", sim. pr = "+str(pr_fail*1./num_tx)+", th pr = "+str(th_pr_fail*1.));    

    sim_values.append((s,pr_fail/num_tx*1.));
    th_values.append((s,th_pr_fail*1.));

    s+=s_step;
    

#Print results
x = list_plot(sim_values,color = 'red',scale='semilogy',plotjoined=True);
x += list_plot(th_values,color = 'blue',scale='semilogy',plotjoined=True);
x.show();

#Print results for tikz 

print("-------------------------------------------------------------------");    
print(" ");
print(" ");
print("Soundness simulation for ");
print("- m = "+str(m));
print("- gamma = "+str(gamma_val));
print("- n = "+str(n));
print("- max num of test = "+str(num_test))
print("- max num of errors = "+str(max_err))
print(" ");
print(" ");
print("-------------------------------------------------------------------");
print(" ");
print(" ");
print("\\addplot[mark=+,blue, mark size=2pt,line width=0.7pt]coordinates{");    
for i in range(len(s_values)-1):
    print(str(sim_values[i]));
print("};\\addlegendentry{Emp.};");    
print(" ");
print("\\addplot[mark=+,red, mark size=2pt,line width=0.7pt]coordinates{");
for i in range(len(s_values)-1):
    print(str(th_values[i]));
print("};\\addlegendentry{Th.};");


