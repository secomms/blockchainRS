##Script to find optimal parameters for ASBK protocol.
##The script finds the optimal values of s (number of samples per light node) and R' (rate of component RS code), to minimize the amount of
##data each light node downloads.
##Input values:
# - q : finite field size (can be any prime power);
# - block_size : size of uncoded list of transactions (expressed in kB);
# - hash_size : digest length (expressed in bits); 
# - m : number of light nodes connected to each full node;
# - pr_target : target adversarial success probability;
# - s_step : difference between two adiacent tested values of s;
# - num_values_n : number of tested values for n.

##The script considers several values for R' and, for each one of them, computes the minimum s which guarantees that advesaries succeed with probability < pr_target
##For each R', the amount of downloaded data is computed and displayed

###########################################################################################

reset();

#########################################################
##Computes base 2 logarithm

def log2(x):
    return log(x*1.)/log(2.);
#########################################################


#########################################################
##Computes optimal setting for ASBK protocol
##Input values:
# - m : number of light nodes;
# - n_vec : vector with length values of 2D-RS code
# - gamma_ratio_vec : gamma values for 2D-RS code
# - s_step : difference among two adiacent values of s
# - pr_target : target adversarial success probability
# - q : finite field size
# - hash_size : digest length (expressed in bits) 
# - k_prime : dimension of component RS code 

#The function returns two vectors, s_values and data_size, having the same length as n_vec. The i-th entry of s_values is the minimum value of s which
#guarantees an adversarial success probability of pr_target, when the 2D-RS code has length n_vec[i]; the i-th entry of data_size is the corresponding
#amount of downloaded data

def find_optimal_setting(m,n_vec,gamma_ratio_vec,s_step,pr_target,q,hash_size,k_prime):   
   
    data_size = [];
    s_values = [];
   
    for n_index in range(len(n_vec)):
      
        n = n_vec[n_index]; #2D-RS code length
        gamma_times_n = round(n*gamma_ratio_vec[n_index]); #minimum amount of symbols which guarantees decoding
       
       
        #Find minimum s such that average number of recovered symbols is at least gamma_ratio*n
        s_min = ceil( n * (1.- (1.-gamma_times_n/n)^(1./m)) );        
        s = s_min;
                
        
        flag_s = 0; #flas_s becomes 1 when we have found the minimum s that yields the desired probability
        while (flag_s==0):                        
            
            x_star = round(n * (1 - (1-s/n)^m));
            d = x_star-gamma_times_n+1;
            
            pr_success = 1-(1-binomial(gamma_times_n-1,s)/binomial(x_star,s))^m;
             
            if pr_success < pr_target:
                flag_s = 1;
            else:
                s += s_step;
       
  
        s_values.append(s);
  
        data_size.append(2*sqrt(n)*hash_size+s*(ceil(log2(q))+hash_size*ceil(log2(sqrt(n)))));
   
   
   ##Uncomment the following lines to print data for tikz figure  
   
  #  print("\\addplot[mark=+,red, mark size=2pt,line width=0.7pt]coordinates{");
  #  for i in range(len(n_vec)):
  #      print("("+str(k_prime/sqrt(n_vec[i])*1.)+", "+str(s_values[i])+")");
  #  print("};\\addlegendentry{Values of $s$, $m = "+str(m)+"$};");
   

  #  print("\\addplot[mark=+,blue, mark size=2pt,line width=0.7pt]coordinates{");
  #  for i in range(len(n_vec)):
  #      print("("+str(k_prime/sqrt(n_vec[i])*1.)+", "+str(data_size[i]/8000.)+")");
  #  print("};\\addlegendentry{Amount of download data, $m = "+str(m)+"$};");
    return s_values, data_size;


##############################################################################################


##Specify input parameters

q = 2^256; #finite field size
block_size = 75; #expressed in kB
hash_size = 256; #digest binary length
m = 1000; #number of light nodes
pr_target = 0.01; #target adversarial success probability
s_step = 1; #difference among two adiacent values of s
num_values_n = 150; #number of tested values for n


#Compute RS codes parameters

k_prime = ceil(sqrt(block_size*8000/ceil(log2(q)))); #RS code dimension (value of k' in the paper)

n_prime_vec = range(ceil(k_prime/0.97),round(k_prime/0.05),max(2,ceil( ( round(k_prime/0.05)-ceil(k_prime/0.97) )/num_values_n)) ); #values of n' that are tested by the script

#For fixed k' and several n', compute the resulting parameters of 2D-RS code
gamma_ratio_vec = vector(RR,len(n_prime_vec));  #vector with values of gamma
real_n_vec = vector(ZZ,len(n_prime_vec)); #vector with lengths of 2D-RS code (n = n'^2)

for i in range(len(n_prime_vec)):
    n_prime = n_prime_vec[i];
    real_n_vec[i] = n_prime^2;
    gamma_ratio_vec[i] = (n_prime^2-(n_prime-k_prime+1)^2+1)/(n_prime^2); #Eq. (4) in the paper
   

#Launch program to find minimum values of s
s_values, download_data = find_optimal_setting(m,real_n_vec,gamma_ratio_vec,s_step,pr_target,q,hash_size,k_prime);


##Print results that minimize ell_D
ell_D = min(download_data);
pos = download_data.index(ell_D);
s = s_values[pos];
best_n_prime = n_prime_vec[pos];
   
print("For R' = k'/n' = "+str(k_prime/n_prime_vec[pos]*1.)+", n' = "+str(n_prime_vec[pos]*1.)+", min download data (ell_D) in kB = "+str(ell_D/8000.)+", s = "+str(s));
