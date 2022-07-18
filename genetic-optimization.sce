// This file is released under the MIT license
//
// sharks and sardins: genetic-optimization
//
// Feel free to use this code for your combinatoric problems
// Please, refer the author : Israël Tankam Chedjou (https://github.com/israeltankam)

    clear // Clears variables
    clc // Cleans the console

//******************************** Input parameters for the optimization *************************************************
      Tmax = 4000; P0 = 100 ; tau = 37;

//----------------------- Epidemiological parameters ----------------------//
      var_beta = 0.1; //  var_beta because beta is a reserved name
      var_gamma = 0.5; // var_gamma because gamma is a reserved name 
      var_d = 210; // var_d because another variable (properly d) will be modified over the seasons
      rho = 0.025; K = 150; a = 0.0002; alpha = 400; delta = 60; mu = 0.04; Omega = 0.0495;
      var_tf = 330; // var_tf because another variable (properly tf) will be modified over the seasons.
                    // The variable tf is named capital D in the paper
      q = 0.05; r = 1/3;
      var_S0 = 60; // var_S0 because another variable (properly S0) will be modified over the seasons
      var_X0 = 0; // Same comment
      m = 0.3; 
      sucker_cost = 1800; // Cost of a healthy sucker
      t0 = 0;
      step = 0.1; // The numerical integration step. Scilab requires it unlike Matlab.

//****************************** Parameters for the genetic algorithm *****************************************************
      global n_profiles; global disc_profiles; global alpha_mutation;

      alpha_mutation = 1 // initial mutation rate
      q_mutation = 0.99 // Rate of mutation rate decrease
      n_profiles = 10 // Number of profiles in a population for the genetic algorithm
      disc_profiles = 4 // Number of discriminated profiles, i.e. unselected profiles for the next generation
      smallest_mutation_rate = 0.5 // The smallest mutation rate as stopping criterium


//************************************************** Useful functions ******************************************************

//**************************************************************************************************************************
//                                        Functions for the dynamics and cost computation                                  *
//**************************************************************************************************************************

///////////////////////////////////////// Model equations  /////////////////////////////////////////
      // The state varaibles are in the P,S,X order
      function yprim= fgrowing(t,y) // Root growth period subsystem
        yprim=[-var_beta*y(1)*y(2) + alpha*a*(1-var_gamma)*y(2)*y(3)/(y(2)+delta) - Omega*y(1) //dP/dt 
                rho*y(2)*(1-(y(2)/K))-a*y(2)*y(3)/(y(2)+delta) //dS/dt
                var_beta*y(1)*y(2) + alpha*a*var_gamma*y(2)*y(3)/(y(2)+delta) - mu*y(3)//dX/dt
                0 // This state variable is the one that will calculate the yield so that by integrating this component of the diff equa...
              ] //... over the period of growth of the fruit we have the direct expression of the yield of a season.
      endfunction
      function yprim= fflow(t,y) // Fruit growth period subsystem
        yprim=[-var_beta*y(1)*y(2) + alpha*a*(1-var_gamma)*y(2)*y(3)/(y(2)+delta) - Omega*y(1)
               -a*y(2)*y(3)/(y(2)+delta)
                var_beta*y(1)*y(2) + alpha*a*var_gamma*y(2)*y(3)/(y(2)+delta) - mu*y(3)
                m*y(2); // And here is this additional component again; related to the yield parameter m.
              ]
      endfunction
      function yprim = fsurvival(t,y) // Pest fallow survival equation
        yprim=[-Omega*y(1)
               0
               0
               0
              ]
      endfunction    
       /////////////////////////////////////////////////////////////////////////////////////////////////////
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Here we are going to write a function that takes as input any distribution of chain lengths and returns           //
        // the banana crop profit for that chain distribution. Between each chain, the same fallow period...                 //
        // ...\tau^* (tau parameter at the top) will be deployed.                                                            //
        // The input parameter is a vector that gives the distribution of chain lengths. Example [4 5 2 3] means...          //
        // ...we start with a chain of 4 cropping seasons then we deploy a fallow, we continue with a chain of ...           //
        // ...5 cropping seasons, fallow, 2 cropping seasons, fallow, 3 cropping seasons. //
       ///////////////////////////////////////////////// ///////////////////////////////////////////////// ////////////////////
       function computed_profit = profit(dist)
           number_of_chains = length(dist)
           t_tot =[]; // Will compute the total duration of all seasons in all chains plus all fallows.
           solution =[]; // The model solution over the total time span
           d = var_d;
           tf = var_tf;
           Income0 = - sucker_cost; // When we plant we begin with a negative income (the expenditure on the healthy sucker)
           S0 = var_S0;
           X0 = var_X0;
           t0 = t0+step // We begin at t_0^+
           for i=1:number_of_chains // We go through all the chains
               number_of_seasons = dist(i) // We retrieve the length of each chain
               for j=1:number_of_seasons
                    t = t0:step:d; // We discretize the time on the root growing period
                    t_tot = [t_tot t]; // We add to the total time
                    Sol1 = ode([P0;S0;X0;Income0],t0,t,fgrowing); // We solve the diff equation over the root growth period
                                                               // The fourth state variable, initialized in Income0, will compute the profit
                    t = d+step:step:tf; // We discretize the time on the fruit growing period
                    t_tot = [t_tot t]; // We add to the total time
                    Sol2 = ode([Sol1(1,:)(length(Sol1(1,:)));Sol1(2,:)(length(Sol1(2,:)));Sol1(3,:)(length(Sol1(3,:)));Sol1(4,:)(length(Sol1(4,:)))],d+step,t,fflow);
                    // The initial conditions on this subsystem are the endpoints of the solutions on the first subsystem    
                    solution = [solution Sol1 Sol2]; // We concatenate to the solution over the total time
                    if(j<>number_of_seasons)then // If we are not in the last season of the chain, we continue with natural growth.
                        if ~exists("Sol2","local") then
                        disp(dist);
                        continue;
                        end
                        Pfin = Sol2(1,:);Sfin=Sol2(2,:);Xfin=Sol2(3,:);Income = Sol2(4,:)// We retrieve the solutions
                        P0 = (Pfin(length(Sol2(1,:))) + (1-r)*Xfin(length(Sol2(1,:))))// Switch equation between two consecutive seasons
                                                                                      // The nematodes in the soil are those that were already there plus
                                                                                      // the proportion inherited from the root of the mother plant
                        S0 = r*Sfin(length(Sol2(1,:)));
                        X0 = r*Xfin(length(Sol2(1,:)));
                        Income0 = Income(length(Income));// Income accumulates
                                                // The final value of this variable therefore computes the integral of m*S(t) over the duration of fruit growth over all seasons
                        t0 = t0 + var_tf + step; // The initial time is reset to that of the following season
                        d = var_d + t0 - step;
                        tf = var_tf + t0 - step;
                     end
                end
                if ~exists("Sol2","local") then
                        disp(dist);
                        disp("***");
                        disp(number_of_seasons)
                        continue;
                end
               // And here we come out of the last season of the chain. There will therefore be a transition to fallow, and the switch that goes with it
                Pfin = Sol2(1,:);Xfin=Sol2(3,:);Income = Sol2(4,:); // We retrieve the solutions
                Pinit = Pfin(length(Sol2(1,:)))+q*Xfin(length(Sol2(3,:))); // Only free-living nematodes in the soil survive, plus a proportion q of the infestants
                if tau < step // One makes negligible the values of tau which would be lower than the step of integration; this solves numerical problems
                  Sol3 =[];
                  PP = Pinit;
                else
                  t=tf+step:step:tf+tau;
                  t_tot = [t_tot t];
                  Sol3 = ode([Pinit;0;0;Income(length(Income))],tf+step,t,fsurvival); // We solve the hostless survival equation
                                                                                // Income accumulates
                  PP = Sol3(1,:);
                end
                solution = [solution Sol3]; // We concatenate the solution
                P0 = PP(length(PP)); // The nematodes in the soil are inherited from the previous chain
                S0 = var_S0 // The healthy sucker has its biomass
                X0 = var_X0 // And it is healthy in case var_X0=0
                Income0 = Income(length(Income)) - sucker_cost; // The cost of sucker is subtracted from the cumulative income
                t0 = t0 + var_tf + tau + step; // We put the time at the beginning of the new chain
                d = var_d + t0 - step;
                tf = var_tf + t0 - step;
           end
           computed_profit = solution(4,:)(length(solution(4,:)));
           // Profit is the final value of the income that has accumulated throughout the chains
       end






//**************************************************************************************************************************
//                                             Functions for thegenetic algorithm                                          *
//**************************************************************************************************************************
        function y = rand_fix_sum(x) // Uniformely rebalances a summing-to-N integer vector after that arises from rounding  
                                     // a real vector with entries summing to N 
           int_y = round(x)
           d = sum(int_y)-sum(x)
           index = 1+round((length(x)-1)*rand(1))
           while int_y(index)<=d && d>0
             index = round(length(x)*rand(1))
           end
           int_y(index)=int_y(index)-d 
           y = int_y
        endfunction
        
        function x = randfixedsum(n,s,a,b) // Draws a random real vector of length n with entries taken between a and b
                                           // and summing up to s
           s = (s-n*a)/(b-a);
           k = max(min(floor(s),n-1),0);
           s = max(min(s,k+1),k);
           s1 = s - [k:-1:k-n+1]; 
           s2 = [k+n:-1:k+1] - s;
           w = zeros(n,n+1); w(1,2) = b;
           t = zeros(n-1,n);
           tiny = 2^(-1074); //The smallest positive matlab 'double' no.
           for i = 2:n
               tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
               tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
               w(i,2:i+1) = tmp1 + tmp2;
               tmp3 = w(i,2:i+1) + tiny; 
               tmp4 = (s2(n-i+1:n) > s1(1:i)); 
               t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
           end
           x = zeros(n,1);
           rt = rand(n-1,1);
           rs = rand(n-1,1);
           s = repmat(s,1,1);
           j = repmat(k+1,1,1);
           sm = 0; pr = 1; 
           for i = n-1:-1:1  
               e = (rt(n-i,:)<=t(i,j)); 
               sx = rs(n-i,:).^(1/i); 
               sm = sm + (1-sx).*pr.*s/(i+1);
               pr = sx.*pr;
               x(n-i) = sm + pr.*e;
               s = s - e; j = j - e;
           end
           x(n) = sm + pr.*s;
           x=(b-a)*x+a
           x = x(grand(1, "prm", (1:n)')')
           x = rand_fix_sum(x)
      endfunction

//*************************************************************************************************************

      function prfs = create_profiles(n,s) //Creates n_profiles (global variable) different profiles of size n with sum s
           prfs = zeros(n,n_profiles)
           for i=1:n_profiles
               prfs(:,i)=randfixedsum(n,s,1,s)
               if min(prfs(:,i))==0 then
               disp(prfs(:,i))
               end
           end
            // Use prfs(:,i) for single profiles
      endfunction

      function tab = table_of_costs(profiles)
           tab = zeros(1,size(profiles,2))
           for i=1:length(tab)
             tab(i) = profit(profiles(:,i))
           end
      endfunction

      function prof = sort_profiles(profiles) // Sorts the matrix of profiles from the highest income to the lowest
           prof = zeros(profiles)
           tab_costs = table_of_costs(profiles)
           [sorted,J]=gsort(tab_costs,'g','d')
           for i=1:length(J)
               prof(:,i)=profiles(:,J(i))
           end
      endfunction

      function prof = discriminate(profiles,j) // Slices the profiles table to it first n_profiles-disc_profiles profiles (j is disc_profiles). 
                                               //To be applied on sorted profiles
           s = size(profiles)
           keep = s(2)-j
           prof = profiles(:,1:keep)
      endfunction

      function v = mutate_vector(v) // This random mutation dicrements the size of a random chain longer than 1
                              // and increments the size of a random other chain in counterpart
           if length(v)>=2 then
               w = find(v>1) // a vector with the index of all entries of v that are greater than 1
               if length(w) >=1 then
                  r = 1+round((length(w)-1)*rand(1)) // choose a random index in w
                  pos_in_v = w(r)
                  pos_to_switch = 1+round((length(v)-1)*rand(1)) // choose a random position in v to alternate with
                  while pos_to_switch == pos_in_v //make sure it is not the same position
                     pos_to_switch = 1+round((length(v)-1)*rand(1)) 
                  end
                  v(pos_in_v) = v(pos_in_v) - 1
                  v(pos_to_switch) = v(pos_to_switch) + 1
               end
           end
      endfunction

      function mutants = create_mutants(prof,j,al) // Creates j (global disc_profiles in the call) mutants of kept profiles prof (dim n, n_profiles-disc_profiles) 
                                                   // With mutation rate al (global alpha_mutation)
                                                   // Apply on sorted profiles
           s = size(prof,2)
           l = size(prof,1)
           su = sum(prof(:,1))
           mutants = zeros(l,j)
           if s >= j then
                 choosen = samwr(j,1,1:s) // Create a vector of size j drawn among column indexes of prof
                 for i=1:j
                      thres = rand(1)
                      if thres <= al then
                           mutants(:,i)=mutate_vector(prof(:,choosen(i))) // Mutations occur at rate al%
                      else
                           mutants(:,i)=prof(:,choosen(i))
                      end
                 end
           else // j is greater than s = m-j
                 kk = pmodulo(j,s)
                 for i=1:j-kk 
                      thres = rand(1)
                      if thres <= al then
                            mutants(:,i)=mutate_vector(prof(:,i)) //Mutations occur at rate al%
                      else
                            mutants(:,i)=prof(:,i)
                      end
                 end
                 if kk>0
                       indexes = samwr(kk,1,1:s)
                       for i=j-kk+1:j
                            thres = rand(1)
                            if thres <= al then
                                 mutants(:,i)=mutate_vector(prof(:,indexes(i-j+kk))) //Mutations occur at rate al%
                            else
                                 mutants(:,i)=prof(:,indexes(i-j+kk))
                            end
                       end
                 end
          end
     endfunction

    function prof = profile_next_gen(profiles,j,al) //Create population of next generation of profiles, with j=number of discriminated profiles and al = mutation rate
         sorted_prof = sort_profiles(profiles)
         selected = discriminate(sorted_prof,j)
         mutants = create_mutants(selected,j,al)
         prof = [selected mutants]
    endfunction

    function best = best_profile(profiles) // Gives the highest income profile in population "profiles"
         sorted_prof = sort_profiles(profiles)
         best = sorted_prof(:,1)
    endfunction


//******************************************** Implementation ******************************************** //

pmax = floor((Tmax-var_tf)/(tau+var_tf))
top = [floor(Tmax/var_tf)] // The longer chain, when p = 0
for p=1:pmax
    Nmax = floor((Tmax-p*tau)/var_tf)
        profiles = create_profiles(p+1,Nmax)
        disp(best_profile(profiles))
        disp(profit(best_profile(profiles)))
        stagnation_index = 0;
        i = 0;
        while (stagnation_index<=100 && alpha_mutation >= smallest_mutation_rate)
            i=i+1
            cbp = profit(best_profile(profiles))
            profiles = profile_next_gen(profiles,disc_profiles,alpha_mutation)
            if pmodulo(i,100) == 0 then
                alpha_mutation = q_mutation*alpha_mutation
            end
            if profit(best_profile(profiles)) == cbp then
            stagnation_index = stagnation_index + 1
            end 
        end
        if stagnation_index > 100 then
            disp("Stagnation")
        end
        disp("alpha=",alpha_mutation)
        best = best_profile(profiles)
        disp("The best profile with ",p," fallow period is",best)
        if profit(best)> profit(top)then
            top = best
        end
end
disp("The best is ",top," and its profit is ", profit(top))
