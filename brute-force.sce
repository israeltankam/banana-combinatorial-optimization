    // This file is released under the MIT license
    //
    // sharks and sardins: brute-force
    //
    // Feel free to use this code for your combinatoric problems
    // Please, refer the author : Israël Tankam Chedjou (https://github.com/israeltankam)
    
    ///////////////////////// Moldel parameters ////////////////////////////////
    //----------------- Input parameters for the optimization ----------------//
    text_param=['Time horizon  T_{max}='; '$Initial soil infestation P_0=';'Constant fallow duration tau*='];
    param=x_mdialog('Input parameters of the optimization',text_param,['4000';'100';'37']);
    Tmax = evstr(param(1)); 
    P0 = evstr(param(2));
    tau = evstr(param(3));
    
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
       //////////////////////////////       copy-paste of the all-compositions functions      ///////////////////////////////////////
       //      The all-composition function will compute all the distribution of chain sizes summing up to a desired maximum       //
       function [poss,tab] = get_first_composition(n,k)
        if n<k then
            poss = %F
        end
        for i=1:1:k
            tab(i)=1
        end
        tab(k) = n-k+1
        poss = %T
        tab = tab'
       endfunction

       function [poss,tab] = next_compostisions(n,k,input_tab)
        tab = input_tab
        if tab(1) == n-k+1 then
            poss = %F
        end
        last = k
        while (tab(last) == 1)
            last = last-1
        end
        z = tab(last)
        if(last-1 >= 1)then
            tab(last-1) = tab(last-1) + 1
            tab(last) = 1
            tab(k) = z - 1
            poss = %T
        end
       endfunction
       function tab_compositions = all_compositions(n,k)
        if n==k then
             tab_compositions = ones(1,k)
        else
             tab_compositions = []
             [poss, tab] = get_first_composition(n,k)
             while poss
                tab_compositions = [tab_compositions;tab]
                [poss, tab] = next_compostisions(n,k,tab)
             end
        end
       endfunction
       /////////////////////////////////////     Here we go with the brute force resolution    ///////////////////////////////
       clc
       pmax = floor((Tmax-var_tf)/(tau + var_tf)) // var_tf is D
       optimal = []
       profit_max = -1e30
       for p=0:pmax
           Nmax = floor((Tmax-p*tau)/var_tf)
           matrix_compositions = all_compositions(Nmax,p+1)
           number_of_compositions = size(matrix_compositions,1)
           for i=1:number_of_compositions
               distribution = matrix_compositions(i,:)
               profit_temp = profit(distribution)
               if profit_temp > profit_max then
                   optimal = distribution
                   profit_max = profit_temp
               end
           end           
       end
       disp(optimal)
       disp(profit_max)
