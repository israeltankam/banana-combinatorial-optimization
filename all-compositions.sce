// This file is released under the MIT license
//
// sharks and sardins: all-compositions
//
// Feel free to use this code for your combinatoric problems
// Please, refer the author : Israël Tankam Chedjou (https://github.com/israeltankam)
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

//Example
clc
disp(all_compositions(6,3))
