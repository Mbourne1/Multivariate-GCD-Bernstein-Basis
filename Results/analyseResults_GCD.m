function [] = analyseResults_GCD(ex_num)





T = readtable("Results_o_gcd.dat");



% Filter based on example number
filteredTable = T(T.EX_NUM == ex_num, : )

end