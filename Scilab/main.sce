// code of zain ul abideen
clear; clc;

disp("ICARUS v1.0")

inpthchk = 0
outpthchkm = 0
outpthchkr = 0
tgtchk = 0

while 1==1
    usrinp = input(">> ", "string")
    
    if usrinp == "halt"
        abort
    
    elseif usrinp == "help"
        disp("ICARUS is an s- process Nucleosynthesis code.")
        disp("it was developed by Zain Ul Abideen (BS FYP).")
        disp("it takes (n,g) cross-sections as inputs, and ")
        disp(" gives MACS and Reaction Rates as outputs.   ")
        disp("")
        disp("following keywords are mandatory for running:")
        disp("set inpath <enter> /home/data.txt")
        disp("set outpathm <enter> /home/macs.txt")
        disp("set outpathr <enter> /home/rates.txt")
        disp("set target <enter> 91.9")
        disp("")
        
    elseif usrinp == "set inpath"
        inpth = input("enter the path to the data file: ", "string")
        inpthchk = 1

    elseif usrinp == "set outpathm"
        outpthm = input("enter the path for the macs file: ", "string")
        outpthchkm = 1
        
    elseif usrinp == "set outpathr"
        outpthr = input("enter the path for the rates file: ", "string")
        outpthchkr = 1
        
    elseif usrinp == "set target"
        tgt = input("enter the mass number of target: ")
        tgtchk = 1
        
    elseif usrinp == "execute"
        if (inpthchk == 1 && outpthchkm == 1 && outpthchkr == 1 && tgtchk == 1)
            break
        else
            disp("all mandatory options need to be set first.")
        end
        
    else
        disp("invalid command. type *help* to get some help.")
    end
end

// sheet = readxls(pth)
// data = sheet(1)

data = fscanfMat(inpth)

eg = round(data(:,1) * 1000)
cs = data(:,2)

stepsize = eg(2) - eg(1)

mn = 939.6
mamu = 931.5

// target = 91.9
target = tgt

for kt = 5 : 5 : 100
    mred = ((mn) .* (mamu) .* (target)) ./ ((mn) + ((mamu) .* (target)))
    vt = sqrt((2 .* (kt ./ 1000) .* (3 .* 10 .^ 8) .^ 2) ./ (mred))
    
    const = (2 ./ sqrt(%pi)) .* (1 ./ (kt) .^ 2)
    
    for i = 1 : length(cs)
        tmp(i) = (cs(i) .* eg(i) .* %e .^ (-eg(i) ./ kt)) .* stepsize
    end
    
    mxw = const .* tmp
    mxw = sum(mxw)
    
    marr(kt ./ 5, 1) = kt
    marr(kt ./ 5, 2) = mxw
    
    redavo = 6.022 .* 10 .^ (-2)
    
    rr = mxw .* vt .* redavo
    
    rrarr(kt ./ 5, 1) = kt
    rrarr(kt ./ 5, 2) = rr
end

write(outpthm, marr)
write(outpthr, rrarr)

disparr(:,1) = marr(:,1)
disparr(:,2) = marr(:,2)
disparr(:,3) = rrarr(:,2)

disp(disparr)

disp("ICARUS congratulates you on a successful calculation.")

// THE END .
