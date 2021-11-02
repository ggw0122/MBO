function [InitFunction, CostFunction, FeasibleFunction] = Ackley

InitFunction = @AckleyInit;
CostFunction = @AckleyCost;
FeasibleFunction = @AckleyFeasible;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = AckleyInit(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.1;
MinParValue = -30*ones(1,OPTIONS.numVar);
MaxParValue = 30*ones(1,OPTIONS.numVar);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = MinParValue + (MaxParValue - MinParValue + 1) .* rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = AckleyCost(OPTIONS, Population)

% Compute the cost of each member in Population

global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    sum1 = 0;
    sum2 = 0;
    for i = 1 : p
        x = Population(popindex).chrom(i);
        sum1 = sum1 + x^2;
        sum2 = sum2 + cos(2*pi*x);
    end
    Population(popindex).cost = 20 + exp(1) - 20 * exp(-0.2*sqrt(sum1/p)) - exp(sum2/p);    
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = AckleyFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue(k));
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue(k));
    end
end
return;