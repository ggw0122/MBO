
function Population1 = CombinePopulation(OPTIONS, Population1, Population2)

numButterfly1 = ceil(OPTIONS.partition*OPTIONS.popsize);

for i = 1: OPTIONS.popsize - numButterfly1
    Population1(numButterfly1 + i).chrom = Population2(i).chrom;
    Population1(numButterfly1 + i).cost = Population2(i).cost;
end

