function cp_mixture = getMixtureCp(temperature, massFractions, massFractionNames)
% Calculates the heat capacity of a mixture, assuming an ideal gas

cp_components = zeros(1,length(massFractions));

for i = 1:length(massFractions)
    cp_components(i) = getComponentCp(temperature,massFractionNames(i));
end

cp_mixture = 0;
for i = 1:length(massFractions)
    cp_mixture  = cp_mixture  + cp_components(i)*massFractions(i);
end

end