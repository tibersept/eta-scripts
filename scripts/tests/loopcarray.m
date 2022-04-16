% Example of a loop on a cell array
function ca = loopcarray()

  indices = [101 102 103 104 105 106];
  
  for i=1:length(indices)
    carray{i} = [];
  end
  
  for i=1:length(indices)
    for j=1:100
      carray{i} = [carray{i} randi(10,20,1)];
    end
  end
  
  ca = carray;
end