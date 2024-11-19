function mean_A=getMeanOfNonNaN(A)
mean_A=zeros(size(A,1),1);
for i=1:size(A,1)
nonnan_elements_index = find(~isnan(A(i,:)));
nonnan_elements = A(i,nonnan_elements_index);
mean_A(i) = mean(nonnan_elements);
end
end