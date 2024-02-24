function distances = evaluate_functional_match(responses1, responses2)

distances = zeros(size(responses1,1),1);

for i=1:size(responses1,1)
    %distances(i) = sqrt(sum((responses1(i,:) - responses2(i,:)).^2));
    distances(i) = 1-corr(responses1(i,:)',responses2(i,:)');
end

[sorted,idx] = sort(distances);

figure();
for i=1:length(idx)
    if(i<=8*8)
    subplot(8,8,i);
    plot(responses1(idx(i),:)); hold on;
    plot(responses2(idx(i),:));
    
    title(strcat(num2str(idx(i)),' : ', num2str(round(sorted(i),4))));
    end
end
