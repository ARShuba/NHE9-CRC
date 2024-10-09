clear all
[x,v,s,h,f] = NHE_Bifur; 
a = x(9,:); %bifurcation parameter
b = x(2,:); 
c = a./1000;

%% Based on eigenvalues to judge stable vs. unstable states
ind = zeros(1,4);
snum = size(f);
num = snum(2);
j = 1;
n = 1;

for n = 1:1:(num-1)
    x1 = find(f(:,n) > 0);
    x2 = find(f(:,n+1) > 0);
    if isempty(x1) && ~isempty(x2)
        ind(j) = n + 1;
        j = j + 1;
    elseif ~isempty(x1) && isempty(x2)
        ind(j) = n + 1;
        j = j + 1;
    end
end

%%

array1=[c(ind(1)) c(ind(2)) c(ind(3)) c(ind(4))];
figure1 = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

plot(c(1:ind(1)),b(1:ind(1)),'b');
hold on
plot(c(ind(1)+1:ind(2)),b(ind(1)+1:ind(2)),'r');
plot(c(ind(2)+1:ind(3)),b(ind(2)+1: ind(3)),'b');
plot(c(ind(3)+1:ind(4)),b(ind(3)+1:ind(4)),'r');
plot(c(ind(4)+1:end),b(ind(4)+1:end),'b');

xlim([0 700]);
xlabel('S ext (10^3 molecules)');
ylabel('Zeb mRNA (10^3 molecules)');