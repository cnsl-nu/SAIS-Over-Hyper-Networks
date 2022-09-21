M = readmatrix("senate-committees\hyperedges-senate-committees.txt");%loading the graph parameters
M2 = readmatrix("senate-committees\node-labels-senate-committees.txt");%loading the graph parameters
N = size(M2,1); %number of nodes
Num_iter = 100; %number of simulations
endtime=15;
tau=0.2;
kappa=0.2;
gamma=1;
c=[2,5,10,15];
diff_dt=[0.01,0.001];
Adj_all=zeros(size(M,1),size(M2,1)); %Adj matrix
for i=1:size(M,1)
    V = rmmissing(M(i,:));
    Adj_all(i,V)=1;
end
for dt_num=1:length(diff_dt)
    dt=diff_dt(dt_num);
    time=endtime/dt;
    Result=zeros(Num_iter,time,length(c));
    Result2=zeros(Num_iter,time,length(c));
    breakcount=0;
    for difc = 1:length(c)
        for iter = 1:Num_iter%% correct type of edges below
            dt
            difc
            iter
            Nodes = zeros(1,N); % Creation of the Identity matrix for agent s state,
            infected=randsample(N,3); % pick random infecteed nodes
            Nodes(infected)=1;
            for t=1:time
                old_Nodes=Nodes;
                Result(iter,t,difc)=length(find(Nodes==1));
                %number of infected nodes
                Result2(iter,t,difc)=length(find(Nodes==2));
                %number of alerted nodes
                %             if Result(iter,t,difc)==0 %if there is no infected nodes
                %                 fprintf("Break")
                %                 breakcount=breakcount+1
                %                 break
                %             end
                fx=sum((old_Nodes==1).*Adj_all,2);
                fx_ale=fx;
                fx(fx>c(difc))=c(difc);
                for node=1:N
                    if old_Nodes(node)==0 %if node sus
                        inf_rate=tau*(fx'*Adj_all(:,node));
                        alert_rate=kappa*(fx_ale'*Adj_all(:,node));
                        if rand()<(1-exp(-inf_rate*dt))
                            Nodes(node)=1; %sus=>inf
                        elseif rand()<(1-exp(-alert_rate*dt))
                            Nodes(node)=2; %sus=>alert
                        end
                    elseif old_Nodes(node)==2 %if node alert
                        inf_rate=tau*(fx'*Adj_all(:,node));
                        if rand()<(1-exp(-inf_rate/2*dt))
                            Nodes(node)=1; %alert=>inf
                        end
                    else %if node inf
                        if rand()<(1-exp(-gamma*dt))
                            Nodes(node)=0; %inf=>sus
                        end
                    end
                end
            end
        end
    end
    breakcount=breakcount/Num_iter/length(c)*100
    %% I nodes
    figure
    Res_av=sum(Result,1)/Num_iter;
    hold on
    t=linspace(1,time,length(Res_av(:,:,1)));
    for i=1:length(c)
        plot(t,Res_av(:,:,i)','DisplayName',['c=',num2str(c(i))])
    end
    title({['Density of infected nodes in SAIS'],['Nodes =',num2str(N),'Hyperedges =',num2str(size(Adj_all,1))],[' tau = ',num2str(tau), ...
        '; kappa =',num2str(kappa),', gamma =',num2str(gamma),', dt =',num2str(dt),'.'],['No infection spreading = ', num2str(breakcount), '%.'], ...
        ['Number of simulations: ',num2str(Num_iter)]});
    xlabel('time')
    ylabel('Density of infected nodes')
    legend('Location','southeast')
    hold off
    %% A nodes
    figure
    Res_av=sum(Result2,1)/Num_iter;
    hold on
    t=linspace(1,time,length(Res_av(:,:,1)));
    for i=1:length(c)
        plot(t,Res_av(:,:,i)','DisplayName',['c=',num2str(c(i))])
    end
    title({['Density of alerted nodes in SAIS'],['Nodes =',num2str(N),'Hyperedges =',num2str(size(Adj_all,1))],[' tau = ',num2str(tau), ...
        '; kappa =',num2str(kappa), ', gamma =',num2str(gamma),', dt =',num2str(dt),'.'],['No infection spreading = ', num2str(breakcount), '%.'], ...
        ['Number of simulations: ',num2str(Num_iter)]});
    xlabel('time')
    ylabel('Density of alerted nodes')
    legend('Location','southeast')
    hold off
end
