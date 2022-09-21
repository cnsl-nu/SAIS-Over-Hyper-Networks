N = 300; %number of nodes
W=15;
H=5;
M=N/W+N/H;
Num_iter = 100; %number of simulations
endtime=15;
tau=0.2;
kappa=0.2;
gamma=1;
c=[2,5,10,15];
diff_dt=[0.01];
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
            Adj=randperm(N);
            Adj_H=reshape(Adj,N/H,H); %Home edge matrix
            Adj=randperm(N);
            Adj_W=reshape(Adj,N/W,W); %Work edge matrix
            Nodes = zeros(1,N); % Creation of the Identity matrix for agent s state,
            infected=randsample(N,3); % pick random infecteed node
            Nodes(infected)=1;
            Adj_all=zeros(M,N); %Adj matrix
            for home=1:N/H
                Adj_all(home,Adj_H(home,:))=1;
            end
            m=0;
            for work=N/H+1:N/H+N/W
                m=m+1;
                Adj_all(work,Adj_W(m,:))=1;
            end
            for t=1:time
                old_Nodes=Nodes;
                Result(iter,t,difc)=length(find(Nodes==1));
                %number of infected nodes
                Result2(iter,t,difc)=length(find(Nodes==2));
                %number of infected nodes
                %             if Result(iter,t,difc)==0 %if there is no infected nodes
                %                 fprintf("Break")
                %                 breakcount=breakcount+1
                %                 break
                %             end
                if t==time
                    break
                end
                fx=sum((old_Nodes==1).*Adj_all,2);
                fx_ale=fx;
                fx(fx>c(difc))=c(difc);
                for node=1:N
                    if old_Nodes(node)==0
                        inf_rate=tau*(fx'*Adj_all(:,node));
                        alert_rate=kappa*(fx_ale'*Adj_all(:,node));
                        if rand()<(1-exp(-inf_rate*dt))
                            Nodes(node)=1;
                        elseif rand()<(1-exp(-alert_rate*dt))
                            Nodes(node)=2;
                        end
                    elseif old_Nodes(node)==2
                        inf_rate=tau*(fx'*Adj_all(:,node));
                        if rand()<(1-exp(-inf_rate/2*dt))
                            Nodes(node)=1;
                        end
                    else
                        if rand()<(1-exp(-gamma*dt))
                            Nodes(node)=0;
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
    title({['Density of infected nodes in SAIS'],['Nodes =',num2str(N),' Hyperedges =',num2str(size(Adj_all,1)),...
        ' W =',num2str(W),' H =',num2str(H)],...
        [' tau = ',num2str(tau),'; kappa =',num2str(kappa),', gamma =',num2str(gamma),', dt =',num2str(dt),'.'],...
        ['No infection spreading = ', num2str(breakcount), '%.'], ...
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
    title({['Density of alerted nodes in SAIS'],['Nodes =',num2str(N),' Hyperedges =',num2str(size(Adj_all,1)),...
        ' W =',num2str(W),' H =',num2str(H)],...
        [' tau = ',num2str(tau),'; kappa =',num2str(kappa),', gamma =',num2str(gamma),', dt =',num2str(dt),'.'],...
        ['No infection spreading = ', num2str(breakcount), '%.'], ...
        ['Number of simulations: ',num2str(Num_iter)]});
    xlabel('time')
    ylabel('Density of alerted nodes')
    legend('Location','southeast')
    hold off
end